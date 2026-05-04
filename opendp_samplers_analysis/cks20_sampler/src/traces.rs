/*
    Traces the CKS20 discrete-Laplace sampler
    We run it with three different bernoulli_exp1 implementations (original, taylor, fixed) 
    and dump one CSV per variant 
    Each row records the first outcome observed at every layer of the pipeline plus the final dlap value
*/

use dashu::{
    integer::{IBig, Sign, UBig},
    rational::RBig,
};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::rng::{fill_bytes, Fallible};

const N_ROWS: usize = 1_000_000;

pub struct SampleRow {
    pub unif: Option<u128>,
    pub unif_below: Option<UBig>,
    pub bernoulli: Option<bool>,
    pub bernoulli_rat: Option<bool>,
    pub bernoulli_exp1: Option<bool>,
    pub bernoulli_exp: Option<bool>,
    pub geom_slow: Option<UBig>,
    pub geom_fast: Option<UBig>,
    pub dlap: Option<IBig>,
}

impl SampleRow {
    pub fn new() -> Self {
        Self {
            unif: None,
            unif_below: None,
            bernoulli: None,
            bernoulli_rat: None,
            bernoulli_exp1: None,
            bernoulli_exp: None,
            geom_slow: None,
            geom_fast: None,
            dlap: None,
        }
    }

    pub fn to_csv_line(&self) -> String {
        fn cell<T: Display>(o: &Option<T>) -> String {
            match o {
                Some(v) => v.to_string(),
                None => String::new(),
            }
        }
        format!(
            "{},{},{},{},{},{},{},{},{}",
            cell(&self.unif),
            cell(&self.unif_below),
            cell(&self.bernoulli),
            cell(&self.bernoulli_rat),
            cell(&self.bernoulli_exp1),
            cell(&self.bernoulli_exp),
            cell(&self.geom_slow),
            cell(&self.geom_fast),
            cell(&self.dlap),
        )
    }
}

// Each of these wraps a primitive sampler and stashes the outcome in the trace, but only the first time it's called per row

fn sample_from_uniform_bytes_u128_traced(trace: &mut SampleRow) -> Fallible<u128> {
    let mut buf = [0u8; 16];
    fill_bytes(&mut buf)?;
    let v = u128::from_be_bytes(buf);
    if trace.unif.is_none() {
        trace.unif = Some(v);
    }
    Ok(v)
}

fn sample_uniform_ubig_below_traced(upper: UBig, trace: &mut SampleRow) -> Fallible<UBig> {
    if upper.is_zero() {
        return Err("upper bound must be > 0".into());
    }
    let byte_len = upper.to_be_bytes().len().max(1);
    let max = UBig::from_be_bytes(&vec![u8::MAX; byte_len]);
    let threshold = &max - (&max % &upper);
    let mut buffer = vec![0u8; byte_len];
    Ok(loop {
        fill_bytes(&mut buffer)?;
        let sample = UBig::from_be_bytes(&buffer);
        if sample < threshold {
            let folded = &sample % &upper;
            if trace.unif_below.is_none() {
                trace.unif_below = Some(folded.clone());
            }
            break folded;
        }
    })
}

fn sample_standard_bernoulli_traced(trace: &mut SampleRow) -> Fallible<bool> {
    let mut buffer = [0u8; 1];
    fill_bytes(&mut buffer)?;
    let res = buffer[0] & 1 == 1;
    if trace.bernoulli.is_none() {
        trace.bernoulli = Some(res);
    }
    Ok(res)
}

fn sample_bernoulli_rational_traced(prob: RBig, trace: &mut SampleRow) -> Fallible<bool> {
    let (numer_ib, denom) = prob.into_parts();
    let (sign, numer) = numer_ib.into_parts();
    if sign == Sign::Negative && !numer.is_zero() {
        return Err("numerator must not be negative".into());
    }
    if numer > denom {
        return Err("prob must not be greater than one".into());
    }
    let s = sample_uniform_ubig_below_traced(denom.clone(), trace)?;
    let res = numer > s;
    if trace.bernoulli_rat.is_none() {
        trace.bernoulli_rat = Some(res);
    }
    Ok(res)
}

// Original CKS20 / OpenDP bernoulli_exp1. Hits the dashu division bug
fn sample_bernoulli_exp1_original_traced(x: RBig, trace: &mut SampleRow) -> Fallible<bool> {
    let mut k = UBig::ONE;
    let res = loop {
        if sample_bernoulli_rational_traced(x.clone() / &k, trace)? {
            k += UBig::ONE;
        } else {
            break k.clone() % 2u8 == 1;
        }
    };
    if trace.bernoulli_exp1.is_none() {
        trace.bernoulli_exp1 = Some(res);
    }
    Ok(res)
}

// Alternating Taylor series sampler
fn sample_bernoulli_exp1_taylor_traced(x: RBig, trace: &mut SampleRow) -> Fallible<bool> {
    if x.is_zero() {
        let res = true;
        if trace.bernoulli_exp1.is_none() {
            trace.bernoulli_exp1 = Some(res);
        }
        return Ok(res);
    }

    let mut buf = [0u8; 16];
    fill_bytes(&mut buf)?;
    let u_big = UBig::from_be_bytes(&buf);
    let grid = UBig::ONE << 128;
    let u_rb = RBig::from(u_big) / RBig::from(grid);

    let mut s_prev = RBig::ONE;
    let mut s_curr = RBig::ONE - x.clone();
    let mut k: u32 = 1;
    let mut term = x.clone();

    let res = loop {
        let (low, high) = if s_prev <= s_curr {
            (s_prev.clone(), s_curr.clone())
        } else {
            (s_curr.clone(), s_prev.clone())
        };
        if u_rb < low {
            break true;
        }
        if u_rb >= high {
            break false;
        }
        k += 1;
        term = (term * x.clone()) / RBig::from(UBig::from(k));
        let s_next = if k % 2 == 0 {
            s_curr.clone() + term.clone()
        } else {
            s_curr.clone() - term.clone()
        };
        s_prev = s_curr;
        s_curr = s_next;
    };

    if trace.bernoulli_exp1.is_none() {
        trace.bernoulli_exp1 = Some(res);
    }
    Ok(res)
}

// Same as the original CKS20 version, but we do the correct x/k division ourselves
fn gcd_ubig(mut a: UBig, mut b: UBig) -> UBig {
    while !b.is_zero() {
        let r = &a % &b;
        a = b;
        b = r;
    }
    a
}

fn div_rbig_by_ubig_exact(numer: &UBig, denom: &UBig, k: &UBig) -> RBig {
    assert!(!k.is_zero(), "division by zero");
    if numer.is_zero() {
        return RBig::ZERO;
    }
    let g = gcd_ubig(numer.clone(), k.clone());
    let n_red: UBig = numer / &g;
    let k_red: UBig = k / &g;
    RBig::from_parts(n_red.into(), denom * k_red)
}

fn sample_bernoulli_exp1_fixed_traced(x: RBig, trace: &mut SampleRow) -> Fallible<bool> {
    let (numer_signed, denom) = x.into_parts();
    let (sign, numer) = numer_signed.into_parts();
    if sign == Sign::Negative && !numer.is_zero() {
        return Err("x must be in [0, 1]".into());
    }
    if numer > denom {
        return Err("x must be in [0, 1]".into());
    }

    let mut k = UBig::ONE;
    let res = loop {
        let x_div_k = div_rbig_by_ubig_exact(&numer, &denom, &k);
        if sample_bernoulli_rational_traced(x_div_k, trace)? {
            k += UBig::ONE;
        } else {
            break k.clone() % 2u8 == 1;
        }
    };
    if trace.bernoulli_exp1.is_none() {
        trace.bernoulli_exp1 = Some(res);
    }
    Ok(res)
}

// The rest of the pipeline is parameterised on which bernoulli_exp1 we use
type B1Traced = fn(RBig, &mut SampleRow) -> Fallible<bool>;

fn sample_bernoulli_exp_traced(
    mut x: RBig,
    trace: &mut SampleRow,
    b1: B1Traced,
) -> Fallible<bool> {
    while x > RBig::ONE {
        if b1(RBig::ONE, trace)? {
            x -= RBig::ONE;
        } else {
            if trace.bernoulli_exp.is_none() {
                trace.bernoulli_exp = Some(false);
            }
            return Ok(false);
        }
    }
    let res = b1(x, trace)?;
    if trace.bernoulli_exp.is_none() {
        trace.bernoulli_exp = Some(res);
    }
    Ok(res)
}

fn sample_geometric_exp_slow_traced(
    x: RBig,
    trace: &mut SampleRow,
    b1: B1Traced,
) -> Fallible<UBig> {
    let mut k = UBig::ZERO;
    loop {
        if sample_bernoulli_exp_traced(x.clone(), trace, b1)? {
            k += UBig::ONE;
        } else {
            if trace.geom_slow.is_none() {
                trace.geom_slow = Some(k.clone());
            }
            return Ok(k);
        }
    }
}

fn sample_geometric_exp_fast_traced(
    x: RBig,
    trace: &mut SampleRow,
    b1: B1Traced,
) -> Fallible<UBig> {
    if x.is_zero() {
        if trace.geom_fast.is_none() {
            trace.geom_fast = Some(UBig::ZERO);
        }
        return Ok(UBig::ZERO);
    }
    let (numer_ib, denom) = x.into_parts();
    let (_sign, numer_mag) = numer_ib.into_parts();

    let mut u = sample_uniform_ubig_below_traced(denom.clone(), trace)?;
    while !sample_bernoulli_exp_traced(
        RBig::from_parts(u.as_ibig().clone(), denom.clone()),
        trace,
        b1,
    )? {
        u = sample_uniform_ubig_below_traced(denom.clone(), trace)?;
    }
    let v2 = sample_geometric_exp_slow_traced(RBig::ONE, trace, b1)?;
    let k = (v2 * denom + u) / numer_mag;
    if trace.geom_fast.is_none() {
        trace.geom_fast = Some(k.clone());
    }
    Ok(k)
}

fn sample_discrete_laplace_traced(
    scale: RBig,
    trace: &mut SampleRow,
    b1: B1Traced,
) -> Fallible<IBig> {
    if scale.is_zero() {
        return Ok(IBig::from(0));
    }

    let _ = sample_from_uniform_bytes_u128_traced(trace)?;

    let (numer_ib, denom) = scale.into_parts();
    let (_sign, numer_mag) = numer_ib.into_parts();
    let inv_scale = RBig::from_parts(denom.as_ibig().clone(), numer_mag);

    let _ = sample_bernoulli_rational_traced(inv_scale.clone(), trace)?;

    loop {
        let positive = sample_standard_bernoulli_traced(trace)?;
        let magnitude_ubig = sample_geometric_exp_fast_traced(inv_scale.clone(), trace, b1)?;
        let magnitude = IBig::from(magnitude_ubig);
        if positive || !magnitude.is_zero() {
            return Ok(if positive { magnitude } else { -magnitude });
        }
    }
}

fn samples_tracer(scale: RBig, b1: B1Traced) -> Fallible<SampleRow> {
    let mut trace = SampleRow::new();
    let dlap = sample_discrete_laplace_traced(scale, &mut trace, b1)?;
    trace.dlap = Some(dlap);
    Ok(trace)
}

pub fn run(out_dir: &Path) -> Fallible<()> {
    std::fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;

    let scale = RBig::from(10u8) / RBig::from(1u8);

    let methods: [(&str, B1Traced); 3] = [
        ("trace_original.csv", sample_bernoulli_exp1_original_traced),
        ("trace_taylor.csv",   sample_bernoulli_exp1_taylor_traced),
        ("trace_fixed.csv",    sample_bernoulli_exp1_fixed_traced),
    ];

    let header = "unif,unif_below,bernoulli,bernoulli_rat,bernoulli_exp1,\
                  bernoulli_exp,geom_slow,geom_fast,dlap";

    let mut writers: Vec<BufWriter<File>> = methods
        .iter()
        .map(|(name, _)| {
            let path = out_dir.join(name);
            BufWriter::new(File::create(&path).expect("create trace csv"))
        })
        .collect();
    for w in &mut writers {
        writeln!(w, "{}", header).map_err(|e| e.to_string())?;
    }

    for row_idx in 0..N_ROWS {
        for i in 0..methods.len() {
            let f = methods[i].1;
            let row = samples_tracer(scale.clone(), f)?;
            writeln!(writers[i], "{}", row.to_csv_line()).map_err(|e| e.to_string())?;
        }
        if (row_idx + 1) % 100_000 == 0 {
            eprintln!("    [traces] {} / {}", row_idx + 1, N_ROWS);
        }
    }

    for w in &mut writers {
        w.flush().map_err(|e| e.to_string())?;
    }

    Ok(())
}