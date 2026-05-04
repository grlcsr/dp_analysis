/*
    Performance and sampling-quality benchmark for the three bernoulli_exp1 variants
    Each iteration draws two discrete laplace samples per variant (centered in 10000 and 10001) and records the timing
    This is mimicking a counting query with sensitivity 1, where the two samples would be the noise added to the true 
    counts of two neighboring datasets 
    We want to see if there is any evident privacy leakage on the faulty original variant vs the fixed and the Taylor one
    We also want to see how the timings of the three variants compare
    The Taylor variant is expected to be slower than the fixed one, but it is interesting to see by how much
    Output is saved in dlap_{original,taylor,fixed}.csv and performance.csv
*/

use dashu::{
    integer::{IBig, Sign, UBig},
    rational::RBig,
};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::rng::{fill_bytes, Fallible};

const N_ROWS: usize = 1_000_000;

fn sample_uniform_ubig_below(upper: UBig) -> Fallible<UBig> {
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
            break &sample % &upper;
        }
    })
}

fn sample_standard_bernoulli() -> Fallible<bool> {
    let mut buffer = [0u8; 1];
    fill_bytes(&mut buffer)?;
    Ok(buffer[0] & 1 == 1)
}

fn sample_bernoulli_rational(prob: RBig) -> Fallible<bool> {
    let (numer_ib, denom) = prob.into_parts();
    let (sign, numer) = numer_ib.into_parts();
    if sign == Sign::Negative && !numer.is_zero() {
        return Err("numerator must not be negative".into());
    }
    if numer > denom {
        return Err("prob must not be greater than one".into());
    }
    let s = sample_uniform_ubig_below(denom.clone())?;
    Ok(numer > s)
}

// Original CKS20 / OpenDP bernoulli_exp1. Hits the dashu division bug.
fn sample_bernoulli_exp1_original(x: RBig) -> Fallible<bool> {
    let mut k = UBig::ONE;
    loop {
        if sample_bernoulli_rational(x.clone() / &k)? {
            k += UBig::ONE;
        } else {
            return Ok(k.clone() % 2u8 == 1);
        }
    }
}

// Alternating Taylor series sampler proposed in the paper.
fn sample_bernoulli_exp1_taylor(x: RBig) -> Fallible<bool> {
    if x.is_zero() {
        return Ok(true);
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

    loop {
        let (low, high) = if s_prev <= s_curr {
            (s_prev.clone(), s_curr.clone())
        } else {
            (s_curr.clone(), s_prev.clone())
        };
        if u_rb < low {
            return Ok(true);
        }
        if u_rb >= high {
            return Ok(false);
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
    }
}

// Same as the original CKS20 version, but we do the x/k division ourselves
// so we never hand dashu the input that triggers the bug.
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

fn sample_bernoulli_exp1_fixed(x: RBig) -> Fallible<bool> {
    let (numer_signed, denom) = x.into_parts();
    let (sign, numer) = numer_signed.into_parts();
    if sign == Sign::Negative && !numer.is_zero() {
        return Err("x must be in [0, 1]".into());
    }
    if numer > denom {
        return Err("x must be in [0, 1]".into());
    }

    let mut k = UBig::ONE;
    loop {
        let x_div_k = div_rbig_by_ubig_exact(&numer, &denom, &k);
        if sample_bernoulli_rational(x_div_k)? {
            k += UBig::ONE;
        } else {
            return Ok(k % 2u8 == 1);
        }
    }
}

// The rest of the pipeline is parameterised on which bernoulli_exp1 we use.
type BExp1 = fn(RBig) -> Fallible<bool>;

fn sample_bernoulli_exp(mut x: RBig, b1: BExp1) -> Fallible<bool> {
    while x > RBig::ONE {
        if b1(RBig::ONE)? {
            x -= RBig::ONE;
        } else {
            return Ok(false);
        }
    }
    b1(x)
}

fn sample_geometric_exp_slow(x: RBig, b1: BExp1) -> Fallible<UBig> {
    let mut k = UBig::ZERO;
    loop {
        if sample_bernoulli_exp(x.clone(), b1)? {
            k += UBig::ONE;
        } else {
            return Ok(k);
        }
    }
}

fn sample_geometric_exp_fast(x: RBig, b1: BExp1) -> Fallible<UBig> {
    if x.is_zero() {
        return Ok(UBig::ZERO);
    }
    let (numer_ib, denom) = x.into_parts();
    let (_sign, numer_mag) = numer_ib.into_parts();

    let mut u = sample_uniform_ubig_below(denom.clone())?;
    while !sample_bernoulli_exp(
        RBig::from_parts(u.as_ibig().clone(), denom.clone()),
        b1,
    )? {
        u = sample_uniform_ubig_below(denom.clone())?;
    }

    let v2 = sample_geometric_exp_slow(RBig::ONE, b1)?;
    Ok((v2 * denom + u) / numer_mag)
}

fn sample_discrete_laplace(scale: RBig, b1: BExp1) -> Fallible<IBig> {
    if scale.is_zero() {
        return Ok(IBig::from(0));
    }
    let (numer_ib, denom) = scale.into_parts();
    let (_sign, numer_mag) = numer_ib.into_parts();
    let inv_scale = RBig::from_parts(denom.as_ibig().clone(), numer_mag);

    loop {
        let positive = sample_standard_bernoulli()?;
        let magnitude_ubig = sample_geometric_exp_fast(inv_scale.clone(), b1)?;
        let magnitude = IBig::from(magnitude_ubig);
        if positive || !magnitude.is_zero() {
            return Ok(if positive { magnitude } else { -magnitude });
        }
    }
}

pub fn run(out_dir: &Path) -> Fallible<()> {
    std::fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;

    let scale = RBig::from(10u8) / RBig::from(1u8);

    let methods: [(&str, BExp1); 3] = [
        ("dlap_original.csv", sample_bernoulli_exp1_original),
        ("dlap_taylor.csv",   sample_bernoulli_exp1_taylor),
        ("dlap_fixed.csv",    sample_bernoulli_exp1_fixed),
    ];

    let mut writers: Vec<BufWriter<File>> = methods
        .iter()
        .map(|(name, _)| {
            let path = out_dir.join(name);
            BufWriter::new(File::create(&path).expect("create dlap csv"))
        })
        .collect();
    for w in &mut writers {
        writeln!(w, "d1,d2").map_err(|e| e.to_string())?;
    }

    let perf_path = out_dir.join("performance.csv");
    let mut perf =
        BufWriter::new(File::create(&perf_path).expect("create performance csv"));
    writeln!(perf, "original_ns,taylor_ns,fixed_ns").map_err(|e| e.to_string())?;

    let center_d1 = IBig::from(10_000);
    let center_d2 = IBig::from(10_001);

    for row_idx in 0..N_ROWS {
        let mut times = [0u128; 3];
        let mut samples: [(IBig, IBig); 3] = [
            (IBig::from(0), IBig::from(0)),
            (IBig::from(0), IBig::from(0)),
            (IBig::from(0), IBig::from(0)),
        ];

        for i in 0..3 {
            let f = methods[i].1;
            let t0 = Instant::now();
            let s1 = sample_discrete_laplace(scale.clone(), f)?;
            let s2 = sample_discrete_laplace(scale.clone(), f)?;
            times[i] = t0.elapsed().as_nanos();
            samples[i] = (&center_d1 + s1, &center_d2 + s2);
        }

        for i in 0..3 {
            writeln!(writers[i], "{},{}", samples[i].0, samples[i].1)
                .map_err(|e| e.to_string())?;
        }
        writeln!(perf, "{},{},{}", times[0], times[1], times[2])
            .map_err(|e| e.to_string())?;

        if (row_idx + 1) % 100_000 == 0 {
            eprintln!("    [perf]   {} / {}", row_idx + 1, N_ROWS);
        }
    }

    for w in &mut writers {
        w.flush().map_err(|e| e.to_string())?;
    }
    perf.flush().map_err(|e| e.to_string())?;

    Ok(())
}