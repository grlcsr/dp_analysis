// This code is based on the OpenDP implementation of the discrete Laplace
// mechanism from the CKS20 paper (Canonne, Kamath, Steinke 2020).
// I've added tracing functionality to debug and understand the sampling
// behavior, specifically to track the first call to each primitive
// (uniform, bernoulli, geometric, etc.) and see where issues might arise.
//
// The core sampling algorithms are from OpenDP, the SampleRow struct and
// tracing logic are my additions for diagnostic purposes.


use dashu::{
    integer::{IBig, UBig, Sign},
    rational::RBig,
};
use std::fmt::Display;

type Fallible<T> = Result<T, String>;

// RNG wrapper around OpenSSL's random bytes or rand crate
pub fn fill_bytes(buffer: &mut [u8]) -> Fallible<()> {
    /*use openssl::rand::rand_bytes;

    rand_bytes(buffer).map_err(|e| format!("OpenSSL error: {:?}", e))?;*/
    use rand::Rng;
    rand::thread_rng().try_fill(buffer).map_err(|e| format!("OpenSSL error: {:?}", e))?;
    Ok(())
}

// Holds the first sample from each primitive we call,
// plus the final discrete Laplace result to trace the full sampling process
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

    // Convert to CSV line 
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

/*
 Traced sampling primitives
 Each one fills the trace row the first time it's called
*/

fn sample_from_uniform_bytes_u128_traced(trace: &mut SampleRow) -> Fallible<u128> {
    let mut buf = [0u8; 16];
    fill_bytes(&mut buf)?;
    let v = u128::from_be_bytes(buf);
    if trace.unif.is_none() {
        trace.unif = Some(v);
    }
    Ok(v)
}

// Sample uniform integer 
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
            // Record first sample we see
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

/// Original CKS20 Bernoulli(exp(-x)) for x in [0,1] traced
/*fn sample_bernoulli_exp1_traced(x: RBig, trace: &mut SampleRow) -> Fallible<bool> {
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
}*/

fn sample_bernoulli_exp1_traced(x: RBig, trace: &mut SampleRow) -> Fallible<bool> {
    use dashu::integer::UBig;
    println!("sample_bernoulli_exp1_traced called with x = {}", x);
    
    if x.is_zero() {
        let res = true;
        if trace.bernoulli_exp1.is_none() {
            trace.bernoulli_exp1 = Some(res);
        }
        return Ok(res);
    }

    // Draw uniform
    let mut buf = [0u8; 16];
    fill_bytes(&mut buf)?;
    let u_big = UBig::from_be_bytes(&buf);
    let grid = UBig::ONE << 128;
    let u_rb = RBig::from(u_big) / RBig::from(grid);
    println!("u_rb = {}", u_rb);

    // Compute alternating series for exp(-x):
    // S0 = 1
    // S1 = 1 - x
    // S2 = 1 - x + x^2/2
    // etc.
    // At each step k, we know exp(-x) is between S_(k-1) and S_k

    let mut s_prev = RBig::ONE;
    let mut s_curr = RBig::ONE - x.clone();
    let mut k: u32 = 1;
    let mut term = x.clone();

    let res = loop {
        // Current values
        let (low, high) = if s_prev <= s_curr {
            (s_prev.clone(), s_curr.clone())
        } else {
            (s_curr.clone(), s_prev.clone())
        };

        println!("k = {}, low = {}, high = {}", k, low, high);

        // Compare
        if u_rb < low {
            break true;
        }
        if u_rb >= high {
            break false;
        }

        // Compute next term
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

// Sample from Bernoulli(exp(-x)) for any x >= 0
fn sample_bernoulli_exp_traced(mut x: RBig, trace: &mut SampleRow) -> Fallible<bool> {
    println!("sample_bernoulli_exp_traced called with x = {}", x);
    
    while x > RBig::ONE {
        if sample_bernoulli_exp1_traced(RBig::ONE, trace)? {
            x -= RBig::ONE;
        } else {
            if trace.bernoulli_exp.is_none() {
                trace.bernoulli_exp = Some(false);
            }
            return Ok(false);
        }
    }
    
    let res = sample_bernoulli_exp1_traced(x, trace)?;
    if trace.bernoulli_exp.is_none() {
        trace.bernoulli_exp = Some(res);
    }
    Ok(res)
}

// Geometric distribution with parameter exp(-x) 
fn sample_geometric_exp_slow_traced(x: RBig, trace: &mut SampleRow) -> Fallible<UBig> {
    let mut k = UBig::ZERO;
    loop {
        if sample_bernoulli_exp_traced(x.clone(), trace)? {
            k += UBig::ONE;
        } else {
            if trace.geom_slow.is_none() {
                trace.geom_slow = Some(k.clone());
            }
            return Ok(k);
        }
    }
}

// Fast geometric sampler
fn sample_geometric_exp_fast_traced(x: RBig, trace: &mut SampleRow) -> Fallible<UBig> {
    if x.is_zero() {
        if trace.geom_fast.is_none() {
            trace.geom_fast = Some(UBig::ZERO.clone());
        }
        return Ok(UBig::ZERO);
    }

    let (numer_ib, denom) = x.into_parts();
    let (_sign, numer_mag) = numer_ib.into_parts();

    let mut u = sample_uniform_ubig_below_traced(denom.clone(), trace)?;
    while !sample_bernoulli_exp_traced(
        RBig::from_parts(u.as_ibig().clone(), denom.clone()),
        trace,
    )? {
        u = sample_uniform_ubig_below_traced(denom.clone(), trace)?;
    }

    let v2 = sample_geometric_exp_slow_traced(RBig::ONE, trace)?;

    let k = (v2 * denom + u) / numer_mag;

    if trace.geom_fast.is_none() {
        trace.geom_fast = Some(k.clone());
    }
    Ok(k)
}

/*
  Discrete Laplace sampler (main algorithm from CKS20)
*/

fn sample_discrete_laplace_traced(scale: RBig, trace: &mut SampleRow) -> Fallible<IBig> {
    if scale.is_zero() {
        return Ok(IBig::from(0));
    }

    // Log one uniform u128 at the start
    let _ = sample_from_uniform_bytes_u128_traced(trace)?;

    let (numer_ib, denom) = scale.into_parts();
    let (_sign, numer_mag) = numer_ib.into_parts();
    let inv_scale = RBig::from_parts(denom.as_ibig().clone(), numer_mag);

    let _ = sample_bernoulli_rational_traced(inv_scale.clone(), trace)?;

    loop {
        let positive = sample_standard_bernoulli_traced(trace)?;
        let magnitude_ubig = sample_geometric_exp_fast_traced(inv_scale.clone(), trace)?;
        let magnitude = IBig::from(magnitude_ubig);

        if positive || !magnitude.is_zero() {
            return Ok(if positive { magnitude } else { -magnitude });
        }
    }
}


pub fn samples_tracer(scale: RBig) -> Fallible<SampleRow> {
    let mut trace = SampleRow::new();
    let dlap = sample_discrete_laplace_traced(scale, &mut trace)?;
    trace.dlap = Some(dlap);
    Ok(trace)
}

fn main() {
    // Print CSV header
    println!(
        "unif,unif_below,bernoulli,bernoulli_rat,bernoulli_exp1,bernoulli_exp,geom_slow,geom_fast,dlap"
    );

    // Set up scale parameter
    // For epsilon = 0.1, we use scale = 10
    let scale = RBig::from(10u8) / RBig::from(1u8);
    let n_rows = 1000000;

    for _ in 0..n_rows {
        let row = samples_tracer(scale.clone()).expect("sampling failed");
        println!("{}", row.to_csv_line());
    }
}
