/*
    Unified module that runs both the "traces" and "performance" modules, and loops over multiple runs
    The "traces" module generates the "trace_*.csv" files, while the "performance" module generates the "dlap_*.csv" and "performance.csv" 
    files
    The output is organized into separate directories for each run, under "./csvs/traces/run_X/" and "./csvs/perf/run_X/"
    The number of runs is determined by the constant "N_RUNS". Each run's duration is printed to the console, along with the 
    total duration for all runs at the end

    Compile in release mode for accurate performance measurements
*/

mod rng;
mod traces;
mod performance;

use std::path::PathBuf;
use std::time::Instant;

use rng::Fallible;

const N_RUNS: usize = 1;

fn main() -> Fallible<()> {
    let global_start = Instant::now();

    for run in 1..=N_RUNS {
        eprintln!("Run {} / {}", run, N_RUNS);

        let traces_dir = PathBuf::from(format!("./csvs/traces/run_{}", run));
        eprintln!("[traces] writing into {}", traces_dir.display());
        let t0 = Instant::now();
        traces::run(&traces_dir)?;
        eprintln!(
            "[traces] done in {:.2?}",
            t0.elapsed()
        );

        let perf_dir = PathBuf::from(format!("./csvs/perf/run_{}", run));
        eprintln!("[perf] writing into {}", perf_dir.display());
        let t0 = Instant::now();
        performance::run(&perf_dir)?;
        eprintln!(
            "[perf] done in {:.2?}",
            t0.elapsed()
        );
    }

    eprintln!(
        "\nAll {} runs finished in {:.2?}",
        N_RUNS,
        global_start.elapsed()
    );
    Ok(())
}