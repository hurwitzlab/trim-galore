extern crate run_trim_galore;
use std::process;

fn main() {
    let config = match run_trim_galore::get_args() {
        Ok(c) => c,
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
    };

    if let Err(e) = run_trim_galore::run(config) {
        println!("Error: {}", e);
        process::exit(1);
    }
}
