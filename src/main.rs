mod coset_table;
mod parser;

use parser::parse_ascii_coxeter_diagram;

fn main() {
    println!("Hello, world!");

    let result = parse_ascii_coxeter_diagram("x4o3o");
    match result {
        Err(e) => println!("{}", e),
        Ok(d) => println!("{:?}", d),
    }
}
