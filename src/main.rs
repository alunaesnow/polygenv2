mod coset_table;
mod parser;

use coset_table::CosetTable;
use parser::parse_ascii_coxeter_diagram;

fn main() {
    println!("Hello, world!");

    let result = parse_ascii_coxeter_diagram("x4o3o");
    match result {
        Err(e) => println!("{}", e),
        Ok(d) => println!("{:?}", d),
    }

    // let mut ct = CosetTable::new(
    //     4,
    //     false,
    //     vec![vec![0, 0], vec![2, 2, 2], vec![0, 2, 0, 2, 0, 2]],
    //     vec![vec![0, 2]],
    // );
    // ct.solve(1000);
    // ct.print_table();

    // let mut ct = CosetTable::new(
    //     4,
    //     false,
    //     vec![vec![0, 0, 2, 2], vec![0, 0, 0, 2, 2, 2, 2, 2]],
    //     vec![],
    // );
    // ct.solve(1000);
    // ct.print_table();

    // let mut ct = CosetTable::new(
    //     6,
    //     false,
    //     vec![
    //         vec![0; 11],
    //         vec![2, 2],
    //         vec![4, 4],
    //         vec![0, 2, 0, 2, 0, 2],
    //         vec![0, 4, 0, 4, 0, 4],
    //         vec![2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4],
    //         vec![0, 0, 2, 4, 2, 4, 0, 3, 5, 3, 5],
    //     ],
    //     vec![],
    // );
    // ct.solve(10000000);
    // println!("len is: {}", ct.len());
}
