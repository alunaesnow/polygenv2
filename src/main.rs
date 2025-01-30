mod coset_table;
mod parser;

use std::collections::HashSet;

use coset_table::CosetTable;
use itertools::Itertools;
use parser::{parse_ascii_coxeter_diagram, CombineMethod, NodeType};

// struct Polytope {

// }

fn main() {
    let diagram = "x3o3o3o*b3o";
    println!("Parsing coxeter diagram `{}`...", diagram);

    let desc = match parse_ascii_coxeter_diagram(&diagram) {
        Err(e) => {
            println!("Failed to parse coxeter diagram:\n{}", e);
            return;
        }
        Ok(d) => {
            println!("Parsed data: {:?}", d);
            d
        }
    };

    if desc.combine_method != CombineMethod::None || desc.subpolytopes.len() != 1 {
        println!("Laced polytopes are currently not supported!");
    }
    let polydesc = &desc.subpolytopes[0];

    if polydesc
        .nodes
        .iter()
        .any(|&n| n != NodeType::Inactive && n != NodeType::Standard)
    {
        println!("snum/dual/holosnub nodes are currently not supported!");
    }

    let d = desc.edge_matrix.len();
    println!("Dimension is: {}", d);

    let mut relations = Vec::new();
    for i in 1..d {
        for j in 0..i {
            let mut relation = Vec::new();
            for _ in 0..desc.edge_matrix[i][j].0 {
                relation.push(i);
                relation.push(j);
            }
            relations.push(relation);
        }
    }

    println!("relations: {:?}", relations);

    let mut counts = Vec::new();
    let mut cts: Vec<Vec<CosetTable>> = Vec::new();
    for sub_d in 0..d {
        println!("generating {}d elements...", sub_d);

        let mut subcts = Vec::new();
        let mut count = 0;
        for combo in (0..d).combinations(sub_d) {
            println!("  testing combination: {:?}", combo);
            if sub_d > 0 {
                let sub_d_m1 = sub_d - 1;
                let binomial_coefficient = (sub_d_m1 * (sub_d_m1 + 1)) / 2;
                let delta = binomial_coefficient as i32 - sub_d as i32;

                // In general, if active mirror count <= number of orthogonal pairs, then
                // the n mirrors generate a face/cell/etc
                let mut active_count = 0;
                let mut orthogonal_pair_count = 0;
                for &i in &combo {
                    active_count += (polydesc.nodes[i] != NodeType::Inactive) as i32;
                }
                for pair in combo.iter().combinations(2) {
                    let val = f64::from(desc.edge_matrix[*pair[0]][*pair[1]].0)
                        / f64::from(desc.edge_matrix[*pair[0]][*pair[1]].1);
                    orthogonal_pair_count += ((2.0 - val).abs() < 1e-15) as i32;
                }
                println!(
                    "    active: {}, orthogonal: {}",
                    active_count, orthogonal_pair_count
                );

                if !(orthogonal_pair_count == active_count + delta) {
                    continue;
                }
            }
            println!("    was valid");

            let mut subgens = Vec::new();
            for s in 0..d {
                if combo.iter().all(|&v| {
                    let val =
                        f64::from(desc.edge_matrix[v][s].0) / f64::from(desc.edge_matrix[v][s].1);
                    (2.0 - val).abs() < 1e-15
                }) {
                    if polydesc.nodes[s] == NodeType::Inactive {
                        subgens.push(vec![s]);
                    }
                }
            }
            for i in combo {
                subgens.push(vec![i]);
            }
            println!("    using subgens: {:?}", subgens);

            let mut ct = CosetTable::new(d, true, relations.clone(), subgens);
            ct.solve(10000);
            count += ct.len();
            subcts.push(ct);
        }

        println!("  generated {} elements", count);
        counts.push(count);
        cts.push(subcts);
    }

    println!("counts: {:?}", counts);

    for i in 0..10 {
        let i = i - 1;
        println!("{}: {}", i + 1, (i * (i + 1)) / 2);
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
