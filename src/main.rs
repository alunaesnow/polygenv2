extern crate blas;
extern crate openblas_src;

mod coset_table;
mod parser;

use coset_table::CosetTable;
use itertools::Itertools;
use ndarray::Array2;
use ndarray_linalg::{Norm, Solve};
use parser::{
    parse_ascii_coxeter_diagram, CombineMethod, NodeType, ParserError, PolytopeDescription,
};

pub struct Polytope {
    pub vertices: Array2<f64>,
    pub description: PolytopeDescription,
    pub dimension: usize,
    cts: Vec<Vec<CosetTable>>,
}

impl Polytope {
    pub fn new(diagram: &str, normalize: bool) -> Result<Self, ParserError> {
        let desc = parse_ascii_coxeter_diagram(diagram)?;

        if desc.combine_method != CombineMethod::None || desc.subpolytopes.len() != 1 {
            println!("Laced polytopes are currently not supported!");
        }
        let subdesc = &desc.subpolytopes[0];

        if subdesc
            .nodes
            .iter()
            .any(|&n| n != NodeType::Inactive && n != NodeType::Standard)
        {
            println!("snum/dual/holosnub nodes are currently not supported!");
        }

        let d = desc.symmetry_matrix.shape()[0];

        let mut relations = Vec::new();
        for i in 1..d {
            for j in 0..i {
                let mut relation = Vec::new();
                for _ in 0..desc.coxeter_matrix[[i, j]] {
                    relation.push(i);
                    relation.push(j);
                }
                relations.push(relation);
            }
        }

        let mut cts: Vec<Vec<CosetTable>> = Vec::new();
        for subd in 0..d {
            // generate `subd`-dimensional elements
            let mut subcts = Vec::new();
            for combo in (0..d).combinations(subd) {
                if subd > 0 {
                    // In general, the number of active mirrors should equal the number of orthogonal paris
                    // - the (sub_d - 1)th triangle number, for the mirrors to generate a face/cell/etc.
                    // wierd eh? i'll be honest i dont know how this works, i just discovered it by looking
                    // at the patterns in those two numbers hahaha
                    let subd_m1 = subd as i32 - 1;
                    let triangle_number = (subd_m1 * (subd_m1 + 1)) / 2;
                    let delta: i32 = triangle_number - subd as i32;

                    let active_count = combo
                        .iter()
                        .filter(|&&i| subdesc.nodes[i] != NodeType::Inactive)
                        .count() as i32;
                    let orthogonal_pair_count = combo
                        .iter()
                        .combinations(2)
                        .filter(|pair| {
                            (2.0 - desc.symmetry_matrix[[*pair[0], *pair[1]]]).abs() < 1e-15
                        })
                        .count() as i32;

                    if !(orthogonal_pair_count == active_count + delta) {
                        continue;
                    }
                }

                let mut subgens = Vec::new();
                for s in 0..d {
                    if combo
                        .iter()
                        .all(|&v| (2.0 - desc.symmetry_matrix[[v, s]]).abs() < 1e-15)
                        && subdesc.nodes[s] == NodeType::Inactive
                    {
                        subgens.push(vec![s]);
                    }
                }
                for i in combo {
                    subgens.push(vec![i]);
                }

                let mut ct = CosetTable::new(d, true, relations.clone(), subgens);
                ct.solve(10000);
                subcts.push(ct);
            }

            cts.push(subcts);
        }

        let normals = desc.place_mirrors();
        let mut inital_vertex = normals.solve_into(subdesc.offsets.clone()).unwrap();
        if normalize {
            inital_vertex /= inital_vertex.norm_l2();
        }

        let vct = &cts[0][0];
        let mut vertices = Array2::<f64>::zeros((vct.len(), d));
        let reprs = vct.get_representatives();

        for (i, repr) in reprs.iter().enumerate() {
            let mut v = inital_vertex.clone();
            for &g in repr {
                let n = normals.row(g);
                let dot = v.dot(&n);
                v = v - n.to_owned() * (2.0 * dot);
            }
            vertices.row_mut(i).assign(&v);
        }

        Ok(Self {
            vertices,
            description: desc,
            dimension: d,
            cts,
        })
    }

    pub fn elem_counts(&self) -> Vec<usize> {
        self.cts
            .iter()
            .map(|cts| cts.iter().fold(0, |acc, ct| acc + ct.len()))
            .collect()
    }

    // pub fn elem_constituents(&self, parent: usize, child: usize) {

    // }
}

fn main() {
    // "x3o3o3o*b3o"
    let poly = Polytope::new("x4o3o", true).unwrap();

    println!("Element counts: {:?}", poly.elem_counts());

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
