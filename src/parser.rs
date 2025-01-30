use std::{collections::HashMap, f64::consts::PI, fmt::Display};

const VALID_NODES: [char; 21] = [
    'o', 'v', 'x', 's', 'm', 'p', 'ß', 'q', 'f', 'h', 'k', 'u', 'w', 'F', 'Q', 'd', 'V', 'U', 'A',
    'X', 'B',
];
const INT_CHARS: [char; 10] = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'];
const OTHER_CHARS: [char; 3] = [' ', '/', '*'];

lazy_static::lazy_static! {
    static ref EDGE_LENGTHS: HashMap<char, f64> = {
        let mut m = HashMap::new();
        m.insert('o', 0.0);
        m.insert('v', (5.0_f64.sqrt() - 1.0) / 2.0);
        m.insert('x', 1.0);
        m.insert('s', 1.0);
        m.insert('m', 1.0);
        m.insert('p', 1.0);
        m.insert('ß', 1.0);
        m.insert('q', 2.0_f64.sqrt());
        m.insert('f', (5.0_f64.sqrt() + 1.0) / 2.0);
        m.insert('h', 3.0_f64.sqrt());
        m.insert('k', (2.0_f64.sqrt() + 2.0).sqrt());
        m.insert('u', 2.0);
        m.insert('w', 2.0_f64.sqrt() + 1.0);
        m.insert('F', (5.0_f64.sqrt() + 3.0) / 2.0);
        m.insert('Q', 2.0 * 2.0_f64.sqrt());
        m.insert('d', 3.0);
        m.insert('V', 1.0 + 5.0_f64.sqrt());
        m.insert('U', 2.0 + 2.0_f64.sqrt());
        m.insert('A', (5.0 + 5.0_f64.sqrt()) / 4.0);
        m.insert('X', 1.0 + 2.0 * 2.0_f64.sqrt());
        m.insert('B', 2.0 + 5.0_f64.sqrt());

        m
    };
}

#[derive(PartialEq, Debug)]
pub enum CombineMethod {
    None,
    Compound,
    LacePrism,
    LaceTower,
    LaceTegum,
    LaceRing,
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum NodeType {
    Inactive,
    Standard,
    Snub,
    Dual,
    DualSnub,
    Holosnub,
}

#[derive(Debug)]
/// Returned by `parse_ascii_coxeter_diagram`, this gives imformation needed
/// to construct the polytope described by the diagram
pub struct PolytopeDescription {
    /// The diagram this information comes from
    pub diagram: String,
    /// Specifies how to combine all the polytopes found in subpolytopes
    pub combine_method: CombineMethod,
    /// Edge length matrix
    /// M[i][j] = (n, d) means the i-th and j-th nodes are connected by an edge
    /// with nominator `n` and denumerator `d`
    pub edge_matrix: Vec<Vec<(u32, u32)>>,
    /// Subpolytopes to combine, if `combine_method` is None this will have size 1
    pub subpolytopes: Vec<SubpolytopeDescription>,
}

#[derive(Clone, Debug)]
pub struct SubpolytopeDescription {
    /// Offsets from the mirrors to place the initial vertex v0
    pub offsets: Vec<f64>,
    /// node types
    pub nodes: Vec<NodeType>,
}

#[derive(Debug)]
pub struct ParserError {
    /// diagram that was being parsed
    diagram: String,
    /// position in diagram where the error began
    err_start: usize,
    /// position in diagram where the error ends
    err_end: usize,
    /// error description
    description: String,
}

impl Display for ParserError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (l, tmp) = self.diagram.split_at(self.err_start);
        let (m, r) = tmp.split_at(self.err_end - l.len());

        write!(
            f,
            "Error parsing diagram at: {}[{}]{}\n- {}",
            l, m, r, self.description,
        )
    }
}

impl ParserError {
    fn new(s: usize, e: usize, diagram: &str, description: String) -> Self {
        Self {
            diagram: diagram.to_string(),
            err_start: s,
            err_end: e,
            description,
        }
    }

    fn new_whole(diagram: &str, description: String) -> Self {
        Self::new(0, diagram.len(), diagram, description)
    }
}

/// Parses a coxter diagram in plaintext format (e.g. x4o3o), and returns information
/// used to construct the polytope.
///
/// A description of how the plaintext format works can be found here:
/// https://polytope.miraheze.org/wiki/Coxeter_diagram
///
/// The supported features are:
/// - Active/Inactive nodes using 'x' and 'o' respectively
///   - e.g. x4o3o (cube)
/// - Spaces to denote 2-length edges
///   - e.g. x4o3o = x x x (cube)
/// - Fractional edges
///   - e.g. x5/2o (pentagram)
/// - Virtual nodes
///   - e.g. x3o3o3*a (triangular tiling)
///   - e.g. x3o3o*b3o (hexadecahedron)
/// - Different edge lengths
///   - use o,v,x,q,f,h,k,u,w,F,Q,d,V,U,A,X,B
///   - e.g. d x (3 by 1 rectangle)
/// - Snub nodes using 's'
///   - e.g. s4s3s (snub cube)
/// - Dual nodes using 'm' in place of 'x' and 'p' in place of 's'
///   - e.g. m4m3m (octahedron)
/// - Holosnub nodes using 'ß'
///   - e.g. ß5o (pentagram)
/// - Compounds with matching symmetry
///   - e.g. xo4oo3oq (compound of cube and octahedron)
/// - Lace prisms/simplices by adding '&#x' to the end of the diagram
///   - e.g. xx3ox&#x (triangular cupola)
/// - Lace towers by adding '&#xt' to the end of the diagram
///   - e.g. xxo3oxx&#xt (cubeoctahedron)
/// - Lace tegums by adding '&#m' to the end of the diagram
///   - e.g. TODO
/// - Lace rings by adding '&#xr' to the end of the diagram
///   - e.g. xxxx5oooo&#xr (square-pentagonal duoprism)
///
/// For a description of the returned data see the accompanying `PolytopeDescription` type.
///
/// @param diagram The coxeter diagram in plaintext format
/// @returns See the `PolytopeDescription` type exported along with this function
pub fn parse_ascii_coxeter_diagram(mut diagram: &str) -> Result<PolytopeDescription, ParserError> {
    let original_diagram = diagram.to_owned();

    // Determine if this is a lace diagram
    let mut combine_method: CombineMethod = CombineMethod::None;

    let idx = diagram.find("&#");
    if idx.is_some() {
        let split = diagram.split("&#").collect::<Vec<&str>>();
        if split.len() > 2 {
            return Err(ParserError::new(
                idx.unwrap(),
                diagram.len(),
                diagram,
                "Multiple lace indicators \"&#\" found, only one should be present".into(),
            ));
        }

        combine_method = match split[1] {
            "x" => CombineMethod::LacePrism,
            "xt" => CombineMethod::LaceTower,
            "m" => CombineMethod::LaceTegum,
            "xr" => CombineMethod::LaceRing,
            "" => {
                return Err(ParserError::new(
                    idx.unwrap(),
                    diagram.len(),
                    diagram,
                    "Lace indicator \"&#\" found, but no lace type specified after it".into(),
                ))
            }
            etc => {
                return Err(ParserError::new(
                    idx.unwrap(),
                    diagram.len(),
                    diagram,
                    format!("Unknown lace type \"{}\"", etc),
                ))
            }
        };

        diagram = split[0]; // Remove lace indicator from diagram
    }

    let chars = diagram.chars().collect::<Vec<char>>();

    // Determine how many nodes there are, check all characters are valid, and count
    // the size of each node group (e.g. xo3o has 2 groups xo of size 1 and o of size 3)
    let mut node_group_sizes: Vec<usize> = Vec::new();
    let mut new_node_group = true;
    let mut node_count = 0;
    for i in 0..chars.len() {
        let c = chars[i];
        let prev_c = if i > 0 { chars[i - 1] } else { '#' };
        if prev_c != '*' {
            // any character following '*' is allowed
            if VALID_NODES.contains(&c) {
                if new_node_group {
                    node_group_sizes.push(0);
                    new_node_group = false;
                }
                let n = node_group_sizes.len();
                node_group_sizes[n - 1] += 1;
                node_count += 1;
            } else if !OTHER_CHARS.contains(&c) && !INT_CHARS.contains(&c) {
                return Err(ParserError::new(
                    i,
                    i + 1,
                    diagram,
                    format!(
                        "Invalid character \"{}\". Valid node types are: {}. Valid edge characters are: {},{}, and any character straight after a '*'",
                        c,
                        VALID_NODES
                            .iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                        INT_CHARS
                            .iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                        OTHER_CHARS
                            .iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join(",") ),
                ));
            } else {
                new_node_group = true;
            }
        }
    }

    // Check that the node group sizes are valid
    let subpolytope_count = node_group_sizes[0];
    if !node_group_sizes
        .iter()
        .all(|size| *size == subpolytope_count)
    {
        return Err(ParserError::new_whole(
            diagram,
            format!(
                "All nodes groups must have the same size, sizes were: {}",
                node_group_sizes
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        ));
    }
    let is_laced = combine_method == CombineMethod::LacePrism
        || combine_method == CombineMethod::LaceRing
        || combine_method == CombineMethod::LaceTegum
        || combine_method == CombineMethod::LaceTower;
    if subpolytope_count == 1 && is_laced {
        return Err(ParserError::new_whole(
            diagram,
            "Laced diagrams should have node groups of minimum size 2".into(),
        ));
    } else if subpolytope_count > 1 && !is_laced {
        // As the groups have size > 1, but there was no lace indicator, this is a
        // compound polytope
        combine_method = CombineMethod::Compound;
    }

    node_count = node_count / subpolytope_count;

    // Prepare the data structures
    let mut M = vec![vec![(2u32, 1u32); node_count]; node_count];

    let mut subpolytopes = vec![
        SubpolytopeDescription {
            nodes: Vec::new(),
            offsets: vec![0.0; node_count]
        };
        subpolytope_count
    ];

    // Parse the diagram

    let vnode_target = |i: usize| -> Result<usize, ParserError> {
        let target = u32::from(chars[i + 1]) as usize - 97;
        // case < 0 is handled by an overflow
        if target >= node_count {
            return Err(ParserError::new(
                i,
                i + 2,
                diagram,
                "Virtual node must point to a valid node".into(),
            ));
        }
        Ok(target as usize)
    };

    let get_int = |i: usize| -> Result<(u32, usize), ParserError> {
        let mut i_to = i;
        while chars[i_to] >= '0' && chars[i_to] <= '9' {
            i_to += 1;
        }
        let v = diagram[i..i_to].parse::<u32>();
        Err(ParserError::new(
            i,
            i_to,
            diagram,
            match v {
                Err(_) => "Failed to parse integer",
                Ok(v) => {
                    if v >= 2 {
                        return Ok((v, i_to));
                    }
                    "Integer must be at least 2"
                }
            }
            .into(),
        ))
    };

    if !VALID_NODES.contains(&chars[0]) {
        return Err(ParserError::new(
            0,
            0,
            diagram,
            "First character must be a node".into(),
        ));
    }

    let mut prev_was_edge = true;
    let mut from: usize = 0;
    let mut to: usize = 0;
    let mut i: usize = 0;
    // println!("diagram: {}", diagram);
    while i < diagram.len() {
        // println!("{}", i);
        if !prev_was_edge {
            // println!("parsing edge");
            // If a vnode exists before the edge value, then we need to jump
            // to the vnode's target node
            if chars[i] == '*' {
                from = vnode_target(i)?;
                i += 2;
            }

            // Parse edge value
            let mut n = 2;
            let mut d = 1;
            if chars[i] != ' ' {
                (n, i) = get_int(i)?;

                if chars[i] == '/' {
                    i += 1;
                    (d, i) = get_int(i)?;
                }
            } else {
                i += 1;
            }

            // If a vnode exists after the edge value, then the edge
            // points to the vnode's target node
            if chars[i] == '*' {
                to -= 1;
                from = vnode_target(i)?;
                i += 2;
            }

            if to >= node_count {
                return Err(ParserError::new(
                    i,
                    i + 1,
                    diagram,
                    "Missing a node here".into(),
                ));
            }

            M[from][to] = (n, d);
            M[to][from] = (n, d);
            prev_was_edge = true;
        } else {
            // println!("parsing nodes");
            // We are processing nodes
            for s in 0..subpolytope_count {
                let subpoly = &mut subpolytopes[s];
                let c = chars[i + s];
                if VALID_NODES.contains(&c) {
                    subpoly.offsets[to] = EDGE_LENGTHS.get(&c).unwrap() / 2.0;
                    subpoly.nodes.push(match c {
                        'o' => NodeType::Inactive,
                        's' => NodeType::Snub,
                        'm' => NodeType::Dual,
                        'p' => NodeType::DualSnub,
                        'ß' => NodeType::Holosnub,
                        _ => NodeType::Standard,
                    });
                } else {
                    return Err(ParserError::new(
                        i + s,
                        i + s + 1,
                        diagram,
                        "Expected a node, is this a valid node type?".into(),
                    ));
                }
            }

            i += subpolytope_count;
            from = to;
            to += 1;
            prev_was_edge = false;
        }
    }

    Ok(PolytopeDescription {
        diagram: original_diagram,
        combine_method,
        edge_matrix: M,
        subpolytopes,
    })
}

// /**
//  * Places a set of mirrors given a symmetry matrix S, where S[i][j] = k means
//  * the i-th and j-th mirrors have an angle of PI/k between them. Returns the
//  * normals of the mirrors.
//  *
//  * Taken from: https://github.com/neozhaoliang/pywonderland/blob/master/src/polytopes/polytopes/helpers.py
//  */
// function placeMirrors(symmetryMatrix: number[][]) {
//     const C = symmetryMatrix.map((row) => row.map((x) => -Math.cos(Math.PI / x)));
//     const M = C.map((row) => row.map(() => 0));
//     const n = M.length;
//     // The first normal vector is simply (1, 0, ..., 0)
//     M[0][0] = 1;
//     // In the i-th row, the j-th entry can be computed via the (j, j) entry
//     for (let i = 1; i < n; i++) {
//       for (let j = 0; j < i; j++) {
//         const mj_colonj = M[j].slice(0, j);
//         const mi_colonj = M[i].slice(0, j);
//         M[i][j] = (C[i][j] - vm.dot(mj_colonj, mi_colonj)) / M[j][j];
//       }
//       const mi_coloni = M[i].slice(0, i);
//       M[i][i] = Math.sqrt(1 - vm.dot(mi_coloni, mi_coloni));
//     }
//     return M;
//   }

impl PolytopeDescription {
    /// Places a set of mirrors given a symmetry matrix S, where S[i][j] = k means
    /// the i-th and j-th mirrors have an angle of PI/k between them. Returns the
    /// normals of the mirrors. Note these are (n-1)-dimensional mirrors, so in 2d
    /// space they are lines, in 3d they are planes, in 4d volumes etc...
    ///
    /// Taken from: https://github.com/neozhaoliang/pywonderland/blob/master/src/polytopes/polytopes/helpers.py
    pub fn place_mirrors(&self) -> Vec<Vec<f64>> {
        let n = self.edge_matrix.len();
        let C: Vec<Vec<f64>> = self
            .edge_matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|v| -(PI / (f64::from(v.0) / f64::from(v.1))).cos())
                    .collect()
            })
            .collect();
        let mut M = vec![vec![0f64; n]; n];
        // The first normal vector is simply (1, 0, ..., 0)
        M[0][0] = 1.0;
        // In the i-th row, the j-th entry can be computed via the (j, j) entry
        for i in 0..n {
            for j in 0..i {
                let mut dot = 0.0;
                for k in 0..j {
                    dot += M[j][k] * M[i][k];
                }
                M[i][j] = (C[i][j] - dot) / M[j][j];
            }
            let mut dot = 0.0;
            for k in 0..i {
                dot += M[i][k] * M[i][k];
            }
            M[i][i] = (1.0 - dot).sqrt();
        }

        M
    }
}
