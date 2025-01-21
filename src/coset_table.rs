// C=(t,x,p,n,M)
// M = max number of cosets allowed
// 1 <= n <= M
//
// t: map from coset -> generator ? AHH THE MINIMUM REPRESENTATION
//
// p: map coset -> coset, p(coset) <= coset
// aka p is the equivalence
//
// x: the table?? no?
// map of a^x where a is coset, x is generator. so yes the table! i think?
//
// omega: live cosets set, where p(coset) = coset. always 1 in omega

use std::{collections::VecDeque, rc::Rc};

/// An implementation of the HLT (Hasselgrove, Leech and Trotter) strategy, for solving
/// the coset enumeration problem.
///
/// based on chapter 5 of:
/// Holt, D.F., Eick, B., & O'Brien, E.A. (2005). Handbook of Computational Group Theory (1st ed.).
/// Chapman and Hall/CRC. https://doi.org/10.1201/9781420035216
///
/// most of the symbols have been renamed to be more easily understandable
///
/// Also inspired by the following python implementation:
/// https://github.com/neozhaoliang/pywonderland/blob/master/src/polytopes/polytopes/todd_coxeter.py
pub struct CosetTable {
    /// number of generators (including inverses)
    n_gen: usize,
    /// true if the generators are self-inverse
    is_coxeter: bool,
    /// generator relations represented as integer arrays (R in book)
    /// odd integers are inverses
    /// the inverse of 0 is 1, inverse of 2 is 3, inverse of 4 is 5 etc...
    /// e.g. a^2 = b^3 = (ab)^3 = I
    /// is represented as [[0, 0], [2, 2, 2], [0, 2, 0, 2, 0, 2]]
    relations: Rc<Vec<Vec<usize>>>,
    /// generators of stabilising subgroup H represented as integer arrays (Y in book)
    /// (same format as relations)
    subgens: Rc<Vec<Vec<usize>>>,

    /// holds the equivalence classes of the cosets (p in book). equivalence[i] = j  
    /// means that coset i is equivalent to j. If equivalence[i] = i then coset i is
    /// "alive", otherwise it is considered "dead". note equivalence[i] <= i always holds.
    equivalence: Vec<usize>,
    /// queue holding "dead" cosets for processing. this is needed as when handling
    /// a coincidence, more coincidences may be discovered.
    dead_queue: VecDeque<usize>,
    /// the coset table itself (χ in book)
    table: Vec<Vec<Option<usize>>>,
    /// true if the coset table is solved
    is_solved: bool,
}

impl CosetTable {
    pub fn new(
        n_gen: usize,
        is_coxeter: bool,
        relations: Vec<Vec<usize>>,
        subgens: Vec<Vec<usize>>,
    ) -> Self {
        for relation in relations.iter() {
            for generator in relation {
                if generator > &n_gen {
                    todo!("error handling");
                }
            }
        }

        for subgen in subgens.iter() {
            for generator in subgen {
                if generator > &n_gen {
                    todo!("error handling");
                }
            }
        }

        // TODO: checks
        Self {
            n_gen,
            is_coxeter,
            relations: Rc::new(relations),
            subgens: Rc::new(subgens),
            equivalence: vec![0],
            dead_queue: VecDeque::new(),
            table: vec![vec![None; n_gen]],
            is_solved: false,
        }
    }

    /// returns the inverse of a generator
    fn inv(&self, gen: usize) -> usize {
        if self.is_coxeter {
            gen
        } else {
            gen ^ 1
        }
    }

    /// mark that coset_a · gen = coset_b
    fn deduction(&mut self, coset_a: usize, coset_b: usize, gen: usize) {
        let g_inv = self.inv(gen);
        self.table[coset_a][gen] = Some(coset_b);
        self.table[coset_b][g_inv] = Some(coset_a);
    }

    // create a new coset new_coset, deducing coset_a · gen (x) = new_coset
    fn define(&mut self, coset_a: usize, gen: usize) {
        let new_coset = self.table.len();
        self.table.push(vec![None; self.n_gen]);
        self.equivalence.push(new_coset);
        self.deduction(coset_a, new_coset, gen);
    }

    /// returns true if a coset is alive
    fn is_alive(&self, coset: usize) -> bool {
        self.equivalence[coset] == coset
    }

    /// returns the smallest (index wise) coset that is equivalent to `coset`.
    /// this method will also compress the equivalence chain as an optimisation
    fn rep(&mut self, coset: usize) -> usize {
        let mut equivalent = coset; // (λ = κ)
        while equivalent != self.equivalence[equivalent] {
            equivalent = self.equivalence[equivalent];
        }

        // path compression
        // e.g. if coset is 6, and equivalence is [0, 0, 2, 3, 1, 5, 4]
        // we follow the path 6 => 4 => 1 => 0
        // so we can set 6 and 4 to both point to 0 too
        // with the resulting equivalcence being compressed to [0, 0, 2, 3, 0, 5, 0]
        let mut c = coset;
        while c != equivalent {
            let tmp = self.equivalence[c];
            self.equivalence[c] = equivalent;
            c = tmp;
        }

        equivalent
    }

    /// given two equivalent cosets, declare the larger one to be dead and add
    /// it to the dead queue
    fn merge(&mut self, coset_a: usize, coset_b: usize) {
        let a = self.rep(coset_a);
        let b = self.rep(coset_b);
        if a != b {
            // Remember equivalence[i] <= i must hold
            let min = a.min(b);
            let max = a.max(b);
            self.equivalence[max] = min;
            self.dead_queue.push_back(max);
        }
    }

    /// a coincidence occurs when we discover that two cosets are in fact the same coset.
    fn coincidence(&mut self, coset_a: usize, coset_b: usize) {
        self.merge(coset_a, coset_b);

        // process all dead cosets
        // the queue holds only one coset at this point, but more may be added
        while let Some(dead) = self.dead_queue.pop_front() {
            // for each coset `next` this dead one connects to
            for g in 0..self.n_gen {
                if let Some(next) = self.table[dead][g] {
                    let g_inv = self.inv(g);
                    // delete `next`'s reference to the dead coset
                    self.table[next][g_inv] = None;

                    // copy the existing information to the representative of the dead coset
                    // but watch out for new coincidences. next may be dead already, so
                    // we need to watch out for that too
                    let dead_rep = self.rep(dead);
                    let next_rep = self.rep(next);

                    let dead_existing = self.table[next_rep][g_inv];
                    let next_existing = self.table[dead_rep][g];

                    if next_existing.is_some_and(|v| v != next_rep) {
                        // we have a coincidence
                        self.merge(next_rep, next_existing.unwrap());
                    } else if dead_existing.is_some_and(|v| v != dead_rep) {
                        // we have a coincidence
                        self.merge(dead_rep, dead_existing.unwrap());
                    } else {
                        // No coincidence, so we can just copy the information
                        self.deduction(dead_rep, next_rep, g)
                    }
                }
            }
        }
    }

    fn scan_and_fill(&mut self, coset_a: usize, word: &Vec<usize>) {
        let mut lv = coset_a; // left value (f)
        let mut rv = coset_a; // right value (b)
        let mut lp = 0; // left position (i)
        let mut rp = word.len() - 1; // right position (j)

        let mut iter = 0;
        while iter < word.len() {
            iter += 1;

            // scan forwards as far as possible
            while lp <= rp && self.table[lv][word[lp]].is_some() {
                lv = self.table[lv][word[lp]].unwrap();
                lp += 1;
            }

            // if complete check for coincidence, then stop
            if lp > rp {
                if lv != rv {
                    self.coincidence(lv, rv);
                }
                return;
            }

            // scan backwards as far as possible
            while rp >= lp && self.table[rv][self.inv(word[rp])].is_some() {
                rv = self.table[rv][self.inv(word[rp])].unwrap();
                rp -= 1;
            }

            // if rp and lp overlap, a coincidence has been found
            if rp < lp {
                self.coincidence(lv, rv);
                return;
            }

            // if they are about to meet, a deduction can be made
            if rp == lp {
                self.deduction(lv, rv, word[lp]);
                return;
            }

            // otherwise, create a new coset and continue scanning forwards
            self.define(lv, word[lp]);
        }

        panic!("over scanned, this should never happen");
    }

    // run the HLT strategy with the provided iteration limit
    fn coset_numerator_r(&mut self, max_iter: usize) {
        for word in self.subgens.clone().iter() {
            self.scan_and_fill(0, word);
        }

        let mut latest_coset = 0;
        while latest_coset < self.table.len() && latest_coset < max_iter {
            for relation in self.relations.clone().iter() {
                if !self.is_alive(latest_coset) {
                    break;
                }
                self.scan_and_fill(latest_coset, relation);
            }

            if self.is_alive(latest_coset) {
                for g in 0..self.n_gen {
                    if self.table[latest_coset][g].is_none() {
                        self.define(latest_coset, g);
                    }
                }
            }

            latest_coset += 1;
        }

        if latest_coset >= max_iter {
            todo!("Error handling: Iteration limit exceeded");
        }
    }

    /// delete all dead cosets in the table, renumbering the live ones
    fn compress(&mut self) {
        let mut n_live_cosets = 0;
        for coset in 0..self.table.len() {
            if self.is_alive(coset) {
                if n_live_cosets != coset {
                    for g in 0..self.n_gen {
                        if let Some(next) = self.table[coset][g] {
                            if next == coset {
                                self.table[n_live_cosets][g] = Some(n_live_cosets);
                            } else {
                                self.deduction(n_live_cosets, next, g);
                            }
                        } else {
                            todo!("Error handling: Compress should not be called before the coset table is solved");
                        }
                    }
                }
                n_live_cosets += 1; // make sure is at end
            }
        }

        self.equivalence = (0..n_live_cosets).collect();
        self.table.truncate(n_live_cosets);
    }

    /// swap two cosets in the table, should only be called on a compressed table
    fn switch(&mut self, coset_a: usize, coset_b: usize) {
        for g in 0..self.n_gen {
            let tmp = self.table[coset_a][g];
            self.table[coset_a][g] = self.table[coset_b][g];
            self.table[coset_b][g] = tmp;

            for coset in 0..self.table.len() {
                if self.table[coset][g] == Some(coset_a) {
                    self.table[coset][g] = Some(coset_b);
                } else if self.table[coset][g] == Some(coset_b) {
                    self.table[coset][g] = Some(coset_a);
                }
            }
        }
    }

    /// rearrange the cosets in the table into a standard order. this should only be called on a
    /// solved table
    fn standerdize(&mut self) {
        // this is 2 in the book because they count cosets from 1,
        // but we are counting from 0
        let mut expected_coset = 1;
        for coset in 0..self.table.len() {
            for g in 0..self.n_gen {
                let next_coset = self.table[coset][g].expect("todo");
                if next_coset >= expected_coset {
                    if next_coset > expected_coset {
                        self.switch(next_coset, expected_coset);
                    }
                    expected_coset += 1;
                    if expected_coset == self.table.len() - 1 {
                        return;
                    }
                }
            }
        }
    }

    /// solves the coset table
    pub fn solve(&mut self, max_iter: usize) {
        if !self.is_solved {
            self.coset_numerator_r(max_iter);
            // self.compress();
            // self.standerdize();
            self.is_solved = true;
        }
    }

    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn print_table(&self) {
        println!("{:?}", self.table);
    }
}
