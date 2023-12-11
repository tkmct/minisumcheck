use std::collections::HashMap;

use ark_bls12_381::Fq;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    univariate::DensePolynomial,
    DenseMVPolynomial, DenseUVPolynomial, Polynomial,
};

// This is a sumcheck protocol step by step
// Domain of g is boolean: H = {0,1}
fn main() {
    // Statement declaration
    // Create a multivariate polynomial in 3 variables, with 4 terms:
    // g(x0, x1, x2) = 2*x_0^3 + x_0*x_2 + x_1*x_2 + 5
    let g = SparsePolynomial::from_coefficients_vec(
        3,
        vec![
            (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
            (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
            (Fq::from(5), SparseTerm::new(vec![])),
        ],
    );
    // calculate all the sum.
    let mut sum = Fq::from(0);
    for i in 0..u32::pow(2, g.num_vars as u32) {
        let point: Vec<Fq> = (0..g.num_vars)
            .map(|n| ((i >> n) & 1) as i8)
            .map(Fq::from)
            .collect();

        sum += g.evaluate(&point);
    }

    println!("The sum is {}", sum);

    // Now, protocol starts here.
    // 1. Prover calculates a univariate polynomial q0(x) and sends it to verifier.
    // q0(x) = Σ(w1,w2) g(x0, w1, w2)

    // Stores coefficients for each index
    // get degree of the x0
    let i = 0;

    let mut degree = 0;
    for t in g.terms.iter() {
        for v in t.1.iter() {
            if v.0 == i && v.1 > degree {
                degree = v.1;
            }
        }
    }
    let num_vars = g.num_vars - 1;
    let num_permutation = i32::pow(2, num_vars as u32);
    println!("num permutation:{}", num_permutation);

    let mut coeffs: Vec<Fq> = vec![Fq::from(0); degree + 1];

    for term in g.terms {
        let (c, vars) = term;
        let mut has_0 = 0;

        let idx = if let Some(x0) = vars.iter().find(|v| v.0 == i) {
            has_0 = 1;
            x0.1
        } else {
            0
        };

        let term_len = vars.len() - has_0;
        coeffs[idx] += Fq::from(i32::pow(2, (num_vars - term_len) as u32)) * c;
    }

    let q0 = DensePolynomial::from_coefficients_vec(coeffs);
    println!("Polynomial: {:?}", q0);
}
