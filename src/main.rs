use rand::thread_rng;
use std::collections::HashMap;

use ark_bls12_381::Fq;
use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    univariate::DensePolynomial,
    DenseMVPolynomial, DenseUVPolynomial, Polynomial,
};
use ark_std::UniformRand;

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
    println!("Statement for round0: ");
    println!("g0: {:?}", g);
    // calculate all the sum.
    let mut u0 = Fq::from(0);
    for i in 0..u32::pow(2, g.num_vars as u32) {
        let point: Vec<Fq> = (0..g.num_vars)
            .map(|n| ((i >> n) & 1) as i8)
            .map(Fq::from)
            .collect();

        u0 += g.evaluate(&point);
    }

    println!("u0: {}", u0);

    // Now, protocol starts here.
    //
    // --- Round 0 start---
    //
    // 1. Prover calculates a univariate polynomial q0(x) and sends it to verifier.
    // q0(x) = Σ(w1,w2) g(x0, w1, w2)

    println!("---Round0 start---");
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
    let mut coeffs: Vec<Fq> = vec![Fq::from(0); degree + 1];

    for term in g.terms.iter() {
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

    // 2. Prover sends q0 to Verifier
    let q0 = DensePolynomial::from_coefficients_vec(coeffs);
    println!("q0(x0) = {:?}", q0);

    // 3. Verifier checks if u = q0(0) + q0(1)
    assert!(q0.evaluate(&Fq::from(0)) + q0.evaluate(&Fq::from(1)) == u0);
    println!("checks u0 == q0(0) = q0(1): Passed.");

    // 4. Sample random r0 and send to Prover
    let mut rng = thread_rng();
    // let r0 = Fq::rand(&mut rng);
    // For testing purpose
    let r0 = Fq::from(2);

    println!("r0: {:?}", r0);

    // 5. Prover calculate g' = g(r0, x1, x2)
    let g1 = {
        // Store coefficients for each term.
        let mut term_reg = HashMap::<Vec<(usize, usize)>, Fq>::new();

        g.terms.iter().for_each(|term| {
            let mut c = term.0;

            if let Some(x0) = term.1.iter().find(|v| v.0 == i) {
                c *= r0.pow([x0.1 as u64]);
            }

            let new_vars = term
                .1
                .iter()
                .filter(|v| v.0 != i)
                .copied()
                .map(|(var, c)| (var - 1, c)) // decrement var number by 1
                .collect::<Vec<_>>();
            // println!("New vars for term:{:?}, vars: {:?}", term, new_vars);

            let zero = Fq::from(0);
            let coeff = term_reg.get(&new_vars).unwrap_or(&zero);
            term_reg.insert(new_vars, coeff + &c);
        });

        let terms = term_reg
            .into_iter()
            .map(|(key, val)| (val, SparseTerm::new(key)))
            .collect::<Vec<(Fq, SparseTerm)>>();

        // println!("degree: {}", g.degree() - 1);
        // println!("Terms: {:?}", terms);

        SparsePolynomial::from_coefficients_vec(g.num_vars() - 1, terms)
    };

    println!("---Round0 end---");
    // --- Round 0 end---

    // statement g1, u_1
    println!("Statement for round1");
    println!("g1: {:?}", g1);
    let u1 = q0.evaluate(&r0);
    println!("u1: {}", u1);

    // --- Round 1 start ---
    println!("---Round1 start---");
    // Prover calculates q1
    let mut degree = 0;
    for t in g1.terms.iter() {
        for v in t.1.iter() {
            if v.0 == i && v.1 > degree {
                degree = v.1;
            }
        }
    }
    let num_vars = g1.num_vars - 1;
    let mut coeffs: Vec<Fq> = vec![Fq::from(0); degree + 1];

    for term in g1.terms.iter() {
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

    // 2. Prover sends q1 to Verifier
    let q1 = DensePolynomial::from_coefficients_vec(coeffs);
    println!("q1: {:?}", q1);

    // 3. Verifier checks if u1 = q1(0) + q1(1)
    assert!(q1.evaluate(&Fq::from(0)) + q1.evaluate(&Fq::from(1)) == u1);
    println!("checks u1 == q1(0) = q1(1): Passed.");

    // 4. Verifier pick random challenge
    let r1 = Fq::from(7);

    println!("r1: {:?}", r1);

    // 5. Prover calculate g2 = g(r0, r1, x2)
    let g2 = {
        // Store coefficients for each term.
        let mut term_reg = HashMap::<Vec<(usize, usize)>, Fq>::new();

        g1.terms.iter().for_each(|term| {
            let mut c = term.0;

            if let Some(x0) = term.1.iter().find(|v| v.0 == i) {
                c *= r1.pow([x0.1 as u64]);
            }

            let new_vars = term
                .1
                .iter()
                .filter(|v| v.0 != i)
                .copied()
                .map(|(var, c)| (var - 1, c)) // decrement var number by 1
                .collect::<Vec<_>>();

            let zero = Fq::from(0);
            let coeff = term_reg.get(&new_vars).unwrap_or(&zero);
            term_reg.insert(new_vars, coeff + &c);
        });

        let terms = term_reg
            .into_iter()
            .map(|(key, val)| (val, SparseTerm::new(key)))
            .collect::<Vec<(Fq, SparseTerm)>>();

        SparsePolynomial::from_coefficients_vec(g1.num_vars() - 1, terms)
    };

    println!("---Round1 end---");
    // --- Round 1 end ---
    // statement g2, u_2
    println!("Statement for round2");
    println!("g2: {:?}", g2);
    let u2 = q1.evaluate(&r1);
    println!("u2: {}", u2);

    // --- Round2 start---
    println!("---Round2 start---");
    // Prover calculates q1
    let mut degree = 0;
    for t in g2.terms.iter() {
        for v in t.1.iter() {
            if v.0 == i && v.1 > degree {
                degree = v.1;
            }
        }
    }
    let num_vars = g2.num_vars - 1;
    let mut coeffs: Vec<Fq> = vec![Fq::from(0); degree + 1];

    for term in g2.terms.iter() {
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

    // 2. Prover sends q2 to Verifier
    let q2 = DensePolynomial::from_coefficients_vec(coeffs);
    println!("q2: {:?}", q2);

    // 3. Verifier checks if u = q2(0) + q2(1)
    assert!(q2.evaluate(&Fq::from(0)) + q2.evaluate(&Fq::from(1)) == u2);
    println!("checks u2 == q2(0) = q2(1): Passed.");

    // 4. Verifier pick random challenge
    let r2 = Fq::from(3);

    println!("r2: {:?}", r2);

    // 5. Prover calculate g3 = g(r0, r1, r2)
    let g3 = {
        // Store coefficients for each term.
        let mut term_reg = HashMap::<Vec<(usize, usize)>, Fq>::new();

        g2.terms.iter().for_each(|term| {
            let mut c = term.0;

            if let Some(x0) = term.1.iter().find(|v| v.0 == i) {
                c *= r2.pow([x0.1 as u64]);
            }

            let new_vars = term
                .1
                .iter()
                .filter(|v| v.0 != i)
                .copied()
                .map(|(var, c)| (var - 1, c)) // decrement var number by 1
                .collect::<Vec<_>>();

            let zero = Fq::from(0);
            let coeff = term_reg.get(&new_vars).unwrap_or(&zero);
            term_reg.insert(new_vars, coeff + &c);
        });

        let terms = term_reg
            .into_iter()
            .map(|(key, val)| (val, SparseTerm::new(key)))
            .collect::<Vec<(Fq, SparseTerm)>>();

        SparsePolynomial::from_coefficients_vec(g2.num_vars() - 1, terms)
    };

    println!("---Round2 end---");

    println!("g3: {:?}", g3);
    let u3 = q2.evaluate(&r2);
    println!("Verifier evaluates g with (r0,r1,r2)");
    let eval = g.evaluate(&vec![r0, r1, r2]);
    println!("g(r0,r1,r2) = {:?}", eval);
    assert!(eval == u3);
    println!("Checking g(r0,r1,r2) == u3: Passed.");
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fq;
    use ark_ff::Field;

    #[test]
    fn test_pow() {
        let x = Fq::from(3);
        let p = 2;

        let res = x.pow([p]);
        assert_eq!(res, Fq::from(9));
    }
}
