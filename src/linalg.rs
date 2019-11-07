/// Solves a system of linear equations.
///
/// This function implements the Gaussian elimination.
/// # Examples
/// Solves `a * x = b`.
///
/// ```
/// use fitting::linalg::solve;
///
/// let a = vec![vec![3., 1.], vec![1., 2.]];
/// let b = vec![9., 8.];
///
/// let x = solve(a, b).unwrap();
/// assert_eq!(x, vec![2., 3.]);
/// ```
pub fn solve(a: Vec<Vec<f64>>, b: Vec<f64>) -> Option<Vec<f64>> {
    let mut a = a;
    let mut b = b;
    for i in 0..(a.len() - 1) {
        for j in i + 1..a.len() {
            // println!("for: {},{}", i, j);
            let coefficient = a[j][i] / a[i][i];
            // a[j] -= a[i] * coefficient;
            a[j] = a[j]
                .iter()
                .zip(a[i].iter())
                .map(|(aj, ai)| aj - ai * coefficient)
                .collect();
            b[j] -= b[i] * coefficient;
        }
    }

    for i in (0..a.len()).rev() {
        b[i] /= a[i][i];
        // a[i] /= a[i][i];
        a[i] = a[i].iter().map(|ai| ai / a[i][i]).collect();
        for j in 0..i {
            // println!("for: {},{}", i, j);
            b[j] -= b[i] * a[j][i];
            a[j][i] = 0.;
        }
    }
    // println!("{:?}", b);
    Some(b)
}

#[test]
fn linalg_solve_2x2() {
    let a = vec![vec![3., 1.], vec![1., 2.]];
    let b = vec![9., 8.];
    let x = solve(a, b).unwrap();
    assert_eq!(x, vec![2., 3.]);
}

// #[test]
// fn linalg_solve_3x3() {
//     let a = vec![vec![3., 2., 1.], vec![2., 2., 4.], vec![-2., 1., -2.]];
//     let b = vec![1., -2., 0.];

//     let x = solve(a, b).unwrap();
//     println!("{:?}", x);
//     assert_eq!(x, vec![0.69230769, -0.15384615, -0.76923077]);
// }
