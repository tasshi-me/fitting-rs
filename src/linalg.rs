use crate::error::Error;
use ndarray::{Array1, Array2};

/// Solves a system of linear equations.
///
/// This function implements the Gaussian elimination.
/// # Examples
/// Solves `a * x = b`.
///
/// ```
/// use approx::assert_abs_diff_eq;
/// use fitting::linalg::solve;
/// use ndarray::array;
///
/// let a = array![[3., 2., -1.], [2., -2., 4.], [-2., 1., -2.]];
/// let b = array![1., -2., 0.];
/// let x = solve(a, b).unwrap();
/// assert_abs_diff_eq!(x, array![1., -2., -2.], epsilon = 1e-9);
/// ```
pub fn solve(a: Array2<f64>, b: Array1<f64>) -> Result<Array1<f64>, Error> {
    let mut a = a;
    let mut b = b;
    for i in 0..(a.nrows() - 1) {
        for j in i + 1..a.nrows() {
            let coefficient = a[[j, i]] / a[[i, i]];
            // a[j] -= a[i] * coefficient;
            let a_i = a.row(i).to_owned();
            let mut view = a.row_mut(j);
            view -= &(&a_i * coefficient);
            b[j] -= b[i] * coefficient;
        }
    }

    for i in (0..a.nrows()).rev() {
        b[i] /= &a[[i, i]];
        // a[i] /= a[i][i];
        let a_i_i = a[[i, i]];
        let mut view = a.row_mut(i);
        view /= a_i_i;
        for j in 0..i {
            b[j] -= b[i] * a[[j, i]];
            a[[j, i]] = 0.;
        }
    }
    Ok(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

    #[test]
    fn linalg_solve_2x2() {
        let a = array![[3., 1.], [1., 2.]];
        let b = array![9., 8.];
        let x = solve(a, b).unwrap();
        assert_eq!(x, array![2., 3.]);
    }

    #[test]
    fn linalg_solve_3x3() {
        let a = array![[3., 2., -1.], [2., -2., 4.], [-2., 1., -2.]];
        let b = array![1., -2., 0.];
        let x = solve(a, b).unwrap();
        assert_abs_diff_eq!(x, array![1., -2., -2.], epsilon = 1e-9);
    }
}
