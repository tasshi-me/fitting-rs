use crate::error::Error;
use approx::{abs_diff_eq, abs_diff_ne};
use ndarray::{s, Array1, Array2, Axis};

use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use ndarray::NdFloat;
use std::iter::Sum;

pub trait Float: NdFloat + Sum + AbsDiffEq + RelativeEq + UlpsEq {}
impl<F: NdFloat + Sum + AbsDiffEq + RelativeEq + UlpsEq> Float for F {}

/// Solves a system of linear equations.
///
/// This function implements the Gaussian elimination.
/// # Examples
/// Solves `a * x = b`.
///
/// ```
/// use fitting::approx::assert_abs_diff_eq;
/// use fitting::linalg;
/// use fitting::ndarray::array;
///
/// let a = array![[3., 2., -1.], [2., -2., 4.], [-2., 1., -2.]];
/// let b = array![1., -2., 0.];
/// let x = linalg::solve(a, b).unwrap();
/// assert_abs_diff_eq!(x, array![1., -2., -2.], epsilon = 1e-9);
/// ```
pub fn solve<F: Float>(a: Array2<F>, b: Array1<F>) -> Result<Array1<F>, Error> {
    let mut a = a;
    let mut b = b;

    // forward elimination
    for i in 0..(a.nrows() - 1) {
        // partial pivoting
        let (pivot_index, _) = a.column(i).iter().enumerate().skip(i).fold(
            (i, a[[i, i]]),
            |(max_index, max), (index, val)| {
                if val.abs() > max {
                    (index, *val)
                } else {
                    (max_index, max)
                }
            },
        );
        if i != pivot_index {
            let (mut a1, mut a2) = a.view_mut().split_at(Axis(0), pivot_index);
            ndarray::Zip::from(a1.row_mut(i))
                .and(a2.row_mut(0))
                .apply(::std::mem::swap);
            b.swap(i, pivot_index);
        }
        for j in i + 1..a.nrows() {
            let coefficient = a[[j, i]] / a[[i, i]];
            // a[j] -= a[i] * coefficient;
            let a_i = a.row(i).to_owned();
            let b_i = b[i];
            let mut view = a.row_mut(j);
            view -= &(&a_i * coefficient);
            b[j] -= b_i * coefficient;
        }
    }

    // check rank of matrix
    // rank_coef: rank of coefficient matrix (given a)
    // rank_aug: rank of augmented matrix
    let mut rank_coef = a.nrows();
    for index in (0..a.nrows()).rev() {
        if a.row(index)
            .iter()
            .all(|val| abs_diff_eq!(*val, F::zero()) || val.is_nan())
        {
            rank_coef -= 1;
        } else {
            break;
        }
    }
    let rank_coef = rank_coef;

    let mut rank_aug = rank_coef;
    for index in ((rank_coef - 1)..a.nrows()).rev() {
        if abs_diff_ne!(b[index], F::zero()) && !b[index].is_nan() {
            rank_aug = index + 1;
            break;
        }
    }
    let rank_aug = rank_aug;

    // no solutions
    if rank_coef != rank_aug {
        return Err(Error::LinalgSolveNoSolutions);
    }

    // infinite solutions
    if rank_coef != a.ncols() {
        return Err(Error::LinalgSolveInfSolutions);
    }

    // backward substitution
    for i in (0..rank_coef).rev() {
        b[i] /= a[[i, i]];
        // a[i] /= a[i][i];
        let a_i_i = a[[i, i]];
        let mut view = a.row_mut(i);
        view /= a_i_i;
        for j in 0..i {
            let b_i = b[i];
            b[j] -= b_i * a[[j, i]];
            a[[j, i]] = F::zero();
        }
    }
    Ok(b.slice(s![0..rank_coef]).to_owned())
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

    #[test]
    fn linalg_solve_pivoting() {
        let a = array![[2., 4., -2.], [1., 2., 1.], [1., 3., 2.],];
        let b = array![8., 6., 9.];
        let x = solve(a, b).unwrap();
        println!("{:?}", x);
        assert_abs_diff_eq!(x, array![1., 2., 1.], epsilon = 1e-9);
    }

    #[test]
    // one solution (4x4)
    // rank_coef: 4
    // rank_aug: 4
    // coef ncols: 4
    fn linalg_solve_has_one_solution() {
        let a = array![
            [2., 1., -3., -2.],
            [2., -1., -1., 3.],
            [1., -1., -2., 2.],
            [-1., 1., 3., -2.]
        ];
        let b = array![-4., 1., -3., 5.];
        let x = solve(a, b).unwrap();
        assert_abs_diff_eq!(x, array![1., 2., 2., 1.], epsilon = 1e-9);
    }

    #[test]
    // one solution (4x3)
    // rank_coef: 3
    // rank_aug: 3
    // coef ncols: 3
    fn linalg_solve_has_one_solution_2() {
        let a = array![[2., 1., -3.], [2., -1., -1.], [1., -1., -2.], [-1., 1., 3.]];
        let b = array![-2., -2., -5., 7.];
        let x = solve(a, b).unwrap();
        assert_abs_diff_eq!(x, array![1., 2., 2.], epsilon = 1e-9);
    }

    #[test]
    // inf solutions (3x4)
    // rank_coef: 3
    // coef ncols: 4
    fn linalg_solve_has_inf_solutions() {
        let a = array![[2., 1., -3., -2.], [2., -1., -1., 3.], [1., -1., -2., 2.]];
        let b = array![4., 1., -3.];
        let err = solve(a, b).unwrap_err(); //panic
        assert_eq!(err, Error::LinalgSolveInfSolutions);
    }

    #[test]
    // inf solutions (4x4)
    // rank_coef: 2
    // coef ncols: 4
    fn linalg_solve_has_inf_solutions_2() {
        let a = array![
            [2., 1., 3., 4.],
            [2., -3., -1., -4.],
            [1., -2., -1., -3.],
            [-1., 2., 1., 3.]
        ];
        let b = array![2., -6. / 5., -1., 1.];
        let err = solve(a, b).unwrap_err(); //panic
        assert_eq!(err, Error::LinalgSolveInfSolutions);
    }

    #[test]
    // no solutions (3x2)
    // rank_coef: 2
    // rank_aug: 3
    fn linalg_solve_has_no_solutions() {
        let a = array![[-2., 3.], [4., 1.], [1., -3.],];
        let b = array![1., 5., -1.];
        let err = solve(a, b).unwrap_err(); //panic
        assert_eq!(err, Error::LinalgSolveNoSolutions);
    }

    #[test]
    // no solutions (3x3)
    // rank_coef: 2
    // rank_aug: 3
    fn linalg_solve_has_no_solutions_2() {
        let a = array![[1., 3., -2.], [-1., 2., -3.], [2., -1., 3.],];
        let b = array![2., -2., 3.];
        let err = solve(a, b).unwrap_err(); //panic
        assert_eq!(err, Error::LinalgSolveNoSolutions);
    }
}
