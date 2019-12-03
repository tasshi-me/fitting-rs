use crate::error::{Error, ErrorKind};
use approx::{abs_diff_eq, abs_diff_ne};
use ndarray::{s, Array1, Array2, Axis};

pub struct LUDecompositioned {
    pub lu: Array2<f64>,
    pub permutation: Array2<f64>,
}

pub enum LUDecompositionMethod {
    RightLooking,
    LeftLooking,
    CroutsMethod,
    BlockGaussian,
}

impl LUDecompositioned {
    pub fn new(a: Array2<f64>) -> Result<LUDecompositioned, Error> {
        let (lu, permutation) = lu_decomposition(a, LUDecompositionMethod::CroutsMethod)?;
        Ok(LUDecompositioned {
            lu: lu,
            permutation: permutation,
        })
    }

    pub fn new_with_method(
        a: Array2<f64>,
        method: LUDecompositionMethod,
    ) -> Result<LUDecompositioned, Error> {
        let (lu, permutation) = lu_decomposition(a, method)?;
        Ok(LUDecompositioned {
            lu: lu,
            permutation: permutation,
        })
    }

    // pub fn solve(b: Array1<f64>) -> Result<Array1<f64>, Error>){
    //     Ok()
    // }
}

fn lu_decomposition(
    a: Array2<f64>,
    method: LUDecompositionMethod,
) -> Result<(Array2<f64>, Array2<f64>), Error> {
    match method {
        LUDecompositionMethod::RightLooking => lu_decomposition_right_looking(a),
        LUDecompositionMethod::LeftLooking => lu_decomposition_left_looking(a),
        LUDecompositionMethod::CroutsMethod => lu_decomposition_crouts_method(a),
        LUDecompositionMethod::BlockGaussian => lu_decomposition_block_gaussian(a),
    }
}

fn lu_decomposition_right_looking(a: Array2<f64>) -> Result<(Array2<f64>, Array2<f64>), Error> {
    let mut a = a;
    let mut lu = Array2::<f64>::zeros(a.dim());
    let mut p = Array2::eye(a.nrows());

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
            let mut view = a.row_mut(j);
            view -= &(&a_i * coefficient);
            b[j] -= b[i] * coefficient;
        }
    }

    // check rank of matrix
    // rank_coef: rank of coefficient matrix (given a)
    // rank_aug: rank of augmented matrix
    let mut rank_coef = a.nrows();
    for index in (0..a.nrows()).rev() {
        if a.row(index)
            .iter()
            .all(|val| abs_diff_eq!(*val, 0.) || val.is_nan())
        {
            rank_coef -= 1;
        } else {
            break;
        }
    }
    let rank_coef = rank_coef;

    let mut rank_aug = rank_coef;
    for index in ((rank_coef - 1)..a.nrows()).rev() {
        if abs_diff_ne!(b[index], 0.) && !b[index].is_nan() {
            rank_aug = index + 1;
            break;
        }
    }
    let rank_aug = rank_aug;

    // no solutions
    if rank_coef != rank_aug {
        return Err(Error::from(ErrorKind::LinalgSolveNoSolutions));
    }

    // infinite solutions
    if rank_coef != a.ncols() {
        return Err(Error::from(ErrorKind::LinalgSolveInfSolutions));
    }

    // backward substitution
    for i in (0..rank_coef).rev() {
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
    Ok(b.slice(s![0..rank_coef]).to_owned())
}

fn lu_decomposition_left_looking(a: Array2<f64>) -> Result<(Array2<f64>, Array2<f64>), Error> {
    Ok((a, Array2::zeros((1, 1))))
}
fn lu_decomposition_crouts_method(a: Array2<f64>) -> Result<(Array2<f64>, Array2<f64>), Error> {
    Ok((a, Array2::zeros((1, 1))))
}
fn lu_decomposition_block_gaussian(a: Array2<f64>) -> Result<(Array2<f64>, Array2<f64>), Error> {
    Ok((a, Array2::zeros((1, 1))))
}
