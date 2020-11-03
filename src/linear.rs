use crate::error::Error;
use crate::linalg::solve;
use ndarray::{array, Array1};
use std::f64;

/// Estimates the parameters of linear function for generic data.
/// The return value is `(a, b)`
///
/// This function implements the [Linear Least Squares](http://mathworld.wolfram.com/LeastSquaresFitting.html).
///
/// # Examples
/// Estimates the parameters for sample data.
///
/// ```
/// use approx::assert_abs_diff_eq;
/// use fitting::linear::{fit, val};
/// use ndarray::{array, Array, Array1};
///
/// let (a, b): (f64, f64) = (5., 3.);
/// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
/// let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, a, b)).collect();
/// let estimated = fit(x_vec, y_vec).unwrap();
/// assert_abs_diff_eq!(
///     &array![estimated.0, estimated.1],
///     &array![a, b],
///     epsilon = 1e-9
/// );
/// ```
///
/// # References
/// [1] [E. Pastuchov ́a and M. Z ́akopˇcan, ”Comparison of Algorithms for Fitting a Gaussian Function used in Testing Smart Sensors”, Journal of Electrical Engineering, vol. 66, no. 3, pp. 178-181, 2015.](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors)
pub fn fit(x_vec: Array1<f64>, y_vec: Array1<f64>) -> Result<(f64, f64), Error> {
    Ok(linear_least_squares(x_vec, y_vec)?)
}

/// Returns a value of linear function.
///
/// # Examples
/// Returns a value of linear function.
///
/// ```
/// use fitting::linear::val;
///
/// // y = 1.2 * x + 3
/// let (a, b): (f64, f64) = (1.2, 3.);
/// let x = 5.;
/// let y = val(x, a, b);
/// assert_eq!(y, 9.);
///
/// ```
///
/// # Examples
/// Returns the values of linear function.
///
/// ```
/// use approx::assert_abs_diff_eq;
/// use fitting::linear::val;
/// use ndarray::{array, Array, Array1};
///
/// // y = 1.2 * x + 3
/// let (a, b): (f64, f64) = (1.2, 3.);
/// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
/// let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, a, b)).collect();
/// let expected_ans = array![4.2, 5.4, 6.6, 7.8, 9., 10.2, 11.4, 12.6, 13.799999999999999];
/// assert_abs_diff_eq!(&y_vec, &expected_ans, epsilon = 1e-9);
///
/// ```
pub fn val(x: f64, a: f64, b: f64) -> f64 {
    a * x + b
}

fn linear_least_squares(x_vec: Array1<f64>, y_vec: Array1<f64>) -> Result<(f64, f64), Error> {
    let sum_y_pow2 = y_vec.iter().map(|y| y.powi(2)).sum();
    let sum_x_y_pow2 = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| x * y.powi(2))
        .sum();
    let sum_x_pow2_y_pow2 = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| x.powi(2) * y.powi(2))
        .sum();
    let sum_x_pow3_y_pow2 = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| x.powi(3) * y.powi(2))
        .sum();
    let sum_x_pow4_y_pow2 = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| x.powi(4) * y.powi(2))
        .sum();

    let a = array![
        [sum_y_pow2, sum_x_y_pow2, sum_x_pow2_y_pow2],
        [sum_x_y_pow2, sum_x_pow2_y_pow2, sum_x_pow3_y_pow2],
        [sum_x_pow2_y_pow2, sum_x_pow3_y_pow2, sum_x_pow4_y_pow2],
    ];

    let sum_y_pow2_log_y = y_vec.iter().map(|y| y.powi(2) * y.ln()).sum();
    let sum_x_y_pow2_log_y = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| x * y.powi(2) * y.ln())
        .sum();
    let sum_x_pow2_y_pow2_log_y = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| x.powi(2) * y.powi(2) * y.ln())
        .sum();
    let b = array![
        sum_y_pow2_log_y,
        sum_x_y_pow2_log_y,
        sum_x_pow2_y_pow2_log_y,
    ];

    let ans_x = solve(a, b)?;
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);

    let mu = -b / (2. * c);
    let sigma = (-1. / (2. * c)).sqrt();
    let _a = (a - (b.powi(2) / (4. * c))).exp();

    Ok((mu, sigma))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::{array, Array};

    #[test]
    fn linear_val_positive_slope() {
        // y = 1.2 * x + 3
        let (a, b): (f64, f64) = (1.2, 3.);
        let x = 5.;
        let y = val(x, a, b);
        assert_eq!(y, 9.);
    }

    #[test]
    fn linear_val_negative_slope() {
        // y = -1.2 * x + 3
        let (a, b): (f64, f64) = (-1.2, 3.);
        let x = 5.;
        let y = val(x, a, b);
        assert_eq!(y, -3.);
    }

    #[test]
    fn linear_val_zero_slope() {
        // y = 0 * x + 3
        let (a, b): (f64, f64) = (0., 3.);
        let x = 5.;
        let y = val(x, a, b);
        assert_eq!(y, 3.);
    }

    #[test]
    fn linear_vals_positive_slope() {
        // y = 1.2 * x + 3
        let (a, b): (f64, f64) = (1.2, 3.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, a, b)).collect();
        let expected_ans = array![4.2, 5.4, 6.6, 7.8, 9., 10.2, 11.4, 12.6, 13.799999999999999];
        assert_abs_diff_eq!(&y_vec, &expected_ans, epsilon = 1e-9);
    }

    #[test]
    fn linear_vals_negative_slope() {
        // y = -1.2 * x + 3
        let (a, b): (f64, f64) = (-1.2, 3.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, a, b)).collect();
        let expected_ans = array![
            1.8,
            0.6000000000000001,
            -0.5999999999999996,
            -1.7999999999999998,
            -3.,
            -4.199999999999999,
            -5.4,
            -6.6,
            -7.799999999999999
        ];
        assert_abs_diff_eq!(&y_vec, &expected_ans, epsilon = 1e-9);
    }

    #[test]
    fn linear_vals_zero_slope() {
        // y = 0 * x + 3
        let (a, b): (f64, f64) = (0., 3.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, a, b)).collect();
        let expected_ans = array![3., 3., 3., 3., 3., 3., 3., 3., 3.];
        assert_abs_diff_eq!(&y_vec, &expected_ans, epsilon = 1e-9);
    }

    #[test]
    fn linear_fit_linear_least_squares() {
        let (a, b): (f64, f64) = (1.2, 3.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, a, b)).collect();
        let estimated = linear_least_squares(x_vec, y_vec).unwrap();
        assert_abs_diff_eq!(
            &array![estimated.0, estimated.1],
            &array![a, b],
            epsilon = 1e-9
        );
    }
}
