use crate::error::Error;
use crate::linalg::solve;
use ndarray::{array, Array1};
use std::f64;

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
}
