use crate::linalg::solve;
use ndarray::{array, Array1};
use std::f64;

/// Estimates the parameters of gaussian function for generic data.
/// The return value is `(mu, sigma, a)`
///
/// This function implements the [Guos Algorithm](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors).
///
/// # Examples
/// Estimates the parameters for sample data.
///
/// ```
/// use approx::assert_abs_diff_eq;
/// use fitting::gaussian::{fit, val};
/// use ndarray::{array, Array, Array1};
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
/// let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
/// let estimated = fit(x_vec, y_vec);
/// assert_abs_diff_eq!(
///     &array![estimated.0, estimated.1, estimated.2],
///     &array![mu, sigma, a],
///     epsilon = 1e-9
/// );
/// ```
///
/// # References
/// [1] [E. Pastuchov ́a and M. Z ́akopˇcan, ”Comparison of Algorithms for Fitting a Gaussian Function used in Testing Smart Sensors”, Journal of Electrical Engineering, vol. 66, no. 3, pp. 178-181, 2015.](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors)
pub fn fit(x_vec: Array1<f64>, y_vec: Array1<f64>) -> (f64, f64, f64) {
    guos(x_vec, y_vec)
}

/// Returns a value of gaussian function.
///
/// # Examples
/// Returns a value of gaussian function.
///
/// ```
/// use fitting::gaussian::val;
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x = 5.;
/// let y = val(x, mu, sigma, a);
/// assert_eq!(y, a);
///
/// ```
///
/// # Examples
/// Returns the values of gaussian function.
///
/// ```
/// use approx::assert_abs_diff_eq;
/// use fitting::gaussian::val;
/// use ndarray::{array, Array, Array1};
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
/// let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
/// let expected_ans = array![
///     0.41111229050718745,
///     0.6065306597126334,
///     0.8007374029168081,
///     0.9459594689067654,
///     1.,
///     0.9459594689067654,
///     0.8007374029168081,
///     0.6065306597126334,
///     0.41111229050718745
/// ];
/// assert_abs_diff_eq!(&y_vec, &expected_ans, epsilon = 1e-9);
/// ```
pub fn val(x: f64, mu: f64, sigma: f64, a: f64) -> f64 {
    a * (-(x - mu).powi(2) / (2. * sigma.powi(2))).exp()
}

fn caruanas(x_vec: Array1<f64>, y_vec: Array1<f64>) -> (f64, f64, f64) {
    let len_x_vec = x_vec.len() as f64;
    let sum_x = x_vec.iter().sum();
    let sum_x_pow2 = x_vec.iter().map(|x| x.powi(2)).sum();
    let sum_x_pow3 = x_vec.iter().map(|x| x.powi(3)).sum();
    let sum_x_pow4 = x_vec.iter().map(|x| x.powi(4)).sum();
    let a = array![
        [len_x_vec, sum_x, sum_x_pow2],
        [sum_x, sum_x_pow2, sum_x_pow3],
        [sum_x_pow2, sum_x_pow3, sum_x_pow4],
    ];

    let sum_log_y = y_vec.iter().map(|y| y.ln()).sum();
    let sum_x_log_y = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| y.ln() * x)
        .sum();
    let sum_x_pow2_log_y = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| y.ln() * x.powi(2))
        .sum();
    let b = array![sum_log_y, sum_x_log_y, sum_x_pow2_log_y];

    let ans_x = solve(a, b).unwrap();
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);

    let mu = -b / (2. * c);
    let sigma = (-1. / (2. * c)).sqrt();
    let a = (a - (b.powi(2) / (4. * c))).exp();

    (mu, sigma, a)
}

fn guos(x_vec: Array1<f64>, y_vec: Array1<f64>) -> (f64, f64, f64) {
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

    let ans_x = solve(a, b).unwrap();
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);

    let mu = -b / (2. * c);
    let sigma = (-1. / (2. * c)).sqrt();
    let a = (a - (b.powi(2) / (4. * c))).exp();

    (mu, sigma, a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::{array, Array};

    #[test]
    fn gaussian_val() {
        let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
        let x = 5.;
        let y = val(x, mu, sigma, a);
        assert_eq!(y, a);
    }

    #[test]
    fn gaussian_vals() {
        let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
        let expected_ans = array![
            0.41111229050718745,
            0.6065306597126334,
            0.8007374029168081,
            0.9459594689067654,
            1.,
            0.9459594689067654,
            0.8007374029168081,
            0.6065306597126334,
            0.41111229050718745
        ];
        assert_abs_diff_eq!(&y_vec, &expected_ans, epsilon = 1e-9);
    }

    #[test]
    fn gaussian_fit_caruanas() {
        let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
        let estimated = caruanas(x_vec, y_vec);
        assert_abs_diff_eq!(
            &array![estimated.0, estimated.1, estimated.2],
            &array![mu, sigma, a],
            epsilon = 1e-9
        );
    }

    #[test]
    fn gaussian_fit_guos() {
        let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
        let estimated = guos(x_vec, y_vec);
        assert_abs_diff_eq!(
            &array![estimated.0, estimated.1, estimated.2],
            &array![mu, sigma, a],
            epsilon = 1e-9
        );
    }
}
