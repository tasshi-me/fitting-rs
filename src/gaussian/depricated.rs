use crate::error::Error;
use crate::gaussian::operations;
use crate::linalg::Float;
use ndarray::Array1;

/// Returns a value of gaussian function.
///
/// # Examples
/// Returns a value of gaussian function.
///
/// ```
/// use fitting::gaussian;
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x = 5.;
/// let y = gaussian::val(x, mu, sigma, a);
/// assert_eq!(y, a);
///
/// ```
///
/// # Examples
/// Returns the values of gaussian function.
///
/// ```
/// use fitting::approx::assert_abs_diff_eq;
/// use fitting::gaussian;
/// use fitting::ndarray::{array, Array, Array1};
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
/// let y_vec: Array1<f64> = x_vec.iter().map(|x| gaussian::val(*x, mu, sigma, a)).collect();
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
#[deprecated(
    since = "0.3.0",
    note = "Please use the Gaussian::val function instead"
)]
pub fn val<F: Float>(x: F, mu: F, sigma: F, a: F) -> F {
    operations::value(x, mu, sigma, a)
}

/// Estimates the parameters of gaussian function for generic data.
/// The return value is `(mu, sigma, a)`
///
/// This function implements the [Guos Algorithm](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors).
///
/// # Examples
/// Estimates the parameters for sample data.
///
/// ```
/// use fitting::approx::assert_abs_diff_eq;
/// use fitting::gaussian;
/// use fitting::ndarray::{array, Array, Array1};
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
/// let y_vec: Array1<f64> = x_vec.iter().map(|x| gaussian::val(*x, mu, sigma, a)).collect();
/// let estimated = gaussian::fit(x_vec, y_vec).unwrap();
/// assert_abs_diff_eq!(
///     &array![estimated.0, estimated.1, estimated.2],
///     &array![mu, sigma, a],
///     epsilon = 1e-9
/// );
/// ```
///
/// # References
/// [1] [E. Pastuchov ́a and M. Z ́akopˇcan, ”Comparison of Algorithms for Fitting a Gaussian Function used in Testing Smart Sensors”, Journal of Electrical Engineering, vol. 66, no. 3, pp. 178-181, 2015.](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors)
#[deprecated(
    since = "0.3.0",
    note = "Please use the Gaussian::fit function instead"
)]
pub fn fit<F: Float>(x_vec: Array1<F>, y_vec: Array1<F>) -> Result<(F, F, F), Error> {
    operations::fitting_guos(x_vec, y_vec)
}
