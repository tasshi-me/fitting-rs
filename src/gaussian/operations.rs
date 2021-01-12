use crate::error::Error;
use crate::linalg;
use crate::linalg::Float;
use ndarray::{array, Array1};

pub fn value<F: Float>(x: F, mu: F, sigma: F, a: F) -> F {
    a * (-(x - mu).powi(2) / (F::from(2).unwrap() * sigma.powi(2))).exp()
}

pub fn values<F: Float>(x_vec: Array1<F>, mu: F, sigma: F, a: F) -> Array1<F> {
    x_vec.iter().map(|x| value(*x, mu, sigma, a)).collect()
}

#[allow(dead_code)]
pub fn fitting_caruanas<F: Float>(x_vec: Array1<F>, y_vec: Array1<F>) -> Result<(F, F, F), Error> {
    let len_x_vec = F::from(x_vec.len()).ok_or(Error::Optional)?;
    let sum_x = x_vec.sum();
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
        .map(|(y, x)| y.ln() * *x)
        .sum();
    let sum_x_pow2_log_y = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| y.ln() * x.powi(2))
        .sum();
    let b = array![sum_log_y, sum_x_log_y, sum_x_pow2_log_y];

    let ans_x = linalg::solve(a, b)?;
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);

    let mu = -b / (F::from(2).unwrap() * c);
    let sigma = (-F::one() / (F::from(2).unwrap() * c)).sqrt();
    let a = (a - (b.powi(2) / (F::from(4).unwrap() * c))).exp();

    Ok((mu, sigma, a))
}

pub fn fitting_guos<F: Float>(x_vec: Array1<F>, y_vec: Array1<F>) -> Result<(F, F, F), Error> {
    let sum_y_pow2: F = y_vec.iter().map(|y| y.powi(2)).sum();
    let sum_x_y_pow2 = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| *x * y.powi(2))
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

    let sum_y_pow2_log_y: F = y_vec.iter().map(|y| y.powi(2) * y.ln()).sum();
    let sum_x_y_pow2_log_y = y_vec
        .iter()
        .zip(x_vec.iter())
        .map(|(y, x)| *x * y.powi(2) * y.ln())
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

    let ans_x = linalg::solve(a, b)?;
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);

    let mu = -b / (F::from(2).unwrap() * c);
    let sigma = ((-F::one() / (F::from(2).unwrap() * c)) as F).sqrt();
    let a = ((a - (b.powi(2) / (F::from(4).unwrap() * c))) as F).exp();

    Ok((mu, sigma, a))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::{array, Array};

    #[test]
    fn gaussian_value() {
        let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
        let x = 5.;
        let y = value(x, mu, sigma, a);
        assert_eq!(y, a);
    }

    #[test]
    fn gaussian_values() {
        let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
        let x_vec: Array1<f64> = Array::range(1., 10., 1.);
        let y_vec: Array1<f64> = values(x_vec, mu, sigma, a);
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
        let y_vec: Array1<f64> = values(x_vec.clone(), mu, sigma, a);
        let estimated = fitting_caruanas(x_vec, y_vec).unwrap();
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
        let y_vec: Array1<f64> = values(x_vec.clone(), mu, sigma, a);
        let estimated = fitting_guos(x_vec, y_vec).unwrap();
        assert_abs_diff_eq!(
            &array![estimated.0, estimated.1, estimated.2],
            &array![mu, sigma, a],
            epsilon = 1e-9
        );
    }
}
