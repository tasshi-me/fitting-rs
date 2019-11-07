use crate::linalg::solve;
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
/// use fitting::gaussian::{fit, val};
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x_vec: Vec<f64> = (1..10).map(|x| x as f64).collect();
/// let y_vec: Vec<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
/// let estimated = fit(x_vec, y_vec);
/// ```
///
/// # References
/// [1] [E. Pastuchov ́a and M. Z ́akopˇcan, ”Comparison of Algorithms for Fitting a Gaussian Function used in Testing Smart Sensors”, Journal of Electrical Engineering, vol. 66, no. 3, pp. 178-181, 2015.](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors)
pub fn fit(x_vec: Vec<f64>, y_vec: Vec<f64>) -> (f64, f64, f64) {
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
/// use fitting::gaussian::val;
///
/// let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
/// let x_vec: Vec<f64> = (1..10).map(|x| x as f64).collect();
/// let y_vec: Vec<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
/// ```
pub fn val(x: f64, mu: f64, sigma: f64, a: f64) -> f64 {
    a * (-(x - mu).powi(2) / (2. * sigma.powi(2))).exp()
}

fn caruanas(x_vec: Vec<f64>, y_vec: Vec<f64>) -> (f64, f64, f64) {
    // let x_vec: Vec<f64> = x_vec
    //     .iter()
    //     .zip(y_vec.iter())
    //     .filter_map(|(x, y)| if *y > 0. { Some(*x) } else { None })
    //     .collect();
    // let y_vec: Vec<f64> = y_vec
    //     .iter()
    //     .filter_map(|y| if *y > 0. { Some(*y) } else { None })
    //     .collect();
    // println!("x: {:?}", x_vec);
    // println!("y: {:?}", y_vec);
    let len_x_vec = x_vec.len() as f64;
    let sum_x = x_vec.iter().sum();
    let sum_x_pow2 = x_vec.iter().map(|x| x.powi(2)).sum();
    let sum_x_pow3 = x_vec.iter().map(|x| x.powi(3)).sum();
    let sum_x_pow4 = x_vec.iter().map(|x| x.powi(4)).sum();
    let a = vec![
        vec![len_x_vec, sum_x, sum_x_pow2],
        vec![sum_x, sum_x_pow2, sum_x_pow3],
        vec![sum_x_pow2, sum_x_pow3, sum_x_pow4],
    ];

    let sum_log_y = y_vec
        .iter()
        .map(|y| {
            // println!("{}", y.ln());
            y.ln()
        })
        .sum();
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
    let b = vec![sum_log_y, sum_x_log_y, sum_x_pow2_log_y];
    // println!("a: {:?}", a);
    // println!("b: {:?}", b);
    let ans_x = solve(a, b).unwrap();
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);
    // println!("a: {}, b: {}, c: {}", a, b, c);

    let mu = -b / (2. * c);
    let sigma = (-1. / (2. * c)).sqrt();
    let a = (a - (b.powi(2) / (4. * c))).exp();

    return (mu, sigma, a);
}

fn guos(x_vec: Vec<f64>, y_vec: Vec<f64>) -> (f64, f64, f64) {
    // let x_vec: Vec<f64> = x_vec
    //     .iter()
    //     .zip(y_vec.iter())
    //     .filter_map(|(x, y)| if *y > 0. { Some(*x) } else { None })
    //     .collect();
    // let y_vec: Vec<f64> = y_vec
    //     .iter()
    //     .filter_map(|y| if *y > 0. { Some(*y) } else { None })
    //     .collect();
    // println!("x: {:?}", x_vec);
    // println!("y: {:?}", y_vec);
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

    let a = vec![
        vec![sum_y_pow2, sum_x_y_pow2, sum_x_pow2_y_pow2],
        vec![sum_x_y_pow2, sum_x_pow2_y_pow2, sum_x_pow3_y_pow2],
        vec![sum_x_pow2_y_pow2, sum_x_pow3_y_pow2, sum_x_pow4_y_pow2],
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
    let b = vec![
        sum_y_pow2_log_y,
        sum_x_y_pow2_log_y,
        sum_x_pow2_y_pow2_log_y,
    ];
    // println!("a: {:?}", a);
    // println!("b: {:?}", b);
    let ans_x = solve(a, b).unwrap();
    let (a, b, c) = (ans_x[0], ans_x[1], ans_x[2]);
    // println!("a: {}, b: {}, c: {}", a, b, c);

    let mu = -b / (2. * c);
    let sigma = (-1. / (2. * c)).sqrt();
    let a = (a - (b.powi(2) / (4. * c))).exp();

    return (mu, sigma, a);
}

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
    let x_vec: Vec<f64> = (1..10).map(|x| x as f64).collect();
    let y_vec: Vec<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
    //assert_eq!(y_vec, a);
}

#[test]
fn gaussian_fit_caruanas() {
    let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
    let x_vec: Vec<f64> = (1..10).map(|x| x as f64).collect();
    let y_vec: Vec<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
    let estimated = caruanas(x_vec, y_vec);
    println!("{:?}, {:?}, {:?}", mu, sigma, a);
    println!("{:?}", estimated);
}

#[test]
fn gaussian_fit_guos() {
    let (mu, sigma, a): (f64, f64, f64) = (5., 3., 1.);
    let x_vec: Vec<f64> = (1..10).map(|x| x as f64).collect();
    let y_vec: Vec<f64> = x_vec.iter().map(|x| val(*x, mu, sigma, a)).collect();
    let estimated = guos(x_vec, y_vec);
    println!("{:?}, {:?}, {:?}", mu, sigma, a);
    println!("{:?}", estimated);
}
