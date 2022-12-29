//! <br>
//!
//! This library provides [`fitting::Gaussian`][crate::Gaussian], which represents gaussian function
//!
//! <br>
//!
//! # Details
//!
//! Use [`Gaussian::new`] to get struct that represents gaussian function.
//!
//! ```
//! use fitting::Gaussian;
//!
//! let gaussian = Gaussian::new(5., 3., 1.);
//! let x = 5.;
//! let y = gaussian.value(x);
//! assert_eq!(&y, gaussian.a());
//!
//! ```
//!
//! Use [`Gaussian::value`] and [`Gaussian::values`] to get value(s) of the function.
//!
//! ```
//! use fitting::Gaussian;
//!
//! let gaussian = Gaussian::new(5., 3., 1.);
//! let x = 5.;
//! let y = gaussian.value(x);
//! assert_eq!(&y, gaussian.a());
//!
//! ```
//!
//! Use [`Gaussian::fit`] to fitting arrays to the gaussian function.
//!
//! ```
//! use fitting::approx::assert_abs_diff_eq;
//! use fitting::Gaussian;
//! use fitting::ndarray::{array, Array, Array1};
//!
//! let gaussian = Gaussian::new(5., 3., 1.);
//! let x_vec: Array1<f64> = Array::range(1., 10., 1.);
//! let y_vec: Array1<f64> = gaussian.values(x_vec.clone());
//! let estimated = Gaussian::fit(x_vec, y_vec).unwrap();
//! assert_abs_diff_eq!(gaussian, estimated, epsilon = 1e-9);
//! ```
//!

pub mod gaussian;
pub mod linalg;

#[doc(inline)]
pub use self::gaussian::Gaussian;

// external dependencies
pub use approx;
pub use ndarray;
pub use serde;
