pub mod error;
pub mod gaussian;
pub mod linalg;

#[doc(inline)]
pub use self::gaussian::Gaussian;

// external dependencies
pub use approx;
pub use ndarray;
pub use serde;
