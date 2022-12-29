mod deprecated;
mod error;
mod gaussian;
mod operations;

#[doc(inline)]
pub use self::deprecated::*;

#[doc(inline)]
pub use self::gaussian::Gaussian;

#[doc(inline)]
pub use self::error::GaussianError;
