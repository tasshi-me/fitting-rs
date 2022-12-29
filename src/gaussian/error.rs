use crate::linalg::LinalgError;
use thiserror::Error;

#[derive(Error, Debug, Eq, PartialEq)]
pub enum GaussianError {
    /// Given y_vec contains a negative value
    #[error("Given y_vec contains a negative value")]
    GivenYVecContainsNegativeValue,
    /// Given y_vec contains a negative value
    #[error("Given x_vec has no element")]
    /// Given x_vec has no element
    GivenXVecHasNoElement,
    /// Error from [`crate::linalg::LinalgError`]
    #[error("Linalg error: {0:?}")]
    Linalg(#[from] LinalgError),
}
