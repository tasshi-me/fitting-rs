use thiserror::Error;

//#[derive(Copy, Clone, Eq, PartialEq, Debug, Fail)]
#[derive(Error, Debug, Eq, PartialEq)]
pub enum Error {
    #[error("None error")]
    Optional,
    #[error("Line Algebra error: Equations have no solutions")]
    LinalgSolveNoSolutions,
    #[error("Line Algebra error: Equations have infinite solutions")]
    LinalgSolveInfSolutions,
    #[error("Fitting error")]
    Fitting,
    #[error("Fitting error: Given y_vec contains a negative value")]
    FittingGivenYVecContainsNegativeValue,
}
