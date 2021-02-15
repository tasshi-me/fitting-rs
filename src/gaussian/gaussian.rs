use std::convert::TryFrom;

use crate::error::Error;
use crate::gaussian::operations;
use crate::linalg::Float;
use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use ndarray::{array, Array1};
use serde::{Deserialize, Serialize};

/// Gaussian Distribution
#[derive(Copy, Clone, Eq, PartialEq, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub struct Gaussian<F: Float> {
    mu: F,
    sigma: F,
    a: F,
}

impl<F: Float> Gaussian<F> {
    /// Create a new `Gaussian` with given parameters.
    pub fn new(mu: F, sigma: F, a: F) -> Gaussian<F> {
        Gaussian { mu, sigma, a }
    }

    /// Return a reference to `mu`.
    pub fn mu(&self) -> &F {
        &self.mu
    }

    /// Return a mutable reference to `mu`.
    pub fn mu_mut(&mut self) -> &mut F {
        &mut self.mu
    }

    /// Return a reference to `sigma`.
    pub fn sigma(&self) -> &F {
        &self.sigma
    }

    /// Return a mutable reference to `sigma`.
    pub fn sigma_mut(&mut self) -> &mut F {
        &mut self.sigma
    }

    /// Return a reference to `a`.
    pub fn a(&self) -> &F {
        &self.a
    }

    /// Return a mutable reference to `a`.
    pub fn a_mut(&mut self) -> &mut F {
        &mut self.a
    }

    /// Return a reference to `self` as a tuple `(mu, sigma, a)`.
    pub fn as_tuple(&self) -> (&F, &F, &F) {
        (&self.mu, &self.sigma, &self.a)
    }

    /// Return a mutable reference to `self` as a tuple `(mu, sigma, a)`.
    pub fn as_mut_tuple(&mut self) -> (&mut F, &mut F, &mut F) {
        (&mut self.mu, &mut self.sigma, &mut self.a)
    }

    /// Returns a copy of `self` as a new tuple `(mu, sigma, a)`.
    pub fn to_tuple(&self) -> (F, F, F) {
        (self.mu, self.sigma, self.a)
    }

    /// Consume `self` and returns a tuple `(mu, sigma, a)`.
    pub fn into_tuple(self) -> (F, F, F) {
        (self.mu, self.sigma, self.a)
    }

    /// Returns a value of gaussian function.
    ///
    /// # Examples
    /// Returns a value of gaussian function.
    ///
    /// ```
    /// use fitting::Gaussian;
    ///
    /// let gaussian = Gaussian::new(5., 3., 1.);
    /// let x = 5.;
    /// let y = gaussian.value(x);
    /// assert_eq!(&y, gaussian.a());
    ///
    /// ```
    pub fn value(&self, x: F) -> F {
        operations::value(x, self.mu, self.sigma, self.a)
    }

    /// Returns the values of gaussian function.
    ///
    /// # Examples
    /// Returns the values of gaussian function.
    ///
    /// ```
    /// use fitting::approx::assert_abs_diff_eq;
    /// use fitting::Gaussian;
    /// use fitting::ndarray::{array, Array, Array1};
    ///
    /// let gaussian = Gaussian::new(5., 3., 1.);
    /// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
    /// let y_vec: Array1<f64> = gaussian.values(x_vec);
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
    pub fn values(&self, x_vec: Array1<F>) -> Array1<F> {
        operations::values(x_vec, self.mu, self.sigma, self.a)
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
    /// use fitting::Gaussian;
    /// use fitting::ndarray::{array, Array, Array1};
    ///
    /// let gaussian = Gaussian::new(5., 3., 1.);
    /// let x_vec: Array1<f64> = Array::range(1., 10., 1.);
    /// let y_vec: Array1<f64> = gaussian.values(x_vec.clone());
    /// let estimated = Gaussian::fit(x_vec, y_vec).unwrap();
    /// assert_abs_diff_eq!(gaussian, estimated, epsilon = 1e-9);
    /// ```
    ///
    /// # References
    /// [1] [E. Pastuchov ́a and M. Z ́akopˇcan, ”Comparison of Algorithms for Fitting a Gaussian Function used in Testing Smart Sensors”, Journal of Electrical Engineering, vol. 66, no. 3, pp. 178-181, 2015.](https://www.researchgate.net/publication/281907940_Comparison_of_Algorithms_For_Fitting_a_Gaussian_Function_Used_in_Testing_Smart_Sensors)
    pub fn fit(x_vec: Array1<F>, y_vec: Array1<F>) -> Result<Gaussian<F>, Error> {
        let (mu, sigma, a) = operations::fitting_guos(x_vec, y_vec)?;
        Ok(Gaussian::<F>::new(mu, sigma, a))
    }
}

impl<F: Float> From<(F, F, F)> for Gaussian<F> {
    fn from(tuple: (F, F, F)) -> Self {
        Gaussian::new(tuple.0, tuple.1, tuple.2)
    }
}

impl<F: Float> From<Gaussian<F>> for (F, F, F) {
    fn from(gaussian: Gaussian<F>) -> Self {
        (gaussian.mu, gaussian.sigma, gaussian.a)
    }
}

impl<F: Float> TryFrom<Array1<F>> for Gaussian<F> {
    type Error = &'static str;

    fn try_from(arr: Array1<F>) -> Result<Self, Self::Error> {
        let error = "The index out of bounds.";
        Ok(Gaussian::new(
            *arr.get(0).ok_or(error)?,
            *arr.get(1).ok_or(error)?,
            *arr.get(2).ok_or(error)?,
        ))
    }
}

impl<F: Float> From<Gaussian<F>> for Array1<F> {
    fn from(gaussian: Gaussian<F>) -> Self {
        array![gaussian.mu, gaussian.sigma, gaussian.a]
    }
}

impl<F: Float> AbsDiffEq for Gaussian<F>
where
    F::Epsilon: Copy,
{
    type Epsilon = F::Epsilon;

    fn default_epsilon() -> F::Epsilon {
        F::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: F::Epsilon) -> bool {
        F::abs_diff_eq(&self.mu, &other.mu, epsilon)
            && F::abs_diff_eq(&self.sigma, &other.sigma, epsilon)
            && F::abs_diff_eq(&self.a, &other.a, epsilon)
    }
}

impl<F: Float> RelativeEq for Gaussian<F>
where
    F::Epsilon: Copy,
{
    fn default_max_relative() -> F::Epsilon {
        F::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: F::Epsilon, max_relative: F::Epsilon) -> bool {
        F::relative_eq(&self.mu, &other.mu, epsilon, max_relative)
            && F::relative_eq(&self.sigma, &other.sigma, epsilon, max_relative)
            && F::relative_eq(&self.a, &other.a, epsilon, max_relative)
    }
}

impl<F: Float> UlpsEq for Gaussian<F>
where
    F::Epsilon: Copy,
{
    fn default_max_ulps() -> u32 {
        F::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: F::Epsilon, max_ulps: u32) -> bool {
        F::ulps_eq(&self.mu, &other.mu, epsilon, max_ulps)
            && F::ulps_eq(&self.sigma, &other.sigma, epsilon, max_ulps)
            && F::ulps_eq(&self.a, &other.a, epsilon, max_ulps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::TryInto;

    #[test]
    fn constructer() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let gaussian = Gaussian::new(mu, sigma, a);
        assert_eq!(gaussian.mu, mu);
        assert_eq!(gaussian.sigma, sigma);
        assert_eq!(gaussian.a, a);
    }

    #[test]
    fn getter() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let gaussian = Gaussian { mu, sigma, a };
        assert_eq!(gaussian.mu(), &mu);
        assert_eq!(gaussian.sigma(), &sigma);
        assert_eq!(gaussian.a(), &a);
    }

    #[test]
    fn setter() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let mut gaussian = Gaussian { mu, sigma, a };
        let (mu, sigma, a): (f64, f64, f64) = (4., 5., 6.);
        *gaussian.mu_mut() = mu;
        *gaussian.sigma_mut() = sigma;
        *gaussian.a_mut() = a;
        assert_eq!(gaussian.mu, mu);
        assert_eq!(gaussian.sigma, sigma);
        assert_eq!(gaussian.a, a);
    }

    #[test]
    fn as_tuple() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let gaussian = Gaussian { mu, sigma, a };
        assert_eq!(gaussian.as_tuple(), (&mu, &sigma, &a));
    }

    #[test]
    fn as_mut_tuple() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let mut gaussian = Gaussian { mu, sigma, a };
        let (mu, sigma, a): (f64, f64, f64) = (4., 5., 6.);
        let (mu_mut, sigma_mut, a_mut) = gaussian.as_mut_tuple();
        *mu_mut = mu;
        *sigma_mut = sigma;
        *a_mut = a;
        assert_eq!(gaussian.mu, mu);
        assert_eq!(gaussian.sigma, sigma);
        assert_eq!(gaussian.a, a);
    }

    #[test]
    fn to_tuple() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let gaussian = Gaussian { mu, sigma, a };
        assert_eq!(gaussian.to_tuple(), (mu, sigma, a));
    }

    #[test]
    fn into_tuple() {
        let (mu, sigma, a): (f64, f64, f64) = (1., 2., 3.);
        let gaussian = Gaussian { mu, sigma, a };
        assert_eq!(gaussian.into_tuple(), (mu, sigma, a));
    }

    #[test]
    fn from_into_tuple() {
        let input: (f64, f64, f64) = (1., 2., 3.);
        // From tuple
        let gaussian = Gaussian::from(input);
        assert_eq!(input, (gaussian.mu, gaussian.sigma, gaussian.a));

        // Into tuple
        let output: (_, _, _) = gaussian.into();
        assert_eq!(input, output);

        // From Gaussian
        // Skip

        // Into Gaussian
        let gaussian: Gaussian<_> = input.into();
        assert_eq!(input, (gaussian.mu, gaussian.sigma, gaussian.a));
    }

    #[test]
    fn from_into_ndarray1() {
        let input = array![1., 2., 3.];
        // TryFrom Array1 success
        let gaussian = Gaussian::try_from(input.clone()).unwrap();
        assert_eq!(input, array![gaussian.mu, gaussian.sigma, gaussian.a]);

        // TryFrom Array1 failure
        let input_failure = array![1.,];
        Gaussian::try_from(input_failure).unwrap_err();

        // Into Array1
        let output: Array1<_> = gaussian.into();
        assert_eq!(input, output);

        // From Gaussian
        let output = Array1::from(gaussian);
        assert_eq!(input, output);

        // TryInto Gaussian success
        let gaussian: Gaussian<_> = input.clone().try_into().unwrap();
        assert_eq!(input, array![gaussian.mu, gaussian.sigma, gaussian.a]);

        // TryInto Gaussian failure
        let input_failure = array![1.,];
        let result: Result<Gaussian<_>, _> = input_failure.try_into();
        result.unwrap_err();
    }
}
