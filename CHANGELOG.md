# Changelog

## [0.3.0](https://github.com/mshrtsr/fitting-rs/compare/0.2.1...0.3.0) (2020-05-16)

- Migrate from the `failure` crate to `thiserror`.
  - https://crates.io/crates/thiserror
- Refactor some tests.


## [0.2.1](https://github.com/mshrtsr/fitting-rs/compare/0.2.0...0.2.1) (2019-12-04)

- Error handing changed. Some functions returns Result instead of Option.
- linalg.solve() is improved. Now it can solve NxM array with pivoting.


## [0.2.0](https://github.com/mshrtsr/fitting-rs/compare/0.1.0...0.2.0) (2019-11-08)

- Using [ndarray](https://crates.io/crates/ndarray) instead of nested Vec
- Improvement of unit test
- Add status badges


## 0.1.0 (2019-11-08)

- Implements linalg solve and gaussian fit.
