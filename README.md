# fitting-rs

[![crates.io](https://img.shields.io/crates/v/fitting.svg)](https://crates.io/crates/fitting)
[![docs.rs](https://docs.rs/fitting/badge.svg)](https://docs.rs/fitting)
[![Build and Test](https://github.com/mshrtsr/fitting-rs/workflows/Build%20and%20Test/badge.svg)](https://github.com/mshrtsr/fitting-rs/actions?query=workflow%3A%22Build+and+Test%22)
[![Format and Lint](https://github.com/mshrtsr/fitting-rs/workflows/Format%20and%20Lint/badge.svg)](https://github.com/mshrtsr/fitting-rs/actions?query=workflow%3A%22Format+and+Lint%22)
[![codecov](https://codecov.io/gh/mshrtsr/fitting-rs/branch/master/graph/badge.svg)](https://codecov.io/gh/mshrtsr/fitting-rs)

Curve fitting library for Rust

## Updates

### 0.3.0

- Migrate from the `failure` crate to `thiserror`.
- Refactor some tests.

### 0.2.1

- Error handing changed. Some functions returns Result instead of Option.
- linalg.solve() is improved. Now it can solve NxM array with pivoting.

### 0.2.0

- Using ndarray instead of nested Vec
- Improvement of unit test
- Add status badges

### 0.1.0

- Implements linalg solve and gaussian fit.

## License

This project is licensed under the MIT license.
