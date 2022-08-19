# Changelog

## [0.4.0](https://github.com/mshrtsr/fitting-rs/compare/0.3.0...0.4.0) (2022-08-19)


### ⚠ BREAKING CHANGES

* **gaussian:** Redesign API (#15)

### Features

* **gaussian:** Redesign API ([#15](https://github.com/mshrtsr/fitting-rs/issues/15)) ([14c9340](https://github.com/mshrtsr/fitting-rs/commit/14c9340b046c3e47086ae685705acb72faf25a50))


### Bug Fixes

* coverage using cargo-kcov ([d1ce724](https://github.com/mshrtsr/fitting-rs/commit/d1ce724c8482288ca4e98b0bf30b76531cec505a))


### Chores

* add CHANGELOG.md ([adefab3](https://github.com/mshrtsr/fitting-rs/commit/adefab34cd8171e54d37172ebeba8cccf93b13f7))
* bootstrap releases for path: . ([#32](https://github.com/mshrtsr/fitting-rs/issues/32)) ([7096e8c](https://github.com/mshrtsr/fitting-rs/commit/7096e8c4aa13e7c66980df713b34cc9e8a4e5b43))
* **deps:** update actions/cache action to v3 ([#28](https://github.com/mshrtsr/fitting-rs/issues/28)) ([aae4a4b](https://github.com/mshrtsr/fitting-rs/commit/aae4a4bf3c2e4bfa4c49a6acbeb4119e8f0c2b5b))
* **deps:** update actions/checkout action to v3 ([#29](https://github.com/mshrtsr/fitting-rs/issues/29)) ([896bf16](https://github.com/mshrtsr/fitting-rs/commit/896bf16d27812a98b1e55a252ba6ae34c25c1921))
* updated Cargo.toml (exclude section) ([fe52125](https://github.com/mshrtsr/fitting-rs/commit/fe52125da9fd3312a4053b9d2a47c864a238d56d))
* updated README.md ([7422682](https://github.com/mshrtsr/fitting-rs/commit/7422682a6c08f57e191c037fdacf0554ec52de4f))

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
