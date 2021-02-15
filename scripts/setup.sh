set -eux
rustup toolchain install nightly
rustup run nightly rustup component add llvm-tools-preview
cargo +nightly install cargo-binutils
cargo +nightly install grcov
cargo +nightly install rustfilt
