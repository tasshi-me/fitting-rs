set -eux
export RUSTFLAGS="-Zinstrument-coverage"

export TARGET_BINARY_DIR="target/debug"
export COVERAGE_HTML_DIR=${TARGET_BINARY_DIR}"/coverage"
export LLVM_PROFILE_DIR=${COVERAGE_HTML_DIR}"/profraw"
export LLVM_PROFILE_FILE=${LLVM_PROFILE_DIR}"/tasshi-me-%p-%m.profraw"
export LLVM_PROFDATA_FILE=${LLVM_PROFILE_DIR}"/merged.profdata"

export TARGET_BINARY_FILE=${TARGET_BINARY_DIR}"/fitting"

rm -rf ${LLVM_PROFILE_DIR}
rm -rf ${COVERAGE_HTML_DIR}
mkdir -p ${LLVM_PROFILE_DIR}
mkdir -p ${COVERAGE_HTML_DIR}

cargo +nightly clean --verbose
cargo +nightly build --verbose
cargo +nightly run --verbose
cargo +nightly test --verbose

# html
rustup run nightly grcov . --binary-path ${TARGET_BINARY_DIR} -s . -t html --branch --ignore-not-existing --ignore "/*" -o ${COVERAGE_HTML_DIR}

# stdout
cargo +nightly profdata -- merge -sparse ${LLVM_PROFILE_DIR}/*  -o ${LLVM_PROFDATA_FILE}
cargo +nightly cov -- show ${TARGET_BINARY_FILE} \
-use-color \
-Xdemangler=rustfilt \
-instr-profile=${LLVM_PROFDATA_FILE} \
-show-line-counts-or-regions \
-show-instantiations \
--ignore-filename-regex="(.cargo|rustc)"
