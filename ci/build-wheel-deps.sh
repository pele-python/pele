#!/usr/bin/env bash
set -euo pipefail

prefix=/tmp/pele-wheel-deps
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cd "$work"

fetch() {
    local url=$1 filename=$2 sha=$3
    curl --fail --location --retry 3 "$url" -o "$filename"
    test "$(openssl dgst -sha256 "$filename" | awk '{print $NF}')" = "$sha"
    tar -xf "$filename"
}

fetch https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.34/OpenBLAS-0.3.34.tar.gz OpenBLAS.tar.gz cd7e129868320cc2d033afa920e31202dfe0b8066a5b66661900ccc0f197dfed
fetch https://github.com/LLNL/sundials/releases/download/v7.9.0/sundials-7.9.0.tar.gz sundials.tar.gz 13f898a27b48fe3449483f9e438a800ed545abf93bc2e2ceec2d1e00ae8db5ef
fetch https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz eigen.tar.gz 8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72

if [[ $(uname -s) == Darwin ]]; then
    fetch https://github.com/llvm/llvm-project/releases/download/llvmorg-19.1.7/openmp-19.1.7.src.tar.xz openmp.tar.xz bd7e6901ab086fd268750363017935fd4a717c153dad3c2aab86cb0140d9e3fe
    fetch https://github.com/llvm/llvm-project/releases/download/llvmorg-19.1.7/cmake-19.1.7.src.tar.xz cmake.tar.xz 11c5a28f90053b0c43d0dec3d0ad579347fc277199c005206b963c19aae514e3
    mv cmake-19.1.7.src cmake
    # Bootstrap libomp before enabling OpenMP flags for the other libraries.
    CFLAGS= CXXFLAGS= LDFLAGS= cmake -S openmp-19.1.7.src -B openmp-build \
        -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" \
        -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 -DLIBOMP_ENABLE_SHARED=ON \
        -DOPENMP_ENABLE_LIBOMPTARGET=OFF -DOPENMP_ENABLE_OMPT_TOOLS=OFF \
        -DLIBOMP_OMPD_SUPPORT=OFF
    cmake --build openmp-build --parallel 2
    cmake --install openmp-build
fi

# Dynamic dispatch selects supported kernels; the baseline must remain generic.
target=GENERIC
if [[ $(uname -m) == arm64 || $(uname -m) == aarch64 ]]; then target=ARMV8; fi
cmake -S OpenBLAS-0.3.34 -B openblas-build -G 'Unix Makefiles' \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_INSTALL_LIBDIR=lib \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DBUILD_TESTING=OFF \
    -DDYNAMIC_ARCH=ON -DTARGET="$target" -DINTERFACE64=OFF -DUSE_OPENMP=OFF \
    -DBUILD_WITHOUT_LAPACK=OFF -DBUILD_WITHOUT_LAPACKE=OFF
cmake --build openblas-build --parallel 2
cmake --install openblas-build
# Keep LAPACK discovery on the pinned provider, including on macOS.
if [[ $(uname -s) == Darwin ]]; then
    ln -s libopenblas.dylib "$prefix/lib/liblapack.dylib"
fi

cmake -S sundials-7.9.0 -B sundials-build -G 'Unix Makefiles' \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_INSTALL_LIBDIR=lib \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DSUNDIALS_PRECISION=DOUBLE \
    -DSUNDIALS_ENABLE_CVODE=ON -DSUNDIALS_ENABLE_CVODES=OFF \
    -DSUNDIALS_ENABLE_ARKODE=OFF -DSUNDIALS_ENABLE_IDA=OFF \
    -DSUNDIALS_ENABLE_IDAS=OFF -DSUNDIALS_ENABLE_KINSOL=OFF \
    -DSUNDIALS_ENABLE_FORTRAN=OFF -DSUNDIALS_ENABLE_C_EXAMPLES=OFF \
    -DSUNDIALS_ENABLE_CXX_EXAMPLES=OFF -DSUNDIALS_ENABLE_EXAMPLES_INSTALL=OFF
cmake --build sundials-build --parallel 2
cmake --install sundials-build

mkdir -p "$prefix/include/eigen3"
cp -R eigen-3.4.0/Eigen eigen-3.4.0/unsupported "$prefix/include/eigen3/"
