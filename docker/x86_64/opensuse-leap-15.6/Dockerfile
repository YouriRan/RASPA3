FROM opensuse/leap:15.6

RUN \
  zypper update -y && \
  zypper install -y --no-recommends wget git tar ca-certificates libstdc++-devel gcc gcc-c++ binutils binutils-devel openssl-devel && \
  zypper install -y --no-recommends python3 python3-pybind11-common-devel python3-pybind11 rpm-build python3-devel && \
  zypper install -y --no-recommends lapack-devel blas-devel && \
  zypper install -y --no-recommends cmake ninja libhdf5_cpp103 libhdf5-103 && \
  wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.12.1.tar.gz && \
  tar -zxvf v1.12.1.tar.gz && \
  cd ninja-1.12.1 && \
  cmake -Bbuild-cmake && \
  cmake --build build-cmake && \
  cmake --install build-cmake && \
  cd .. && \
  rm -rf ninja-1.12.1 && \
  rm  v1.12.1.tar.gz && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-18.1.8.tar.gz && \
  tar -zxvf llvmorg-18.1.8.tar.gz && \
  cd llvm-project-llvmorg-18.1.8 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build && \
  ninja -C build install && \
  cd .. && \
  rm -rf llvm* && \
  ldconfig && \
  wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.4.3/hdf5.tar.gz && \
  tar -zxvf hdf5.tar.gz && \
  cd hdf5-1.14.4-3 &&\
  cmake -B build -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DCMAKE_INSTALL_PREFIX=/usr -DHDF5_BUILD_CPP_LIB=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++23 -stdlib=libc++"  -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF && \
  cd build && \
  make -j 16 && \
  make install && \
  cd ../.. && \
  rm -rf hdf5-1.14.4-3 && \
  ldconfig && \
  zypper clean
