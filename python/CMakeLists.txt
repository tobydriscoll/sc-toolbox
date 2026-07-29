cmake_minimum_required(VERSION 3.16)
project(sctoolbox_python CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Reuse the C++ library target from ../cpp. Turn its Catch2 test target off so
# the Python build needs only Eigen (fetched automatically if not installed).
set(SCTOOLBOX_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cpp
                 ${CMAKE_CURRENT_BINARY_DIR}/cpp)

# Locate Python and nanobind. When built through scikit-build-core these are
# provided automatically; the find_package calls also allow a plain
# `cmake -S python -B build` for development.
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)

if(NOT nanobind_ROOT)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
endif()
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_sctoolbox STABLE_ABI NB_STATIC
  src/sctoolbox_ext.cpp)
target_link_libraries(_sctoolbox PRIVATE sctoolbox)

# Install into the `sctoolbox` package (scikit-build-core places the wheel
# contents relative to this install prefix).
install(TARGETS _sctoolbox LIBRARY DESTINATION sctoolbox)
