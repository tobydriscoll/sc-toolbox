cmake_minimum_required(VERSION 3.16)
project(sctoolbox_cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Tests need Catch2; consumers who add_subdirectory() this project (e.g. the
# Python bindings in ../python) do not, so default the suite OFF in that case.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(SCTOOLBOX_BUILD_TESTS_DEFAULT ON)
else()
  set(SCTOOLBOX_BUILD_TESTS_DEFAULT OFF)
endif()
option(SCTOOLBOX_BUILD_TESTS "Build the golden-value test suite" ${SCTOOLBOX_BUILD_TESTS_DEFAULT})

# Eigen: prefer a system/Homebrew install (the documented macOS workflow), but
# fall back to fetching a pinned release so the library builds on machines
# without Eigen installed (e.g. CI, or the Python wheel build).
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found; fetching a pinned copy via FetchContent")
  include(FetchContent)
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72)
  FetchContent_MakeAvailable(eigen)
  add_library(eigen_headers INTERFACE)
  target_include_directories(eigen_headers INTERFACE ${eigen_SOURCE_DIR})
  add_library(Eigen3::Eigen ALIAS eigen_headers)
endif()

add_library(sctoolbox
  src/gaussj.cpp
  src/scqdata.cpp
  src/scangle.cpp
  src/neqrdcmp.cpp
  src/nechdcmp.cpp
  src/nefdjac.cpp
  src/nesolve.cpp
  src/scaddvtx.cpp
  src/scfix.cpp
  src/isinpoly.cpp
  src/polygon.cpp
  src/dquad.cpp
  src/hpquad.cpp
  src/dequad.cpp
  src/stderiv.cpp
  src/stquad.cpp
  src/crquad.cpp
  src/dderiv.cpp
  src/hpderiv.cpp
  src/dederiv.cpp
  src/stquadh.cpp
  src/dmap.cpp
  src/hpmap.cpp
  src/demap.cpp
  src/stmap.cpp
  src/scinvopt.cpp
  src/ode45.cpp
  src/findz0.cpp
  src/dinvmap.cpp
  src/hpinvmap.cpp
  src/deinvmap.cpp
  src/stinvmap.cpp
  src/ellipkkp.cpp
  src/ellipjc.cpp
  src/r2strip.cpp
  src/rcorners.cpp
  src/rectmap_internal.cpp
  src/rderiv.cpp
  src/rmap.cpp
  src/rinvmap.cpp
  src/moebius3.cpp
  src/crembed.cpp
  src/crspread.cpp
  src/crmap0.cpp
  src/crderiv.cpp
  src/crmap.cpp
  src/crgather.cpp
  src/crimap0.cpp
  src/crinvmap.cpp
  src/dabsquad.cpp
  src/dparam.cpp
  src/hpparam.cpp
  src/deparam.cpp
  src/stparam.cpp
  src/rparam.cpp
  src/crpsdist.cpp
  src/crossrat.cpp
  src/crtriang.cpp
  src/crcdt.cpp
  src/crqgraph.cpp
  src/crsplit.cpp
  src/craffine.cpp
  src/crfixwc.cpp
  src/crparam.cpp
  src/diskmap.cpp
  src/hplmap.cpp
  src/extermap.cpp
  src/stripmap.cpp
  src/rectmap.cpp
  src/crdiskmap.cpp
  src/annulusmap_internal.cpp
  src/annulusmap.cpp
  src/moebius.cpp
  src/composite.cpp
)
target_include_directories(sctoolbox PUBLIC include)
target_link_libraries(sctoolbox PUBLIC Eigen3::Eigen)
# Position-independent code so the static library can be linked into the
# Python extension module (a shared object).
set_target_properties(sctoolbox PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(SCTOOLBOX_BUILD_TESTS)
  find_package(Catch2 3 REQUIRED)

  add_executable(sctoolbox_tests
    tests/golden_reader.hpp
    tests/test_gaussj.cpp
    tests/test_scqdata.cpp
    tests/test_scangle.cpp
    tests/test_neqrdcmp.cpp
    tests/test_nechdcmp.cpp
    tests/test_nefdjac.cpp
    tests/test_nesolve.cpp
    tests/test_scfix.cpp
    tests/test_isinpoly.cpp
    tests/test_polygon.cpp
    tests/test_dquad.cpp
    tests/test_hpquad.cpp
    tests/test_dequad.cpp
    tests/test_stderiv.cpp
    tests/test_stquad.cpp
    tests/test_crquad.cpp
    tests/test_dderiv.cpp
    tests/test_hpderiv.cpp
    tests/test_dederiv.cpp
    tests/test_dmap.cpp
    tests/test_hpmap.cpp
    tests/test_demap.cpp
    tests/test_stmap.cpp
    tests/test_dinvmap.cpp
    tests/test_hpinvmap.cpp
    tests/test_deinvmap.cpp
    tests/test_stinvmap.cpp
    tests/test_ellipkkp.cpp
    tests/test_ellipjc.cpp
    tests/test_rderiv.cpp
    tests/test_rmap.cpp
    tests/test_rinvmap.cpp
    tests/test_crderiv.cpp
    tests/test_crmap.cpp
    tests/test_crinvmap.cpp
    tests/test_dabsquad.cpp
    tests/test_dparam.cpp
    tests/test_hpparam.cpp
    tests/test_deparam.cpp
    tests/test_stparam.cpp
    tests/test_rparam.cpp
    tests/test_crparam.cpp
    tests/test_diskmap.cpp
    tests/test_hplmap.cpp
    tests/test_extermap.cpp
    tests/test_stripmap.cpp
    tests/test_rectmap.cpp
    tests/test_crdiskmap.cpp
    tests/test_annulusmap.cpp
    tests/test_moebius.cpp
    tests/test_composite.cpp
  )
  target_link_libraries(sctoolbox_tests PRIVATE sctoolbox Catch2::Catch2WithMain)
  target_compile_definitions(sctoolbox_tests PRIVATE
    GOLDENS_TEXT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/goldens_text")

  include(CTest)
  include(Catch)
  enable_testing()
  catch_discover_tests(sctoolbox_tests)
endif()
