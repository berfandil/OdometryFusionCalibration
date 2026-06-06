cmake_minimum_required(VERSION 3.16)
project(OdometryFusionCalibration LANGUAGES CXX VERSION 0.0.1)

include(FetchContent)

# ---- Options ---------------------------------------------------------------
option(OFC_BUILD_ADAPTERS "Build middleware adapters (relaxed edge)" OFF)
option(OFC_BUILD_SIM      "Build the simulation rig (relaxed edge)"  OFF)
option(OFC_BUILD_TESTS    "Build unit/observability/consistency tests" ON)
option(OFC_STRICT_CORE    "Compile the core with strict warning flags" ON)

# ---- Eigen (header-only) ---------------------------------------------------
# SO(3)/SE(3) ops are in-house (no Sophus). Use system Eigen if present, else fetch.
find_package(Eigen3 3.3 QUIET)
if (NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found — fetching 3.4.0")
    FetchContent_Declare(eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING   OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
endif()

# ---- Core library (STRICT) -------------------------------------------------
# C++14 by constraint (AUTOSAR). Strict warning flags below. True no-exceptions/
# no-RTTI enforcement is the cross-compile target's job; this dev build keeps
# exceptions on for Eigen/MSVC ergonomics (strict core, relaxed dev toolchain).
add_library(ofc_core STATIC
    src/core/lie.cpp
    src/core/buffer.cpp
    src/core/median.cpp
    src/core/eskf.cpp
    src/core/estimator.cpp
    src/core/histogram.cpp
    src/core/timesync.cpp
)
target_include_directories(ofc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ofc_core PUBLIC Eigen3::Eigen)

target_compile_features(ofc_core PUBLIC cxx_std_14)
set_target_properties(ofc_core PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

if (OFC_STRICT_CORE)
    if (MSVC)
        target_compile_options(ofc_core PRIVATE /W4 /permissive-)
    else()
        target_compile_options(ofc_core PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
endif()

# ---- Relaxed edges ---------------------------------------------------------
if (OFC_BUILD_ADAPTERS)
    add_subdirectory(adapters)
endif()
# The sim rig is the validation backbone (D24): the tests link it, so build it
# whenever tests are on, even if OFC_BUILD_SIM was not requested explicitly.
if (OFC_BUILD_SIM OR OFC_BUILD_TESTS)
    add_subdirectory(sim)
endif()
if (OFC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
