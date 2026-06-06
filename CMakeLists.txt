cmake_minimum_required(VERSION 3.16)
project(OdometryFusionCalibration LANGUAGES CXX VERSION 0.0.1)

# ---- Options ---------------------------------------------------------------
option(OFC_BUILD_ADAPTERS "Build middleware adapters (relaxed edge)" OFF)
option(OFC_BUILD_SIM      "Build the simulation rig (relaxed edge)"  OFF)
option(OFC_BUILD_TESTS    "Build unit/observability/consistency tests" OFF)
option(OFC_STRICT_CORE    "Compile the core to the AUTOSAR-leaning strict profile" ON)

# ---- Core library (STRICT) -------------------------------------------------
# C++14 by constraint (AUTOSAR). The core must build without exceptions/RTTI and
# must not allocate after init(). Strict warning flags applied below.
add_library(ofc_core STATIC
    src/core/estimator.cpp
    # src/core/buffer.cpp        # Slice 1
    # src/core/median.cpp        # Slice 2
    # src/core/eskf.cpp          # Slice 2
    # src/core/histogram.cpp     # Slice 4
    # src/core/timesync.cpp      # Slice 5
    # src/core/calibration.cpp   # Slices 6-8
    # src/core/lifecycle.cpp     # Slice 3
)
target_include_directories(ofc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_features(ofc_core PUBLIC cxx_std_14)
set_target_properties(ofc_core PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Eigen (header-only) — linear algebra. SO(3)/SE(3) ops are in-house (no Sophus).
find_package(Eigen3 QUIET)
if (Eigen3_FOUND)
    target_link_libraries(ofc_core PUBLIC Eigen3::Eigen)
else()
    message(WARNING "Eigen3 not found — ofc_core will not build until it is available.")
endif()

if (OFC_STRICT_CORE)
    if (MSVC)
        target_compile_options(ofc_core PRIVATE /W4 /permissive-)
    else()
        target_compile_options(ofc_core PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow
            -fno-exceptions -fno-rtti)
    endif()
endif()

# ---- Relaxed edges ---------------------------------------------------------
if (OFC_BUILD_ADAPTERS)
    add_subdirectory(adapters)
endif()
if (OFC_BUILD_SIM)
    add_subdirectory(sim)
endif()
if (OFC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
