# adapters — middleware integrations (RELAXED EDGE, Slice 13). Never ships to the ECU as
# the strict core does: std containers / std::fstream / std::thread / exceptions / heap are
# all fine here. Depends only on the PUBLIC ofc/core headers (via ofc_core) + Eigen. Built
# when OFC_BUILD_ADAPTERS is ON (the gate sets it; see scripts/dev.ps1).
#
# Three production adapters ship here (file-persistence double-buffer, a threading wrapper,
# a dependency-free config loader). The ROS node is a DOCUMENTED, compile-guarded header
# SKETCH only (adapters/include/ofc_adapters/ros_node.hpp) — NOT compiled into this library
# (no ROS on this box); the real ROS node + recorded-bag round-trip is deferred to Slice 13b.
add_library(ofc_adapters STATIC
    src/file_persistence.cpp
    src/threaded_estimator.cpp
    src/config_loader.cpp
    src/gps_correction.cpp
)
target_include_directories(ofc_adapters PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ofc_adapters PUBLIC ofc_core Eigen3::Eigen)

# std::thread needs the platform thread library where one exists (no-op on MSVC).
find_package(Threads QUIET)
if (Threads_FOUND)
    target_link_libraries(ofc_adapters PUBLIC Threads::Threads)
endif()

# Relaxed edge: C++14 (matches the core) but NO strict-core warning flags.
target_compile_features(ofc_adapters PUBLIC cxx_std_14)
set_target_properties(ofc_adapters PROPERTIES
    CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# ---- Adapters doctest target -----------------------------------------------------------
# Self-contained (links only ofc_core + ofc_adapters + doctest — NOT ofc_sim), so it builds
# regardless of whether the sim/tests subdirs were configured. doctest is fetched the same
# way tests/ does; FetchContent_MakeAvailable is idempotent (keyed by name) so declaring it
# here too is harmless when tests/ also declares it.
if (OFC_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(doctest
        URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.11.tar.gz)
    FetchContent_MakeAvailable(doctest)

    add_executable(ofc_adapters_tests
        tests/test_main.cpp
        tests/test_file_persistence.cpp
        tests/test_threaded_estimator.cpp
        tests/test_config_loader.cpp
        tests/test_gps_correction.cpp
    )
    target_link_libraries(ofc_adapters_tests PRIVATE ofc_core ofc_adapters doctest::doctest)
    set_target_properties(ofc_adapters_tests PROPERTIES
        CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON)

    # The root adds this subdir BEFORE its own enable_testing() (the adapters block precedes the
    # tests block in the root CMakeLists), so register testing here too — enable_testing() is
    # idempotent and add_test() needs it active in this scope for CTest to pick up the test.
    enable_testing()
    add_test(NAME adapters COMMAND ofc_adapters_tests)
endif()
