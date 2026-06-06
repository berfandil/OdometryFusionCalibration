# sim — the simulation rig (relaxed edge). Never ships to the ECU: std containers /
# std::mt19937 / exceptions are fine. Depends only on the public ofc/core headers
# (via ofc_core) + Eigen. Built when OFC_BUILD_SIM is ON, and pulled into the test
# build by the root CMakeLists when OFC_BUILD_TESTS is ON (so the green gate covers it).
add_library(ofc_sim STATIC
    src/trajectory.cpp
    src/synthetic_source.cpp
    src/rig.cpp
)
target_include_directories(ofc_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ofc_sim PUBLIC ofc_core Eigen3::Eigen)

# Relaxed edge: C++14 (matches the core) but no strict-core warning flags.
target_compile_features(ofc_sim PUBLIC cxx_std_14)
set_target_properties(ofc_sim PROPERTIES
    CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
