# tests — relaxed edge. doctest via FetchContent.
include(FetchContent)
FetchContent_Declare(doctest
    URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.11.tar.gz)
FetchContent_MakeAvailable(doctest)

add_executable(ofc_tests
    test_main.cpp
    test_lie.cpp
    test_config.cpp
    test_buffer.cpp
    test_median.cpp
    test_median_split.cpp
    test_split_median.cpp
    test_eskf.cpp
    test_fusion.cpp
    test_lifecycle.cpp
    test_histogram.cpp
    test_timesync.cpp
    test_calib_phase1.cpp
    test_calib_phase2.cpp
    test_calib_rot3d.cpp
    test_calib_joint_lever_scale.cpp
    test_calib_translation_only.cpp
    test_calib_vote_guard.cpp
    test_calib_feedback.cpp
    test_smoother.cpp
    test_smoother_calib.cpp
    test_weights.cpp
    test_weights_split.cpp
    test_heading_monitor.cpp
    test_sim.cpp
    test_validation.cpp
    test_correction.cpp
    test_chi2_gate.cpp
    test_bias.cpp
    test_multi_bias.cpp
    test_multi_bias_sim.cpp
    test_persistence.cpp
    test_cov_calibration.cpp
    test_cov_calibration_split.cpp
    test_placeholder_sweep.cpp
)
target_link_libraries(ofc_tests PRIVATE ofc_core ofc_sim doctest::doctest)
set_target_properties(ofc_tests PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON)

add_test(NAME unit COMMAND ofc_tests)
