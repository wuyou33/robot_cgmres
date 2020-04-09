add_executable(
    time_varying_smooth_horizon_test
    time_varying_smooth_horizon_test.cpp
    ${SRC_DIR}/ocp/time_varying_smooth_horizon.cpp
)
target_include_directories(
    time_varying_smooth_horizon_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
)
target_link_libraries(
    time_varying_smooth_horizon_test
    PRIVATE
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    time_varying_smooth_horizon_test
    googletest
)
add_test(
    NAME time_varying_smooth_horizon_test
    COMMAND $<TARGET_FILE:time_varying_smooth_horizon_test>
)
