add_executable(
    fdgmres_test
    fdgmres_test.cpp
    ${SRC_DIR}/common/memory_manager.cpp
    ${SRC_DIR}/common/linear_algebra.cpp
)
target_include_directories(
    fdgmres_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
)
target_link_libraries(
    fdgmres_test
    PRIVATE
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    fdgmres_test
    googletest
)
add_test(
    NAME fdgmres_test
    COMMAND $<TARGET_FILE:fdgmres_test>
)