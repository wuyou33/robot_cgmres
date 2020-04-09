add_executable(
    memory_manager_test
    memory_manager_test.cpp
    ${SRC_DIR}/common/memory_manager.cpp
)
target_include_directories(
    memory_manager_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
)
target_link_libraries(
    memory_manager_test
    PRIVATE
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    memory_manager_test
    googletest
)
add_test(
    NAME memory_manager_test 
    COMMAND $<TARGET_FILE:memory_manager_test>
)

add_executable(
    linear_algebra_test
    linear_algebra_test.cpp
    ${SRC_DIR}/common/linear_algebra.cpp
    ${SRC_DIR}/common/memory_manager.cpp
)
target_include_directories(
    linear_algebra_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
)
target_link_libraries(
    linear_algebra_test
    PRIVATE
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    linear_algebra_test
    googletest
)
add_test(
    NAME linear_algebra_test
    COMMAND $<TARGET_FILE:linear_algebra_test>
)