# test/CMakeLists.txt

# Define the test executable
add_executable(test_mixedRadixFFT
    test_mixedRadixFFT.cpp
)

# Include directories
target_include_directories(test_mixedRadixFFT PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${fff_SOURCE_DIR}
    ${SRC_DIR}/src  # Include src directory for headers
)

# Link libraries
target_link_libraries(test_mixedRadixFFT
    highspeedFFT
    gtest_main
    gmock_main
)

# Register tests
include(GoogleTest)
gtest_discover_tests(test_mixedRadixFFT)