# src/CMakeLists.txt

# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project name and specify the language
project(highspeedFFT C)

# Set the C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Define the source directory
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Add all source files from the src directory
set(SOURCES
    ${SRC_DIR}/highspeedFFT.c
    ${SRC_DIR}/convolve.c
    ${SRC_DIR}/real.c
)

# Add all header files (optional, but good practice for IDEs)
set(HEADERS
    ${SRC_DIR}/highspeedFFT.h
    ${SRC_DIR}/real.h
)

# Create a static library instead of an executable
add_library(highspeedFFT STATIC ${SOURCES} ${HEADERS})

# Add the src directory to the include path so headers can be found
target_include_directories(highspeedFFT PUBLIC ${SRC_DIR})

# Link math library and ensure C library linkage
target_link_libraries(highspeedFFT PRIVATE m)

# (Optional) Add compiler flags if needed
# Example: target_compile_options(highspeedFFT PRIVATE -Wall -Wextra)