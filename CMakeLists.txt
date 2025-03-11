# Top-level CMakeLists.txt

# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project name
project(highspeedFFT_project LANGUAGES C CXX)

# Set C and C++ standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1  # Updated to a valid tag
)
FetchContent_MakeAvailable(googletest)

# Fetch FFF (Fake Function Framework)
FetchContent_Declare(
  fff
  GIT_REPOSITORY https://github.com/meekrosoft/fff.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fff)

# Include subdirectories
add_subdirectory(src)
add_subdirectory(test)

# Enable testing
enable_testing()