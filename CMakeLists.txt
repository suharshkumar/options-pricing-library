cmake_minimum_required(VERSION 3.14)
project(options_pricer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra)

find_package(Threads REQUIRED)

# --- Core library -----------------------------------------------------------
add_library(pricer src/BlackScholes.cpp src/Binomial.cpp src/MonteCarlo.cpp)
target_include_directories(pricer PUBLIC include)
target_link_libraries(pricer PUBLIC Threads::Threads)

# --- Demo executable --------------------------------------------------------
add_executable(demo src/main.cpp)
target_link_libraries(demo PRIVATE pricer)

# --- Tests (GoogleTest, fetched automatically) ------------------------------
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(pricer_tests tests/test_pricer.cpp)
target_link_libraries(pricer_tests PRIVATE pricer GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(pricer_tests)
