cmake_minimum_required(VERSION 3.17)
project(lzw_pattern_matching VERSION "1.0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(main
  algo.cpp
)

find_package(spdlog CONFIG REQUIRED)
target_link_libraries(main PRIVATE spdlog::spdlog)

find_package(GTest CONFIG REQUIRED)

add_executable(
  suffix_tree_test
  suffix_tree_test.cpp
)

target_link_libraries(
  suffix_tree_test
  GTest::gtest_main
)
