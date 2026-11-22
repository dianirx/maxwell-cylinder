cmake_minimum_required(VERSION 3.14)
project(hdg-maxwell-cylinder CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(hdg_cylinder src/main.cpp)
target_include_directories(hdg_cylinder PRIVATE include)

# No external dependencies, C++ standard library only.
