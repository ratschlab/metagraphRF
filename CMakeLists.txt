cmake_minimum_required(VERSION 3.16)

project(metagraphRF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add metagraph as a subdirectory (from git submodule)
add_subdirectory(metagraph/metagraph)

# Create example executables
add_executable(basic_query src/pybindings.cpp)