cmake_minimum_required(VERSION 3.22)
project(MetagraphRF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(extern/parlaylib/include)
add_subdirectory(extern/metagraph/metagraph)
add_subdirectory(extern/pybind11)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")

set(PYBIND11_FINDPYTHON ON)
pybind11_add_module(_core src/pybindings.cpp)

target_link_libraries(_core PUBLIC metagraph-cli metagraph-core pybind11::headers)
install(TARGETS _core DESTINATION rfmetagraph)
