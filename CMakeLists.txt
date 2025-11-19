cmake_minimum_required(VERSION 3.16)

project(metagraphRF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add metagraph as a subdirectory (from git submodule)
add_subdirectory(metagraph/metagraph)

# Create example executables
add_executable(basic_query src/pybindings.cpp)

# Link against metagraph libraries - includes are inherited automatically
target_link_libraries(basic_query PRIVATE
        metagraph-cli
        metagraph-core
)

# Set alphabet (DNA or Protein)
if(NOT DEFINED CMAKE_DBG_ALPHABET)
    set(CMAKE_DBG_ALPHABET "DNA" CACHE STRING "Alphabet for MetaGraph (DNA or Protein)")
endif()