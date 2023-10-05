cmake_minimum_required(VERSION 3.16)
project(remote_testing)

set(_components C CXX)
find_package(MPI COMPONENTS ${_components})
find_package(ADIOS2  REQUIRED COMPONENTS ${_components})


add_executable(test_vars test_vars.cpp)
target_link_libraries(test_vars adios2::adios2 MPI::MPI_C)
install(TARGETS test_vars
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
