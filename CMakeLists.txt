cmake_minimum_required(VERSION 3.5)

project(mpi_lab LANGUAGES CXX)

find_package(Boost COMPONENTS mpi serialization REQUIRED)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mpi_lab main.cpp)

target_link_libraries(mpi_lab Boost::mpi Boost::serialization)
