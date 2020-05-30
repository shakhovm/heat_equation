cmake_minimum_required(VERSION 3.5)

project(mpi_lab LANGUAGES CXX)

find_package(Boost COMPONENTS mpi serialization REQUIRED)
find_package(JPEG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC array2d.cpp)
set(INC array2d.h)

add_executable(mpi_lab main.cpp ${SRC} ${INC})
add_executable(mpi_lab2 main2.cpp ${SRC} ${INC})
target_link_libraries(mpi_lab Boost::mpi Boost::serialization JPEG::JPEG)
target_link_libraries(mpi_lab2 Boost::mpi Boost::serialization JPEG::JPEG)
