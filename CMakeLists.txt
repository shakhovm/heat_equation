cmake_minimum_required(VERSION 3.5)

project(mpi_lab LANGUAGES CXX)

find_package(Boost COMPONENTS mpi serialization REQUIRED)
find_package(ImageMagick COMPONENTS Magick++ MagickWand convert MagickCore REQUIRED)
find_package(Threads REQUIRED)
include_directories(${ImageMagick_INCLUDE_DIRS})
link_directories(${ImageMagick_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC array2d.cpp concur_queue.cpp heat_equation.cpp )
set(INC array2d.h concur_queue.h heat_equation.h)

add_executable(mpi_lab main2.cpp ${SRC} ${INC})

target_link_libraries(mpi_lab Boost::mpi Boost::serialization ${ImageMagick_LIBRARIES} Threads::Threads)
