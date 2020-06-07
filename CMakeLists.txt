cmake_minimum_required(VERSION 3.5)

project(mpi_lab LANGUAGES CXX)

find_package(Boost COMPONENTS mpi serialization REQUIRED)
find_package(ImageMagick COMPONENTS Magick++ MagickWand MagickCore REQUIRED)
find_package(Threads REQUIRED)
include_directories(${ImageMagick_INCLUDE_DIRS})
link_directories(${ImageMagick_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC src/array2d/array2d.cpp src/heat_equation/heat_equation.cpp src/conf_handler/confhandler.cpp src/heat_equation/heat_equation_functions.cpp)
set(INC include/array2d/array2d.h include/heat_equation/heat_equation.h include/conf_handler/confhandler.h include/heat_equation/heat_equation_functions.h)

add_executable(mpi_lab main.cpp ${SRC} ${INC})

target_link_libraries(mpi_lab Boost::mpi Boost::serialization ${ImageMagick_LIBRARIES} Threads::Threads)
