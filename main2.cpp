#include <boost/mpi.hpp>
#include <iostream>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <exception>
#include <Magick++.h>
#include <cstdio>
#include <string>
#include <fstream>
#include "array2d.h"
#include <sstream>
#include "concur_queue.h"
#include <thread>
#include "heat_equation.h"
#include "include/conf_handler/confhandler.h"


inline bool von_neumann_criterion(const ConfHandler::ConfParams& params)
{
    return params.delta_t > 0.25 * pow(std::max(params.delta_x, params.delta_y), 2)
                      * params.density * params.temp_capacity / params.temp_conduct;
}

Array2D file_handler(const std::string& file_name) {
    size_t rows, cols;

    std::ifstream input_stream(file_name);
    if (!input_stream.is_open()) {
        throw std::runtime_error("No File with matrix found");
    }
    size_t buffer;
    input_stream >> rows >> cols;
    Array2D plate_matrix(cols, rows);
    for (size_t row = 0; row < rows; ++row)
    {
        for (size_t col = 0; col < cols; ++col)
        {
            input_stream >> buffer;
            plate_matrix(row, col) = buffer;
        }
    }
    return plate_matrix;
}

eqution_params_t params_init(const ConfHandler::ConfParams& params) {

    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 400, density = 8'900, temp_capacity = 1;
    double delta_x_sq = params.delta_x * params.delta_x,
           delta_y_sq = params.delta_y * params.delta_y,
           phys_params = params.temp_conduct / (params.density * params.temp_capacity);
    return eqution_params_t{1, 1, 1, 1, params.delta_t,
                phys_params, delta_x_sq, delta_y_sq,
                static_cast<size_t>(static_cast<double>(params.time_to_save)/params.delta_t),
                static_cast<size_t>(static_cast<double>(params.max_time)/params.delta_t)};
}


int main(int argc, char* argv[])
{
    std::string conf_file = argc > 1 ? argv[1] : "conf.dat";
    ConfHandler conf_handler;
    try {
       conf_handler = ConfHandler(conf_file);
    } catch (std::runtime_error& e) {
        std::cout << "Failed to handle conf file" << std::endl;
        return 1;
    }

    auto conf_params = conf_handler.getConfParams();
//    if (!von_neumann_criterion(conf_params)) {
//        std::cout << "Bad von-neumann criterion" << std::endl;
//        return -1;
//    }

    boost::mpi::environment env{argc, argv};
    boost::mpi::communicator world;
    eqution_params_t params = params_init(conf_params);
    Array2D plate_matrix;
    try {
        plate_matrix = file_handler(conf_params.matrix_file);
    } catch (std::runtime_error& e) {
        std::cout << "No File with matrix found" << std::endl;
        return -2;
    }

    if (world.size() == 1)
        return 0;
    else
        heat_equation(params, plate_matrix, world, conf_params.out_file);
}

//template <typename T>
//Magick::ColorRGB heatmap_color(const T& value, const T& min_val, const T& max_val)
//{
//    double red, green, blue;
//    double ratio = 2 * static_cast<double>(value - min_val) / (max_val - min_val);
//    blue = std::max(0., 1 - ratio);
//    red = std::max(0., ratio - 1);
//    green = 1 - blue - red;
//    return Magick::ColorRGB(red, green, blue);
//}
