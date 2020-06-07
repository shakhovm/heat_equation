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
#include "include/array2d/array2d.h"
#include <sstream>
#include "include/heat_equation/heat_equation.h"
#include "include/conf_handler/confhandler.h"
#include <chrono>


inline bool von_neumann_criterion(const ConfHandler::ConfParams& params)
{
    return params.delta_t > 0.25 * pow(std::max(params.delta_x, params.delta_y), 2)
                      * params.density * params.temp_capacity / params.temp_conduct;
}

eqution_params_t params_init(const ConfHandler::ConfParams& params) {

    double delta_x_sq = params.delta_x * params.delta_x,
           delta_y_sq = params.delta_y * params.delta_y,
           phys_params = params.temp_conduct / (params.density * params.temp_capacity);
    return eqution_params_t{1, params.rows - 1, 1, params.cols - 1, params.delta_t,
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
    auto start = std::chrono::high_resolution_clock::now();
    if (world.size() == 1) {
        one_process(params, conf_params.matrix_file, conf_params.out_file);
    }
    else if (world.size() == 2) {
        two_processes(params, conf_params.matrix_file, conf_params.out_file, world);
    }
    else
    {
        heat_equation(params, conf_params.matrix_file, world, conf_params.out_file, conf_params.rows, conf_params.cols);
    }
    if (world.rank() == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        size_t x = 0;

        std::cout << "The program has taken " << x << " millisecond";
    }
}
