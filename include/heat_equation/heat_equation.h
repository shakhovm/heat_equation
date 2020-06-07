#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H
#include "heat_equation_functions.h"


void heat_equation(eqution_params_t params, const std::string& matrix_file, boost::mpi::communicator& world,
                   const std::string& out_file, size_t rows, size_t cols);

void main_process(eqution_params_t params, Array2D& plate_matrix, const std::string& out_file,
                  boost::mpi::communicator& world, size_t matrix_height);

void one_process(eqution_params_t params, const std::string& matrix_file, const std::string& out_file);

void two_processes(eqution_params_t params, const std::string& matrix_file, const std::string& out_file,
                   boost::mpi::communicator& world);

template <typename T>
Magick::ColorRGB heatmap_color(const T& value, const T& min_val, const T& max_val);

#endif // HEAT_EQUATION_H
