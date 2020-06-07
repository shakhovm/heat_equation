#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H
#include "array2d.h"
#include <boost/mpi.hpp>
#include <Magick++.h>
#include "concur_queue.h"
struct eqution_params_t {
    size_t row_start,
           row_end,
           col_start,
           col_end;

    double delta_t,
           phys_params,
           delta_x_sq,
           delta_y_sq;

    size_t iteration_to_save,
           max_iterations;

};

void bound_calculate(size_t row, size_t col, const Array2D& plate_matrix, Array2D& plate_buffer, const eqution_params_t& params,
       size_t lower_bound, size_t upper_bound);

void calculate_table(const eqution_params_t& params,
                     const Array2D& plate_matrix, Array2D& plate_buffer);

Array2D mpi_redistribute_heat(Array2D &plate_matrix, eqution_params_t params,
                              boost::mpi::communicator& world);

template <typename T>
Magick::ColorRGB heatmap_color(const T& value, const T& min_val, const T& max_val);

void heat_equation(eqution_params_t params, Array2D& plate_matrix, boost::mpi::communicator& world,
                   const std::string& out_file);

void main_process(eqution_params_t params, Array2D& plate_matrix, const std::string& out_file,
                  boost::mpi::communicator& world, size_t matrix_height);
#endif // HEAT_EQUATION_H
