#ifndef HEAT_EQUATION_FUNCTIONS_H
#define HEAT_EQUATION_FUNCTIONS_H

#include "../array2d/array2d.h"
#include <boost/mpi.hpp>
#include <Magick++.h>

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
       double lower_bound, double upper_bound);

void calculate_table(const eqution_params_t& params,
                     const Array2D& plate_matrix, Array2D& plate_buffer);

Array2D mpi_redistribute_heat(Array2D &plate_matrix, eqution_params_t params,
                              boost::mpi::communicator& world);


Array2D file_handler(const std::string& file_name);

void copy_vertical_bounds(Array2D &to, const Array2D &from);


#endif // HEAT_EQUATION_FUNCTIONS_H
