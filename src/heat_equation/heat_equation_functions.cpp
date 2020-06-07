#include "../../include/heat_equation/heat_equation_functions.h"
#include <fstream>

void copy_vertical_bounds(Array2D &to, const Array2D &from)
{
    bool is_flat_v = to.get_width() == 1, is_flat_h = to.get_height() == 1;
    size_t rows = from.get_height(), cols = from.get_width();
    size_t last_row = from.get_height() - 1, last_col = from.get_width() - 1;
    for (size_t row = 0; row < rows; ++row)
    {
        to(row, 0) = from(row, 0);
        if (!is_flat_v)
            to(row, last_col) = from(row, last_col);
    }
    for (size_t col = 0; col < cols; ++col)
    {
        to(0, col) = from(0, col);
        if (!is_flat_h)
            to(last_row, col) = from(last_row, col);
    }
}

void bound_calculate(size_t row, size_t col, const Array2D& plate_matrix, Array2D& plate_buffer, const eqution_params_t& params,
       double lower_bound, double upper_bound) {
    double laplasian_x = (plate_matrix(row, col - 1) - 2 * plate_matrix(row, col) +
                          plate_matrix(row, col + 1)) / params.delta_x_sq;

    double laplasian_y = (lower_bound - 2 * plate_matrix(row, col) +
                          upper_bound) / params.delta_y_sq;

    plate_buffer(row, col) = plate_matrix(row, col) + params.delta_t *
            params.phys_params * (laplasian_x + laplasian_y);
}


void calculate_table(const eqution_params_t& params,
                     const Array2D& plate_matrix, Array2D& plate_buffer) {

    for (size_t row = params.row_start; row < params.row_end; ++row)
    {
        for (size_t col = params.col_start; col < params.col_end; ++col)
        {
            bound_calculate(row, col, plate_matrix, plate_buffer, params, plate_matrix(row - 1, col), plate_matrix(row + 1, col));
        }
    }
}



Array2D file_handler(const std::string& file_name) {
    size_t rows, cols;

    std::ifstream input_stream(file_name);
    if (!input_stream.is_open()) {
        throw std::runtime_error("No File with matrix found");
    }
    double buffer;
    input_stream >> rows >> cols;
    Array2D plate_matrix(rows, cols);
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
