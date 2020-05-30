#include <boost/mpi.hpp>
#include <iostream>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <exception>
#include <jpeglib.h>
#include "array2d.h"

template <typename T>
inline bool von_neumann_criterion(const T& delta_x, const T& delta_y, const T& delta_t,
            const T& conduction, const T& density, const T& capacity)
{
    if (!std::is_arithmetic<T>::value)
        throw std::invalid_argument("Non-arithmetic type passed");
    return delta_t <= 0.25 * pow(std::max(delta_x, delta_y), 2)
                      * density * capacity / conduction;
}

Array2D jpeg_colormap(const Array2D &heat_matrix, double max_value);
template <typename T>
std::tuple<JSAMPLE, JSAMPLE, JSAMPLE> heatmap_color(const T& value,
        const T& min_val, const T& max_val);

int main(int argc, char* argv[])
{
    size_t rows = 5, cols = 5;
    double delta_t = 100.1, delta_x = 1, delta_y = 1,
           temp_conduct = 74, density = 2'700, temp_capacity = 0.46;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    double end_time = 2;

    Array2D plate_matrix(cols, rows);
    Array2D plate_buffer;
    plate_matrix(0, 1) = 100;



    auto color = heatmap_color(0, 0, 100);
    std::cout << +std::get<0>(color) << " " << +std::get<1>(color) << " " << +std::get<2>(color) << std::endl;




    if (!von_neumann_criterion(delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity))
    {
        std::cout << "von Neumann criterion is not satisfied" << std::endl;
        return 0;
    }

    for (double cur_time = 0; cur_time < end_time; cur_time += delta_t)
    {
        plate_buffer = plate_matrix;
        for (size_t row = 1; row < rows - 1; ++row)
        {
            for (size_t col = 1; col < cols - 1; ++col)
            {
                double laplasian_x = (plate_buffer(row, col - 1) - 2 * plate_buffer(row, col) + plate_buffer(row, col + 1)) / delta_x_sq,
                        laplasian_y = (plate_buffer(row - 1, col) - 2 * plate_buffer(row, col) + plate_buffer(row + 1, col)) / delta_y_sq;
                plate_matrix(row, col) = plate_buffer(row, col) + delta_t * phys_params * (laplasian_x + laplasian_y);
            }
        }
        plate_matrix.print();
    }
}

template <typename T>
std::tuple<JSAMPLE, JSAMPLE, JSAMPLE> heatmap_color(const T& value,
        const T& min_val, const T& max_val)
{
    JSAMPLE red, green, blue;
    double ratio = 2 * static_cast<double>(value - min_val) / (max_val - min_val);
    blue = static_cast<JSAMPLE>(std::max(0., 255 * (1 - ratio)));
    red = static_cast<JSAMPLE>(std::max(0., 255 * (ratio - 1)));
    green = 255 - blue - red;
    return std::make_tuple(red, green, blue);
}
