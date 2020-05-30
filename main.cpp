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

template <typename T>
inline bool von_neumann_criterion(const T& delta_x, const T& delta_y, const T& delta_t,
            const T& conduction, const T& density, const T& capacity)
{
    if (!std::is_arithmetic<T>::value)
        throw std::invalid_argument("Non-arithmetic type passed");
    return delta_t <= 0.25 * pow(std::max(delta_x, delta_y), 2)
                      * density * capacity / conduction;
}

template <typename T>
Magick::ColorRGB heatmap_color(const T& value, const T& min_val, const T& max_val);

template <typename T>
Array2D redistribute_heat(Array2D &plate_matrix, const T& delta_x, const T& delta_y, const T& delta_t,
                          const T& conduction, const T& density, const T& capacity)
{
    size_t rows = plate_matrix.get_height(),
            cols = plate_matrix.get_width();
    double delta_x_sq = delta_x * delta_x,
            delta_y_sq = delta_y * delta_y,
            phys_params = conduction / (density * capacity);
    Array2D plate_buffer = plate_matrix;
    for (size_t row = 1; row < rows - 1; ++row)
    {
        for (size_t col = 1; col < cols - 1; ++col)
        {
            double laplasian_x = (plate_matrix(row, col - 1) - 2 * plate_matrix(row, col) + plate_matrix(row, col + 1)) / delta_x_sq,
                    laplasian_y = (plate_matrix(row - 1, col) - 2 * plate_matrix(row, col) + plate_matrix(row + 1, col)) / delta_y_sq;
            plate_buffer(row, col) = plate_matrix(row, col) + delta_t * phys_params * (laplasian_x + laplasian_y);
        }
    }
    return plate_buffer;
}

int main(int argc, char* argv[])
{
    size_t rows = 5, cols = 5;
    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 400, density = 8'900, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    double end_time = 2;

    std::ifstream input_stream("./../../table.txt", std::ifstream::in);
    size_t buffer;
    input_stream >> rows >> cols;
    Array2D plate_matrix(cols, rows), plate_buffer;
    for (size_t row = 0; row < rows; ++row)
    {
        for (size_t col = 0; col < cols; ++col)
        {
            input_stream >> buffer;
            plate_matrix(row, col) = buffer;
        }
    }
    input_stream.close();
    if (!von_neumann_criterion(delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity))
    {
        std::cout << "bad initial conditions" << std::endl;
        return 1;
    }

    std::string output_path = "./output/time_", cur_path;
    for (int i = 0; i < 2000; i++)
    {
        Magick::InitializeMagick("");
        Magick::Image out_img(Magick::Geometry(cols, rows), "white");
        out_img.type(Magick::TrueColorType);
        for (size_t row = 0; row < rows; ++row) {
            for (size_t col = 0; col < cols; ++col) {
                if (i == 99)
                {
                    auto c = heatmap_color(plate_matrix(row, col),0., 100.);
                }
                out_img.pixelColor(col, row, heatmap_color(plate_matrix(row, col), 0., 100.));
            }
        }
        cur_path = output_path + std::to_string(i) + ".bmp";
        if (i % 25 == 0)
        {
            Magick::Image ttt = out_img;
            ttt.scale(Magick::Geometry(cols * 10, rows * 10));

            ttt.write(cur_path);
        }

        plate_matrix = redistribute_heat(plate_matrix, delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity);
    }
}

template <typename T>
Magick::ColorRGB heatmap_color(const T& value, const T& min_val, const T& max_val)
{
    double red, green, blue;
    double ratio = 2 * static_cast<double>(value - min_val) / (max_val - min_val);
    blue = std::max(0., 1 - ratio);
    red = std::max(0., ratio - 1);
    green = 1 - blue - red;
    return Magick::ColorRGB(red, green, blue);
}
