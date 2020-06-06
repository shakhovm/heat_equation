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
#include <thread>
#include "concur_queue.h"


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

Array2D file_handler(const std::string file_name) {
    size_t rows = 5, cols = 5;
    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 400, density = 8'900, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    double end_time = 2;

    std::ifstream input_stream(file_name, std::ifstream::in);
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

        throw std::runtime_error("bad initial conditions");
    }

    return plate_matrix;
}

void f() {

    size_t rows, cols;
    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 4000, density = 8'900, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    double end_time = 2;
    size_t iteration = 100;
    std::ifstream input_stream("./table.txt", std::ifstream::in);
    double buffer;
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
        return;
    }

    std::string output_path = "./../output/time_", cur_path;
    std::stringstream ss;
    ss << rows << "x" << cols;
    std::string size = ss.str();

    std::stringstream ss2;
    ss2 << rows*10 << "x" << cols*10;
    std::string scaled_size = ss.str();
    Magick::InitializeMagick("");
    for (int i = 0; i < iteration; i++)
    {
        plate_matrix = redistribute_heat(plate_matrix, delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity);
    }
    plate_matrix.print();
    Magick::Image out_img(size.c_str(), "white");
    out_img.type(Magick::TrueColorType);
    for (size_t row = 0; row < rows; ++row)
    {
        for (size_t col = 0; col < cols; ++col)
        {
            auto color = heatmap_color(plate_matrix(row, col), 0.,100.);
            out_img.pixelColor(col, row, heatmap_color(plate_matrix(row, col), 0., 100.));
        }
    }
//    out_img.scale(scaled_size.c_str());
    out_img.write("output_mpi.bmp");
}

struct eqution_params_t {
    size_t row_start,
           row_end,
           col_start,
           col_end;

    double delta_t,
           phys_params,
           delta_x_sq,
           delta_y_sq;

};

void calculate_table(const eqution_params_t& params,
                     const Array2D& plate_matrix, Array2D& plate_buffer) {

    for (size_t row = params.row_start; row < params.row_end; ++row)
    {
        for (size_t col = params.col_start; col < params.col_end; ++col)
        {
            double laplasian_x = (plate_matrix(row, col - 1) - 2 * plate_matrix(row, col) + plate_matrix(row, col + 1)) / params.delta_x_sq,
                    laplasian_y = (plate_matrix(row - 1, col) - 2 * plate_matrix(row, col) + plate_matrix(row + 1, col)) / params.delta_y_sq;
            plate_buffer(row, col) = plate_matrix(row, col) + params.delta_t * params.phys_params * (laplasian_x + laplasian_y);
        }
    }
}


Array2D mpi_redistribute_heat(Array2D &plate_matrix, eqution_params_t params,
                              boost::mpi::communicator& world)
{
    size_t rows = plate_matrix.get_height();
    size_t cols = plate_matrix.get_width();
    params.row_end = rows - 1;
    params.col_end = cols - 1;

    Array2D plate_buffer = plate_matrix;

    if (world.rank() == 1) {
        Array2D upper_bound = plate_matrix(rows - 1, rows, 0, cols);
        world.send(2, 0, upper_bound);
        calculate_table(params, plate_matrix, plate_buffer);
        Array2D new_upper_bound(1, cols);
        world.recv(2, 0, new_upper_bound);

        for (size_t col = 1; col < cols - 1; ++col) {
            double laplasian_x = (plate_matrix(rows - 1, col - 1) - 2 * plate_matrix(rows - 1, col) +
                                  plate_matrix(rows - 1, col + 1)) / params.delta_x_sq;

            double laplasian_y = (plate_matrix(rows - 2, col) - 2 * plate_matrix(rows - 1, col) +
                                  new_upper_bound(0, col)) / params.delta_y_sq;

            plate_buffer(rows - 1, col) = plate_matrix(rows - 1, col) + params.delta_t *
                    params.phys_params * (laplasian_x + laplasian_y);
        }

    } else if (world.rank() == 2) {
        Array2D lower_bound = plate_matrix(0, 1, 0, cols);
        world.send(1, 0, lower_bound);

        calculate_table(params, plate_matrix, plate_buffer);
        Array2D new_lower_bound(1, plate_matrix.get_width());

        world.recv(1, 0, new_lower_bound);

        for (size_t col = 1; col < cols - 1; ++col)
        {
            double laplasian_x = (plate_matrix(0, col - 1) - 2 * plate_matrix(0, col) + plate_matrix(0, col + 1)) / params.delta_x_sq;
            double laplasian_y = (new_lower_bound(0, col) - 2 * plate_matrix(0, col) + plate_matrix(1, col)) / params.delta_y_sq;
            plate_buffer(0, col) = plate_matrix(0, col) + params.delta_t * params.phys_params * (laplasian_x + laplasian_y);
        }

    }

    return plate_buffer;
}

eqution_params_t params_init() {
    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 400, density = 8'900, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    return eqution_params_t{1, 1, 1, 1, delta_t, phys_params, delta_x_sq, delta_y_sq};
}

void recv_segment(const size_t rank, concur_queue<SEGMENT>& seg_queue,
                    boost::mpi::communicator& world,
                    const size_t rows, const size_t cols)
{
    Array2D segment_matrix(rows, cols);
    std::cout << " waiting for recv rank " << rank << std::endl;
    world.recv(rank, 0, segment_matrix);
    seg_queue.push(SEGMENT(segment_matrix, rank));
}

void misha_function(Array2D& plate_matrix, boost::mpi::communicator& world,
                    const eqution_params_t& params) {
    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

    double end_time = 2;
    size_t iterations = 500;
    Array2D plate_matrix_first = plate_matrix(0, rows / 2, 0, cols);
    Array2D plate_matrix_second = plate_matrix(rows / 2, rows, 0, cols);
    std::stringstream ss;
    ss << cols << "x" << rows;
    std::string size = ss.str();
    std::stringstream ss2;
    ss2 << cols*10 << "x" << rows*10;
    std::string scaled_size = ss.str();

    if (world.rank() == 0)
    {
        // getting parts and merging
        size_t com_size = world.size(), cur_segment_rank;
        size_t segment_rows_n = rows / (com_size - 1);
        std::vector<std::thread> recv_threads;
        concur_queue<SEGMENT> segments;
        for (size_t com_rank = 1; com_rank < com_size; ++com_rank)
        {
            recv_threads.emplace_back(recv_segment, com_rank, std::ref(segments), std::ref(world), segment_rows_n, cols);
        }
        // merge segments
        SEGMENT cur_segment;
        size_t segments_left = com_size - 1;
        size_t row_offset;
        while (segments_left)
        {
            row_offset = cur_segment.second * segment_rows_n;
            cur_segment = segments.pop();
            std::cout << std::endl << "checkpoint; segment rank " << cur_segment.second << std::endl;
            for (size_t row = 0; row < segment_rows_n; ++row)
            {
                for (size_t col = 0; col < cols; ++col)
                {
                    plate_matrix(row + row_offset, col) = cur_segment.first(row, col);
                }
            }
            --segments_left;
        }
        for (auto& item: recv_threads)
            item.join();

        Magick::InitializeMagick("");
        Magick::Image out_img(size.c_str(), "white");
        out_img.type(Magick::TrueColorType);
        for (size_t row = 0; row < rows; ++row)
        {
            for (size_t col = 0; col < cols; ++col)
            {
                auto color = heatmap_color(plate_matrix(row, col), 0.,100.);
                out_img.pixelColor(col, row, heatmap_color(plate_matrix(row, col), 0., 100.));
            }
        }
        out_img.scale(scaled_size.c_str());
        out_img.write("./output/output_mpi_misha.bmp");

    } else if (world.rank() == 1)
    {
        for (size_t i = 0; i < iterations; ++i)
        {
            plate_matrix_first = mpi_redistribute_heat(plate_matrix_first, params, world);
        }
        world.send(0, 0, plate_matrix_first);
    } else if (world.rank() == 2)
    {
        for (size_t i = 0; i < iterations; ++i)
        {
            plate_matrix_second = mpi_redistribute_heat(plate_matrix_second, params, world);
        }
        world.send(0, 0, plate_matrix_second);
    }

    std::cout << std::endl << world.rank() << " is dying " << std::endl;
}

void sequantial_program() {
    size_t iteration = 1;
    std::cout << "seq program!" << std::endl;
    Array2D plate_matrix = file_handler("./../table.txt");

    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 400, density = 8'900, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    double end_time = 2;

    for(size_t i = 0; i < iteration; ++i) {
        plate_matrix = redistribute_heat(plate_matrix, delta_x, delta_y, delta_t,
                                               temp_conduct, density, temp_capacity);
        plate_matrix.print();
        std::cout << std::endl;
    }


}

int main(int argc, char* argv[])
{
    boost::mpi::environment env{argc, argv};
    boost::mpi::communicator world;
    eqution_params_t params = params_init();
    Array2D plate_matrix = file_handler("./table.txt");

////    if (world.rank() == 0)
////        sequantial_program();

    misha_function(plate_matrix, world, params);

//    f();
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
