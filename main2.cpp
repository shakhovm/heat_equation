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

template <typename T>
inline bool von_neumann_criterion(const T& delta_x, const T& delta_y, const T& delta_t,
            const T& conduction, const T& density, const T& capacity)
{
    if (!std::is_arithmetic<T>::value)
        throw std::invalid_argument("Non-arithmetic type passed");
    return delta_t <= 0.25 * pow(std::max(delta_x, delta_y), 2)
                      * density * capacity / conduction;
}

//template <typename T>
//Magick::ColorRGB heatmap_color(const T& value, const T& min_val, const T& max_val);

//template <typename T>
//Array2D redistribute_heat(Array2D &plate_matrix, const T& delta_x, const T& delta_y, const T& delta_t,
//                          const T& conduction, const T& density, const T& capacity)
//{
//    size_t rows = plate_matrix.get_height(),
//            cols = plate_matrix.get_width();
//    double delta_x_sq = delta_x * delta_x,
//            delta_y_sq = delta_y * delta_y,
//            phys_params = conduction / (density * capacity);
//    Array2D plate_buffer = plate_matrix;
//    for (size_t row = 1; row < rows - 1; ++row)
//    {
//        for (size_t col = 1; col < cols - 1; ++col)
//        {
//            double laplasian_x = (plate_matrix(row, col - 1) - 2 * plate_matrix(row, col) + plate_matrix(row, col + 1)) / delta_x_sq,
//                    laplasian_y = (plate_matrix(row - 1, col) - 2 * plate_matrix(row, col) + plate_matrix(row + 1, col)) / delta_y_sq;
//            plate_buffer(row, col) = plate_matrix(row, col) + delta_t * phys_params * (laplasian_x + laplasian_y);
//        }
//    }
//    return plate_buffer;
//}

Array2D file_handler(const std::string& file_name) {
    size_t rows, cols;

    std::ifstream input_stream(file_name, std::ifstream::in);
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
//    input_stream.close();
//    if (!von_neumann_criterion(delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity))
//    {

//        throw std::runtime_error("bad initial conditions");
//    }

    return plate_matrix;
}

//void f() {

//    size_t rows, cols;
//    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
//           temp_conduct = 400, density = 8'900, temp_capacity = 1;
//    double delta_x_sq = delta_x * delta_x,
//           delta_y_sq = delta_y * delta_y,
//           phys_params = temp_conduct / (density * temp_capacity);
//    double end_time = 2;
//    size_t iteration = 1;
//    std::ifstream input_stream("../table4.txt", std::ifstream::in);
//    double buffer;
//    input_stream >> rows >> cols;
//    Array2D plate_matrix(cols, rows), plate_buffer;

//    for (size_t row = 0; row < rows; ++row)
//    {
//        for (size_t col = 0; col < cols; ++col)
//        {
//            input_stream >> buffer;
//            plate_matrix(row, col) = buffer;
//        }
//    }
//    input_stream.close();
//    if (!von_neumann_criterion(delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity))
//    {
//        std::cout << "bad initial conditions" << std::endl;
//        return;
//    }

//    std::string output_path = "./../output/time_", cur_path;
//    std::stringstream ss;
//    ss << rows << "x" << cols;
//    std::string size = ss.str();

//    std::stringstream ss2;
//    ss2 << rows*10 << "x" << cols*10;
//    std::string scaled_size = ss.str();
//    Magick::InitializeMagick("");
//    for (int i = 0; i < iteration; i++)
//    {
//        plate_matrix = redistribute_heat(plate_matrix, delta_x, delta_y, delta_t, temp_conduct, density, temp_capacity);
//    }
////    plate_matrix.print();
//    Magick::Image out_img(size.c_str(), "white");
//    out_img.type(Magick::TrueColorType);
//    for (size_t row = 0; row < rows; ++row)
//    {
//        for (size_t col = 0; col < cols; ++col)
//        {
//            auto color = heatmap_color(plate_matrix(row, col), 0.,100.);
//            out_img.pixelColor(col, row, heatmap_color(plate_matrix(row, col), 0., 100.));
//        }
//    }
////    out_img.scale(scaled_size.c_str());
//    out_img.write("../output/output_mpi.bmp");
//}


eqution_params_t params_init() {
    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
           temp_conduct = 400, density = 8'900, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    return eqution_params_t{1, 1, 1, 1, delta_t, phys_params, delta_x_sq, delta_y_sq};
}

//void recv_segment(const size_t rank, concur_queue<SEGMENT>& seg_queue,
//                    boost::mpi::communicator& world,
//                    const size_t rows, const size_t cols, std::mutex& m)
//{
//    Array2D segment_matrix(rows, cols);
//    {
//        std::lock_guard lg{m};
//        world.recv(rank, 0, segment_matrix);
//    }
//    seg_queue.push(SEGMENT(segment_matrix, rank));
//}

//void misha_function(Array2D& plate_matrix, boost::mpi::communicator& world,
//                    const eqution_params_t& params) {


//    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

//    double end_time = 2;

//    size_t iterations = 1000;
//    std::stringstream ss;
//    ss << cols << "x" << rows;
//    std::string size = ss.str();
//    std::stringstream ss2;
//    ss2 << cols*10 << "x" << rows*10;
//    std::string scaled_size = ss.str();

//    bool module = rows % (world.size() - 1);
//    size_t matrix_height = rows / (world.size() - 1) + static_cast<size_t>(module != 0);

//    if (world.rank() == 0) {

////        for (size_t i = 0; i < world.size(); ++i) {
////            Array2D part_to_send = plate_matrix(i * matrix_height, (i + 1)*matrix_height, 0, cols);;
////        }
//        std::mutex m;
//        size_t com_size = world.size(), cur_segment_rank;
//        size_t segment_rows_n = rows / (com_size - 1);
//        std::vector<std::thread> recv_threads;
//        concur_queue<SEGMENT> segments;
//        for (size_t com_rank = 1; com_rank < com_size; ++com_rank)
//        {
//            segment_rows_n = com_rank == world.size() - 1 ?
//                        rows - (com_rank - 1)*matrix_height : matrix_height;
////            segment_rows_n = matrix_height;
//            std::cout << segment_rows_n << std::endl;
//            recv_threads.emplace_back(recv_segment, com_rank, std::ref(segments), std::ref(world), segment_rows_n, cols,
//                                      std::ref(m));
//        }
////        return;
//        // merge segments
//        SEGMENT cur_segment;
//        size_t segments_left = com_size - 1;
//        size_t row_offset;
//        while (segments_left > 0)
//        {

//          cur_segment = segments.pop();
////          row_offset = cur_segment.second == world.size() - 1 ?
////                      rows - (cur_segment.second - 1) * matrix_height : (cur_segment.second - 1) * segment_rows_n;
//          row_offset = (cur_segment.second - 1)*matrix_height;
//          segment_rows_n = cur_segment.second == world.size() - 1 ?
//                      rows - (cur_segment.second - 1) * matrix_height : matrix_height;
////          row_offset = (cur_segment.second - 1) * segment_rows_n;
//          std::cout << std::endl << "checkpoint; segment rank " << cur_segment.second << std::endl;
//          std::cout << "OFFSET" << row_offset << std::endl;

//          for (size_t row = 0; row < segment_rows_n; ++row)
//          {
//              for (size_t col = 0; col < cols; ++col)
//              {
//                  plate_matrix(row + row_offset, col) = cur_segment.first(row, col);
//              }
//          }
//          --segments_left;
//          std::cout << "WEOPKDP" << std::endl;
//        }
//        for (auto& item: recv_threads)
//            item.join();
//        std::cout << "JOINT";

//        Magick::InitializeMagick("");
//        Magick::Image out_img(size.c_str(), "white");
//        out_img.type(Magick::TrueColorType);
//          for (size_t row = 0; row < rows; ++row)
//          {
//              for (size_t col = 0; col < cols; ++col)
//              {
//                  auto color = heatmap_color(plate_matrix(row, col), 0.,100.);
//                  out_img.pixelColor(col, row, heatmap_color(plate_matrix(row, col), 0., 100.));
//              }
//          }
//          out_img.scale(scaled_size.c_str());
//          out_img.write("../output/output_mpi_misha.bmp");
//    }
//    else {
//        size_t matrix_bound = world.rank() == world.size() - 1 ? rows : world.rank() * matrix_height;
//        Array2D new_plate_matrix = plate_matrix((world.rank() - 1) * matrix_height, matrix_bound,
//                                                0, cols);
//        for (size_t i = 0; i < iterations; ++i) {
//            new_plate_matrix = mpi_redistribute_heat(new_plate_matrix, params, world);

//        }
////        new_plate_matrix.print();
//        world.send(0, 0, new_plate_matrix);
//        std::cout << "SENT!" << std::endl;
//    }


//}

////void sequantial_program(Array2D& plate_matrix) {
////    size_t iteration = 1;

////    double delta_t = 0.005, delta_x = 0.1, delta_y = 0.1,
////           temp_conduct = 400, density = 8'900, temp_capacity = 1;
////    double delta_x_sq = delta_x * delta_x,
////           delta_y_sq = delta_y * delta_y,
////           phys_params = temp_conduct / (density * temp_capacity);
////    double end_time = 2;

////    for(size_t i = 0; i < iteration; ++i) {
////        plate_matrix = redistribute_heat(plate_matrix, delta_x, delta_y, delta_t,
////                                               temp_conduct, density, temp_capacity);
////        plate_matrix.print();
////        std::cout << std::endl;
////    }


////}


int main(int argc, char* argv[])
{
    boost::mpi::environment env{argc, argv};
    boost::mpi::communicator world;
    eqution_params_t params = params_init();
    Array2D plate_matrix = file_handler("../table4.txt");
    if (world.size() == 1)
        return 0;
    else
        heat_equation(params, plate_matrix, world, 1000, "../output/x.gif");
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
