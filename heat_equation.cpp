#include "heat_equation.h"
#include <thread>
#include <mutex>
#include "concur_queue.h"
#include <Magick++.h>


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
        world.send(world.rank() + 1, 0, upper_bound);
        calculate_table(params, plate_matrix, plate_buffer);
        Array2D new_upper_bound(1, cols);
        world.recv(world.rank() + 1, 0, new_upper_bound);


        for (size_t col = 1; col < cols - 1; ++col) {
            bound_calculate(rows - 1, col, plate_matrix, plate_buffer, params, plate_matrix(rows - 2, col), new_upper_bound(0, col));
        }

    } else if (world.rank() == world.size() - 1) {
        Array2D lower_bound = plate_matrix(0, 1, 0, cols);

        world.send(world.rank() - 1, 0, lower_bound);


        calculate_table(params, plate_matrix, plate_buffer);

        Array2D new_lower_bound(1, plate_matrix.get_width());

        world.recv(world.rank() - 1, 0, new_lower_bound);

        for (size_t col = 1; col < cols - 1; ++col)
        {
            bound_calculate(0, col, plate_matrix, plate_buffer, params, new_lower_bound(0, col), plate_matrix(1, col));
        }

    }



    else {
        Array2D upper_bound = plate_matrix(rows - 1, rows,
                                           0, cols);
        world.send(world.rank() + 1, 0, upper_bound);

        Array2D lower_bound = plate_matrix(0, 1, 0, cols);

        world.send(world.rank() - 1, 0, lower_bound);


        calculate_table(params, plate_matrix, plate_buffer);

        Array2D new_upper_bound(1, cols);
        Array2D new_lower_bound(1, cols);

        world.recv(world.rank() - 1, 0, new_lower_bound);
        world.recv(world.rank() + 1, 0, new_upper_bound);
        for (size_t col = 1; col < cols - 1; ++col)
        {
            bound_calculate(rows - 1, col, plate_matrix, plate_buffer, params,
              plate_matrix(rows - 2, col), new_upper_bound(0, col));

            bound_calculate(0, col, plate_matrix, plate_buffer, params,
              new_lower_bound(0, col), plate_matrix(1, col));
        }
    }
    return plate_buffer;
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

void main_process(eqution_params_t params, Array2D& plate_matrix, const std::string& out_file,
                  boost::mpi::communicator& world,  size_t matrix_height) {
    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

    std::stringstream ss;
    ss << cols << "x" << rows;
    std::string size = ss.str();
    std::stringstream ss2;
    ss2 << cols*10 << "x" << rows*10;
    std::string scaled_size = ss.str();

    size_t com_size = world.size();

    std::vector<std::thread> recv_threads;
    concur_queue<SEGMENT> segments;

    std::vector<Magick::Image> frames;

    // merge segments

    Magick::InitializeMagick("");
    Magick::Image out_img(size.c_str(), "white");
    out_img.type(Magick::TrueColorType);
    out_img.scale(scaled_size.c_str());
    size_t all_iterations = params.max_iterations / params.iteration_to_save;
    for (size_t iteration = 0; iteration <= all_iterations; iteration++)
    {
        for (size_t i = 1; i <= com_size - 1; ++i)
        {

            size_t segment_rows_n = i == com_size - 1 ?
                        rows - (i - 1) * matrix_height : matrix_height;

            Array2D cur_segment(segment_rows_n, cols);

            world.recv(i, 0, cur_segment);

            size_t row_offset = (i - 1) * matrix_height;


            for (size_t row = 0; row < segment_rows_n; ++row)
            {
                for (size_t col = 0; col < cols; ++col)
                {
                    out_img.pixelColor(col, row + row_offset, heatmap_color(cur_segment(row, col), 0., 100.));
                }
            }


            //        frames.push_back(out_img);

        }
        out_img.write(out_file);
    }

    for (auto& item: recv_threads)
        item.join();

//    Magick::writeImages(frames.begin(), frames.end(), out_file);



}
void heat_equation(eqution_params_t params, Array2D& plate_matrix, boost::mpi::communicator& world,
                   const std::string& out_file) {
    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

    bool module = rows % (world.size() - 1);

    size_t matrix_height = rows / (world.size() - 1) + static_cast<size_t>(module != 0);

    if (world.rank() == 0) {
        main_process(params, plate_matrix, out_file, world, matrix_height);
    }
    else {
        size_t matrix_bound = world.rank() == world.size() - 1 ? rows : world.rank() * matrix_height;
        Array2D new_plate_matrix = plate_matrix((world.rank() - 1) * matrix_height, matrix_bound,
                                                0, cols);

        for (size_t i = 0; i < params.max_iterations; ++i) {
            if (i % params.iteration_to_save == 0) {
                world.send(0, 0, new_plate_matrix);
            }
            new_plate_matrix = mpi_redistribute_heat(new_plate_matrix, params, world);
        }
        world.send(0, 0, new_plate_matrix);

    }

}
