#include "../../include/heat_equation/heat_equation.h"

Array2D mpi_redistribute_heat(Array2D &plate_matrix, eqution_params_t params,
                              boost::mpi::communicator& world)
{
    size_t rows = plate_matrix.get_height();
    size_t cols = plate_matrix.get_width();
    params.row_end = rows - 1;
    params.col_end = cols - 1;

    Array2D plate_buffer(rows, cols);
    copy_vertical_bounds(plate_buffer, plate_matrix);
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


        Array2D new_lower_bound(1, plate_matrix.get_width());
        world.recv(world.rank() - 1, 0, new_lower_bound);

        Array2D lower_bound = plate_matrix(0, 1, 0, cols);

        world.send(world.rank() - 1, 0, lower_bound);


        calculate_table(params, plate_matrix, plate_buffer);


        for (size_t col = 1; col < cols - 1; ++col)
        {
            bound_calculate(0, col, plate_matrix, plate_buffer, params, new_lower_bound(0, col), plate_matrix(1, col));
        }

    }

    else {

        Array2D new_upper_bound(1, cols);
        Array2D new_lower_bound(1, cols);

        world.recv(world.rank() - 1, 0, new_lower_bound);

        Array2D upper_bound = plate_matrix(rows - 1, rows,
                                           0, cols);
        world.send(world.rank() + 1, 0, upper_bound);



        Array2D lower_bound = plate_matrix(0, 1, 0, cols);

        world.send(world.rank() - 1, 0, lower_bound);

        world.recv(world.rank() + 1, 0, new_upper_bound);



        calculate_table(params, plate_matrix, plate_buffer);

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


void main_process(eqution_params_t params, Array2D& plate_matrix, const std::string& out_file,
                  boost::mpi::communicator& world,  size_t matrix_height) {
    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

    std::stringstream ss;
    ss << cols << "x" << rows;
    std::string size = ss.str();

    size_t com_size = world.size();

    std::vector<Magick::Image> frames;

    // merge segments

    Magick::InitializeMagick("");
    Magick::Image out_img(size.c_str(), "white");

    out_img.type(Magick::TrueColorType);

    size_t all_iterations = params.max_iterations / params.iteration_to_save;
    for (size_t iteration = 0; iteration <= all_iterations; iteration++)
    {
        for (size_t i = 1; i < com_size; ++i)
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

        }
        out_img.animationDelay(50);
        frames.push_back(out_img);
//        out_img.write("out_file" + std::to_string(iteration) + ".bmp");
    }
    Magick::writeImages(frames.begin(), frames.end(), out_file);
}

void one_process(eqution_params_t params, const std::string& matrix_file, const std::string& out_file) {
    Array2D plate_matrix = file_handler(matrix_file);
    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

    std::stringstream ss;
    ss << cols << "x" << rows;
    std::string size = ss.str();

    std::vector<Magick::Image> frames;


    Magick::InitializeMagick("");
    Magick::Image out_img(size.c_str(), "white");

    out_img.type(Magick::TrueColorType);

    for (size_t iteration = 0; iteration < params.max_iterations; ++iteration) {
        if (iteration % params.iteration_to_save == 0) {
            for (size_t row = 0; row < rows; ++row)
            {
                for (size_t col = 0; col < cols; ++col)
                {
                    out_img.pixelColor(col, row, heatmap_color(plate_matrix(row, col), 0., 100.));
                }
            }
            out_img.animationDelay(50);
            frames.push_back(out_img);
//            out_img.write("out_file" + std::to_string(iteration) + ".bmp");

        }
        Array2D plate_buffer(rows, cols);
        copy_vertical_bounds(plate_buffer, plate_matrix);
        calculate_table(params, plate_matrix, plate_buffer);
        plate_matrix = plate_buffer;

    }
    Magick::writeImages(frames.begin(), frames.end(), out_file);
}

void two_processes(eqution_params_t params, const std::string& matrix_file, const std::string& out_file,
                   boost::mpi::communicator& world) {
    Array2D plate_matrix = file_handler(matrix_file);
    size_t rows = plate_matrix.get_height(), cols = plate_matrix.get_width();

    std::stringstream ss;
    ss << cols << "x" << rows;
    std::string size = ss.str();

    std::vector<Magick::Image> frames;


    Magick::InitializeMagick("");
    Magick::Image out_img(size.c_str(), "white");

    out_img.type(Magick::TrueColorType);

    if (world.rank() == 1) {
        for (size_t iteration = 0; iteration < params.max_iterations; ++iteration) {
            if (iteration % params.iteration_to_save == 0) {
                world.send(0, 0, plate_matrix);
            }
            Array2D plate_buffer(rows, cols);
            copy_vertical_bounds(plate_buffer, plate_matrix);
            calculate_table(params, plate_matrix, plate_buffer);
            plate_matrix = plate_buffer;

        }
        world.send(0, 0, plate_matrix);
    }
    else if (world.rank() == 0)
    {
        size_t all_iterations = params.max_iterations / params.iteration_to_save;
        for (size_t iteration = 0; iteration <= all_iterations; iteration++)
        {
            Array2D cur_segment(rows, cols);
            world.recv(1, 0, cur_segment);
            for (size_t row = 0; row < rows; ++row)
            {
                for (size_t col = 0; col < cols; ++col)
                {
                    out_img.pixelColor(col, row, heatmap_color(cur_segment(row, col), 0., 100.));
                }
            }
            out_img.animationDelay(50);
            frames.push_back(out_img);
//            out_img.write("out_file" + std::to_string(iteration) + ".bmp");

        }
        Magick::writeImages(frames.begin(), frames.end(), out_file);
    }
}



void heat_equation(eqution_params_t params, const std::string& matrix_file, boost::mpi::communicator& world,
                   const std::string& out_file, size_t rows, size_t cols) {

    bool module = rows % (world.size() - 1);

    size_t matrix_height = rows / (world.size() - 1) + static_cast<size_t>(module != 0);

    if (world.rank() == 0) {
        Array2D plate_matrix = file_handler(matrix_file);


        for (size_t i = 1; i < world.size(); ++i) {
            size_t matrix_bound = i == world.size() - 1 ? rows : i * matrix_height;
            Array2D new_plate_matrix = plate_matrix((i - 1) * matrix_height, matrix_bound,
                                                    0, cols);
            world.send(i, 0, new_plate_matrix);
        }
        main_process(params, plate_matrix, out_file, world, matrix_height);
    }
    else {
        size_t segment_rows_n = world.rank() == world.size() - 1 ?
                    rows - (world.rank() - 1) * matrix_height : matrix_height;

        Array2D new_plate_matrix(segment_rows_n, cols);
        world.recv(0, 0, new_plate_matrix);
        for (size_t i = 0; i < params.max_iterations; ++i) {
            if (i % params.iteration_to_save == 0) {
                world.send(0, 0, new_plate_matrix);
            }
            new_plate_matrix = mpi_redistribute_heat(new_plate_matrix, params, world);
        }
        world.send(0, 0, new_plate_matrix);

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




