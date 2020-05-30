#include <boost/mpi.hpp>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
    size_t rows = 5, cols = 5;
    double delta_t = 1, delta_x = 10, delta_y = 10,
           temp_conduct = 1, density = 1, temp_capacity = 1;
    double delta_x_sq = delta_x * delta_x,
           delta_y_sq = delta_y * delta_y,
           phys_params = temp_conduct / (density * temp_capacity);
    double end_time = 3;

    //    Array2D plate_matrix(rows, cols);
//    Array2D plate_buffer;
    std::vector<std::vector<double>> plate_matrix(rows, std::vector<double>(cols, 0)), plate_buffer;
    plate_matrix[0][1] = 10;


    // TODO: von Newman condition
    for (size_t cur_time = 0; cur_time < end_time; cur_time += delta_t)
    {
        plate_buffer = plate_matrix;
        for (size_t row = 1; row < rows - 1; ++row)
        {
            for (size_t col = 1; col < cols - 1; ++col)
            {
                double laplasian_x = (plate_buffer[row][col - 1] - 2 * plate_buffer[row][col] + plate_buffer[row][col + 1]) / delta_x_sq,
                        laplasian_y = (plate_buffer[row - 1][col] - 2 * plate_buffer[row][col] + plate_buffer[row + 1][col]) / delta_y_sq;
                plate_matrix[row][col] = plate_buffer[row][col] + delta_t * phys_params * (laplasian_x + laplasian_y);
            }
        }
        std::cout << std::endl;
        for (int i = 0; i < rows; i++)
        {
          for (int j = 0; j < cols; j++)
            std::cout << plate_matrix[i][j] << " ";
          std::cout << std::endl;
        }


    }
}
