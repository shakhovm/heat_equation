#ifndef ARRAY2D_H
#define ARRAY2D_H
#include <string>
#include <iostream>
#include <boost/mpi.hpp>

class Array2D
{
public:
    double *data_;
    typedef double* iterator;
    size_t width;
    size_t height;

    void array_copy(double* to, const double *from, size_t size);

    Array2D(): data_(nullptr){};
    Array2D(size_t width, size_t height, double value = 0.);
    Array2D(const Array2D& array);

    Array2D &operator=(const Array2D& array);
    ~Array2D();

    void swap(Array2D& array);
    double &operator()(size_t width, size_t height);

    Array2D operator()(size_t row_first, size_t row_last, size_t col_first, size_t col_last) {
        Array2D array(row_last - row_first, col_last - col_first);
        for (size_t i = row_first, i2 = 0; i < row_last; ++i, ++i2) {
            for (size_t j = col_first, j2 = 0; j < col_last; ++j, ++j2) {
                array(i2, j2) = operator()(i, j);
            }
        }
        return array;
    }


    inline iterator begin() { return data_; }

    inline iterator end() { return data_ + width*height + 1 ;}

    size_t get_width() const { return width; }

    size_t get_height() const { return height; }

    inline double* data() { return data_; }

    void print() {
        for (size_t i = 0; i < width; ++i) {
            for (size_t j = 0; j < height; ++j) {
                std::cout << operator()(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & width;
        ar & height;
        for (size_t i = 0; i < width*height; ++i) {
            ar & data_[i];
        }
    }

};

#endif // ARRAY2D_H
