#ifndef ARRAY2D_H
#define ARRAY2D_H
#include <string>
#include <iostream>

class Array2D
{

    double *data_;
    typedef double* iterator;
    size_t width;
    size_t height;

    void array_copy(double* to, const double *from, size_t size);
public:
    Array2D(): data_(nullptr){};
    Array2D(size_t width, size_t height, double value = 0.);
    Array2D(const Array2D& array);
    Array2D &operator=(const Array2D& array);
    ~Array2D();

    void swap(Array2D& array);
    double &operator()(size_t width, size_t height);

    inline iterator begin() { return data_; }

    inline iterator end() { return data_ + width*height + 1 ;}

    size_t& get_width() { return width; }

    size_t& get_height() { return height; }

    inline double* data() { return data_; }

    void print() {
        for (size_t i = 0; i < width; ++i) {
            for (size_t j = 0; j < height; ++j) {
                std::cout << operator()(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
};

#endif // ARRAY2D_H
