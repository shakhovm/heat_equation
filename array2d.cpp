#include "array2d.h"

void Array2D::array_copy(double* to, const double *from, size_t size)  {
    for (size_t i = 0; i < size; i++) {
        new (to + i) double (from[i]);
    }
}

Array2D::Array2D(size_t row, size_t col, double value) :
    data_(new double[row*col + 1]), width(row), height(col)
{
    for (size_t i = 0; i < row*col; ++i) {
        data_[i] = value;
    }
}

Array2D::Array2D(const Array2D &array) {
    data_ = new double[array.width*array.height + 1];
    array_copy(data_, array.data_, array.width*array.height);
    width = array.width;
    height = array.height;
}

Array2D &Array2D::operator=(const Array2D &array)
{
    Array2D new_array(array);
    swap(new_array);
    return *this;
}

Array2D::~Array2D() {
    delete[] data_;
}

void Array2D::swap(Array2D &array)
{
    std::swap(array.data_, data_);
    std::swap(array.width, width);
    std::swap(array.height, height);
}

double &Array2D::operator()(size_t row, size_t col)
{
    return data_[width*col + row];
}
