#ifndef __MATRIX_HEADER__
#define __MATRIX_HEADER__

#include <iostream>

template <typename T>
class Matrix{

  private:
  std::size_t xsize;
  std::size_t ysize;
  T *data;

  public:

  Matrix<T>(size_t x,size_t y) {
    xsize = x;
    ysize = y;
    data = new T[x*y]();
  }

  T& operator()(size_t x, size_t y) {
    return data[ y*xsize + x ];
  }

  const T& operator()(size_t x, size_t y) const {
    return data[ y*xsize + x ];
  }

  ~Matrix(){
    delete[] data;
  }
};

#endif
