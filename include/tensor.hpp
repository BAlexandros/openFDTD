#ifndef __TENSOR_HEADER__
#define __TENSOR_HEADER__

#include <iostream>

template <typename T>
class Tensor{

  private:
  std::size_t xsize;
  std::size_t ysize;
  std::size_t zsize;
  T *data;

  public:

  Tensor<T>(size_t x,size_t y,size_t z){
    xsize = x;
    ysize = y;
    zsize = z;
    data = new T[x*y*z]();
  }

  T& operator()(size_t x, size_t y, size_t z){
    return data[ (z*ysize + y)*xsize + x];
  }

  const T& operator()(size_t x, size_t y, size_t z) const {
    return data[ (z*ysize + y)*xsize + x];
  }

  ~Tensor(){
    delete[] data;
  }
};

#endif
