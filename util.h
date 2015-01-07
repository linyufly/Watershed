// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef UTIL_H_
#define UTIL_H_

#include <cstring>

void report_error(const char *format, ...);

template<class T>
T **create_matrix(int num_of_rows, int num_of_cols) {
  T *data_array = new T[num_of_rows * num_of_cols];
  T **matrix = new T *[num_of_rows];

  for (int i = 0; i < num_of_rows; i++) {
    matrix[i] = data_array + i * num_of_cols;
  }

  return matrix;
}

template<class T>
void delete_matrix(T **matrix) {
  delete [] matrix[0];
  delete [] matrix;
}

template<class T>
T ***create_3d_array(int n1, int n2, int n3) {
  T *arr_123 = new T[n1 * n2 * n3];

  T **arr_12 = new T *[n1 * n2];
  for (int i = 0; i < n1 * n2; i++) {
    arr_12[i] = &arr_123[i * n3];
  }

  T ***arr_1 = new T **[n1];
  for (int i = 0; i < n1; i++) {
    arr_1[i] = &arr_12[i * n2];
  }

  return arr_1;
}

template<class T>
T ****create_4d_array(int n1, int n2, int n3, int n4) {
  T *arr_1234 = new T[n1 * n2 * n3 * n4];

  T **arr_123 = new T *[n1 * n2 * n3];
  for (int i = 0; i < n1 * n2 * n3; i++) {
    arr_123[i] = &arr_1234[i * n4];
  }

  T ***arr_12 = new T **[n1 * n2];
  for (int i = 0; i < n1 * n2; i++) {
    arr_12[i] = &arr_123[i * n3];
  }

  T ****arr_1 = new T ***[n1];
  for (int i = 0; i < n1; i++) {
    arr_1[i] = &arr_12[i * n2];
  }

  return arr_1;
}

template<class T>
void delete_3d_array(T ***array) {
  delete [] array[0][0];
  delete [] array[0];
  delete [] array;
}

template<class T>
void delete_4d_array(T ****array) {
  delete [] array[0][0][0];
  delete [] array[0][0];
  delete [] array[0];
  delete [] array;
}

#endif  // UTIL_H_
