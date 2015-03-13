// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_extractor.h"

#include "util.h"

#include <cmath>

#include <algorithm>
#include <map>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

const int kConnectivity = 26;
const int kLabellingConnectivity = 26;
const double kZeroTolerance = 1e-8;

const int kDirections[26][3] = {
    {0, 0, 1}, {0, 0, -1},
    {0, 1, 0}, {0, -1, 0},
    {1, 0, 0}, {-1, 0, 0},

    {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1},
    {1, 0, 1}, {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
    {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},

    {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
    {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};

struct GridPointData {
  int x_, y_, z_, dist_;
  double height_;

  GridPointData() {}

  GridPointData(int x, int y, int z, int dist, double height)
      : x_(x), y_(y), z_(z), dist_(dist), height_(height) {}

  bool operator < (const GridPointData &other) const {
    if (height_ != other.height_) {
      return height_ < other.height_;
    }

    return dist_ < other.dist_;
  }
};

bool outside(int x, int y, int z, int nx, int ny, int nz) {
  return x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz;
}

void calculate_plateau_distance(double ***scalar_field,
                                int nx, int ny, int nz, int ***dist) {
  int *queue_x = new int[nx * ny * nz];
  int *queue_y = new int[nx * ny * nz];
  int *queue_z = new int[nx * ny * nz];

  int head = 0, tail = -1;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        dist[x][y][z] = -1;

        for (int k = 0; k < kConnectivity; k++) {
          int next_x = x + kDirections[k][0];
          int next_y = y + kDirections[k][1];
          int next_z = z + kDirections[k][2];

          if (outside(next_x, next_y, next_z, nx, ny, nz)) {
            continue;
          }

          if (scalar_field[next_x][next_y][next_z] < scalar_field[x][y][z]) {
            dist[x][y][z] = 0;
            break;
          }
        }

        if (dist[x][y][z] == 0) {
          tail++;
          queue_x[tail] = x;
          queue_y[tail] = y;
          queue_z[tail] = z;
        }
      }
    }
  }

  while (head <= tail) {
    int x = queue_x[head];
    int y = queue_y[head];
    int z = queue_z[head];
    head++;

    for (int k = 0; k < kConnectivity; k++) {
      int next_x = x + kDirections[k][0];
      int next_y = y + kDirections[k][1];
      int next_z = z + kDirections[k][2];

      if (outside(next_x, next_y, next_z, nx, ny, nz)) {
        continue;
      }

      if (dist[next_x][next_y][next_z] == -1
          && scalar_field[next_x][next_y][next_z] == scalar_field[x][y][z]) {
        dist[next_x][next_y][next_z] = dist[x][y][z] + 1;
        tail++;
        queue_x[tail] = next_x;
        queue_y[tail] = next_y;
        queue_z[tail] = next_z;
      }
    }
  }

  delete [] queue_x;
  delete [] queue_y;
  delete [] queue_z;
}

void check_waterproof(int ***basin_index, int nx, int ny, int nz) {
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        if (basin_index[x][y][z] != -1) {
          for (int k = 0; k < kConnectivity; k++) {
            int next_x = x + kDirections[k][0];
            int next_y = y + kDirections[k][1];
            int next_z = z + kDirections[k][2];

            if (outside(next_x, next_y, next_z, nx, ny, nz)) {
              continue;
            }

            if (basin_index[next_x][next_y][next_z] != -1
                && basin_index[next_x][next_y][next_z]
                   != basin_index[x][y][z]) {
              printf("Found non-waterproof.\n");
              return;
            }
          }
        }
      }
    }
  }

  printf("The watershed result is waterproof.\n");
}

double get_value(double ***scalars, const int *index) {
  return scalars[index[0]][index[1]][index[2]];
}

double central_differencing_partial(
    double ***scalars, const double *spacing, const int *dimensions,
    const int *index, int component) {
  int index_1[] = {index[0], index[1], index[2]};
  int index_2[] = {index[0], index[1], index[2]};

  if (index_1[component] > 0) {
    index_1[component]--;
  }

  if (index_2[component] + 1 < dimensions[component]) {
    index_2[component]++;
  }

  return (get_value(scalars, index_2)
         - get_value(scalars, index_1))
         / (spacing[component] * (index_2[component] - index_1[component]));
}

double square(double value) {
  return value * value;
}

double norm_3d(double *vec) {
  return sqrt(square(vec[0]) + square(vec[1]) + square(vec[2]));
}

void normalize_3d(double *coordinates) {
  double norm = norm_3d(coordinates);
  if (norm < kZeroTolerance) {
    coordinates[0] = coordinates[1] = coordinates[2] = 0.0;
    return;
  }

  coordinates[0] /= norm;
  coordinates[1] /= norm;
  coordinates[2] /= norm;
}

double dot_product_3d(double *v_1, double *v_2) {
  return v_1[0] * v_2[0] + v_1[1] * v_2[1] + v_1[2] * v_2[2];
}

void remove_seams_by_gradients(
    double ***scalar_field, int nx, int ny, int nz, double *spacing,
    int ***basin_index, int ***dist_2_valley) {
  int dimensions[] = {nx, ny, nz};

  while (true) {
    bool flag = false;

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          if (basin_index[x][y][z] == -1) {
            double grad[3];
            int index[] = {x, y, z};
            for (int component = 0; component < 3; component++) {
              grad[component] = central_differencing_partial(
                  scalar_field, spacing, dimensions, index, component);
            }
            normalize_3d(grad);

            double max_dot = -2.0;
            int label = -1;

            for (int k = 0; k < kLabellingConnectivity; k++) {
              int next_x = x + kDirections[k][0];
              int next_y = y + kDirections[k][1];
              int next_z = z + kDirections[k][2];

              if (outside(next_x, next_y, next_z, nx, ny, nz)) {
                continue;
              }

              // if (scalar_field[next_x][next_y][next_z]
              //     > scalar_field[x][y][z]) {
              //   continue;
              // }

              if (basin_index[next_x][next_y][next_z] < 0) {
                continue;
              }

              int next_index[] = {next_x, next_y, next_z};
              double next_grad[3];
              for (int component = 0; component < 3; component++) {
                next_grad[component] = central_differencing_partial(
                    scalar_field, spacing, dimensions, next_index, component);
              }
              normalize_3d(next_grad);
              double curr_dot = dot_product_3d(grad, next_grad);
              if (curr_dot > max_dot) {
                max_dot = curr_dot;
                label = basin_index[next_x][next_y][next_z];
              }
            }

            if (label == -1) {
              continue;
            }
           
            // Label the ridge point.
            basin_index[x][y][z] = -2 - label;

            // Update dist_2_valley.
            for (int k = 0; k < kConnectivity; k++) {
              int next_x = x + kDirections[k][0];
              int next_y = y + kDirections[k][1];
              int next_z = z + kDirections[k][2];

              if (outside(next_x, next_y, next_z, nx, ny, nz)) {
                continue;
              }

              if (scalar_field[next_x][next_y][next_z]
                  > scalar_field[x][y][z]) {
                continue;
              }

              if (basin_index[next_x][next_y][next_z] != label) {
                continue;
              }

              if (dist_2_valley[x][y][z] == -1
                  || dist_2_valley[x][y][z]
                     > dist_2_valley[next_x][next_y][next_z] + 1) {
                dist_2_valley[x][y][z] =
                    dist_2_valley[next_x][next_y][next_z] + 1;
              }
            }

            flag = true;
          }
        }
      }
    }

    if (!flag) {
      break;
    }

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          if (basin_index[x][y][z] < -1) {
            basin_index[x][y][z] = -(basin_index[x][y][z] + 2);
          }
        }
      }
    }
  }
}

void remove_seams(double ***scalar_field, int nx, int ny, int nz,
                  int ***basin_index, int ***dist_2_valley) {
  while (true) {
    bool flag = false;

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          /// DEBUG ///
          bool debug = x == 0 && y == 15 && z == 176;

          if (basin_index[x][y][z] == -1) {
            std::map<int, int> labels;
            for (int k = 0; k < kLabellingConnectivity; k++) {
              int next_x = x + kDirections[k][0];
              int next_y = y + kDirections[k][1];
              int next_z = z + kDirections[k][2];

              if (outside(next_x, next_y, next_z, nx, ny, nz)) {
                continue;
              }

              /// DEBUG ///
              // if (debug) {
              //   printf("next: %d, %d, %d\n", next_x, next_y, next_z);
              //   printf("scalar_field[next_x][next_y][next_z] = %lf\n",
              //          scalar_field[next_x][next_y][next_z]);
              //
              // }

              if (scalar_field[next_x][next_y][next_z]
                  > scalar_field[x][y][z]) {
                continue;
              }

              if (basin_index[next_x][next_y][next_z] >= 0) {
                labels[basin_index[next_x][next_y][next_z]]++;
              }
            }

            /// DEBUG ///
            if (debug) {
              printf("labels.size() = %d\n", static_cast<int>(labels.size()));
            }

            if (labels.empty()) {
              continue;
            }
            int max_occ = 0, label = -1;
            for (std::map<int, int>::iterator itr = labels.begin();
                itr != labels.end(); itr++) {
              if (itr->second > max_occ) {
                max_occ = itr->second;
                label = itr->first;
              }
            }
            basin_index[x][y][z] = -2 - label;

            for (int k = 0; k < kConnectivity; k++) {
              int next_x = x + kDirections[k][0];
              int next_y = y + kDirections[k][1];
              int next_z = z + kDirections[k][2];

              if (outside(next_x, next_y, next_z, nx, ny, nz)) {
                continue;
              }

              if (scalar_field[next_x][next_y][next_z]
                  > scalar_field[x][y][z]) {
                continue;
              }

              if (basin_index[next_x][next_y][next_z] != label) {
                continue;
              }

              if (dist_2_valley[x][y][z] == -1
                  || dist_2_valley[x][y][z]
                     > dist_2_valley[next_x][next_y][next_z] + 1) {
                dist_2_valley[x][y][z] =
                    dist_2_valley[next_x][next_y][next_z] + 1;
              }
            }

            flag = true;
          }
        }
      }
    }

    if (!flag) {
      break;
    }

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          if (basin_index[x][y][z] < -1) {
            basin_index[x][y][z] = -(basin_index[x][y][z] + 2);
          }
        }
      }
    }
  }
}

void shuffle_labels(int nx, int ny, int nz, int num_labels,
                    int ***basin_index, std::vector<double> *valley_height) {
  int *permutations = new int[num_labels];

  for (int i = 0; i < num_labels; i++) {
    permutations[i] = i;
  }

  std::random_shuffle(permutations, permutations + num_labels);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        basin_index[x][y][z] = permutations[basin_index[x][y][z]];
      }
    }
  }

  std::vector<double> new_valley_height(num_labels);
  for (int i = 0; i < num_labels; i++) {
    new_valley_height[permutations[i]] = (*valley_height)[i];
  }

  for (int i = 0; i < num_labels; i++) {
    (*valley_height)[i] = new_valley_height[i];
  }

  delete [] permutations;
}

void watershed_process(double ***scalar_field, int nx, int ny, int nz,
                       double *spacing,
                       int ***basin_index, int ***dist_2_valley,
                       std::vector<double> *valley_height) {
  GridPointData *data_array = new GridPointData[nx * ny * nz];
  int ***dist = create_3d_array<int>(nx, ny, nz);  // plateau distance

  calculate_plateau_distance(scalar_field, nx, ny, nz, dist);

  /// DEBUG ///
  // int mini_cnt = 0;
  // for (int x = 0; x < nx; x++) {
  //   for (int y = 0; y < ny; y++) {
  //     for (int z = 0; z < nz; z++) {
  //       if (dist[x][y][z] == -1) {
  //         printf("%d# minima: %d %d %d, height: %.20lf\n", ++mini_cnt, x, y, z, scalar_field[x][y][z]);
  //       }
  //     }
  //   }
  // }

  int count = 0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        data_array[count++] = GridPointData(
            x, y, z, dist[x][y][z], scalar_field[x][y][z]);
      }
    }
  }

  std::sort(data_array, data_array + count);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        basin_index[x][y][z] = -1;
        dist_2_valley[x][y][z] = -1;
      }
    }
  }

  int num_labels = 0; 

  int *queue_x = new int[nx * ny * nz];
  int *queue_y = new int[nx * ny * nz];
  int *queue_z = new int[nx * ny * nz];

  valley_height->clear();

  for (int i = 0; i < nx * ny * nz; i++) {
    int x = data_array[i].x_;
    int y = data_array[i].y_;
    int z = data_array[i].z_;

    if (basin_index[x][y][z] != -1) {  // a local minima previously covered
      continue;
    }

    /// DEBUG ///
    // bool debug = false;
    // if (x == 76 && y == 93 && z == 105) {
    //   debug = true;
    // }

    /// DEBUG ///
    // if (debug) {
    //   printf("dist[%d][%d][%d] = %d\n", x, y, z, dist[x][y][z]);

    //   for (int k = 0; k < kConnectivity; k++) {
    //     int next_x = x + kDirections[k][0];
    //     int next_y = y + kDirections[k][1];
    //     int next_z = z + kDirections[k][2];
    //     if (outside(next_x, next_y, next_z, nx, ny, nz)) {
    //       continue;
    //     }

    //     printf("scalar_field[%d][%d][%d] = %lf\n",
    //            next_x, next_y, next_z, scalar_field[next_x][next_y][next_z]);
    //   }
    // }

    if (dist[x][y][z] == -1) {  // local minima
      /// DEBUG ///
      // printf("label #%d: %d, %d, %d: %.20lf\n",
      //     num_labels, x, y, z, data_array[i].height_);

      int head = 0, tail = 0;
      queue_x[0] = x;
      queue_y[0] = y;
      queue_z[0] = z;
      basin_index[x][y][z] = num_labels;
      dist_2_valley[x][y][z] = 0;
      valley_height->push_back(data_array[i].height_);

      /// DEBUG ///
      // printf("%d %d %d: %lf\n", x, y, z, data_array[i].height_);

      while (head <= tail) {
        int curr_x = queue_x[head];
        int curr_y = queue_y[head];
        int curr_z = queue_z[head];
        head++;

        for (int k = 0; k < kConnectivity; k++) {
          int next_x = curr_x + kDirections[k][0];
          int next_y = curr_y + kDirections[k][1];
          int next_z = curr_z + kDirections[k][2];

          if (outside(next_x, next_y, next_z, nx, ny, nz)) {
            continue;
          }

          if (dist[next_x][next_y][next_z] != -1) {  // not in local minima
            continue;
          }

          if (basin_index[next_x][next_y][next_z] == -1) {  // not covered yet
            basin_index[next_x][next_y][next_z] = num_labels;
            tail++;
            queue_x[tail] = next_x;
            queue_y[tail] = next_y;
            queue_z[tail] = next_z;
            dist_2_valley[next_x][next_y][next_z] = 0;

            /// DEBUG ///
            // printf("%d %d %d\n", next_x, next_y, next_z);
          }
        }
      }

      num_labels++;
      continue;
    }

    int label = -1;
    int min_dist_2_valley = -1;

    for (int k = 0; k < kConnectivity; k++) {
      int next_x = x + kDirections[k][0];
      int next_y = y + kDirections[k][1];
      int next_z = z + kDirections[k][2];

      if (outside(next_x, next_y, next_z, nx, ny, nz)) {
        continue;
      }

      if (basin_index[next_x][next_y][next_z] == -1) {
        continue;
      }

      if (scalar_field[next_x][next_y][next_z] > data_array[i].height_
          || (scalar_field[next_x][next_y][next_z] == data_array[i].height_
          && dist[next_x][next_y][next_z] > data_array[i].dist_)) {
        report_error("Incorrect order found.\n");
      }

      if (scalar_field[next_x][next_y][next_z] == data_array[i].height_
          && dist[next_x][next_y][next_z] == data_array[i].dist_) {
        continue;
      }

      if (scalar_field[next_x][next_y][next_z] > scalar_field[x][y][z]) {
        report_error("Incorrect order (type 2) found.\n");
      }

      if (label == -1) {
        label = basin_index[next_x][next_y][next_z];
        min_dist_2_valley = dist_2_valley[next_x][next_y][next_z];
      } else if (label != basin_index[next_x][next_y][next_z]) {
        label = -2;
      } else if (min_dist_2_valley > dist_2_valley[next_x][next_y][next_z]) {
        min_dist_2_valley = dist_2_valley[next_x][next_y][next_z];
      }
    }

    if (label >= 0) {  // belong to a basin
      basin_index[x][y][z] = label;
      dist_2_valley[x][y][z] = min_dist_2_valley + 1;
    }
  }

  /// DEBUG ///
  // check_waterproof(basin_index, nx, ny, nz);  // not necessarily hold

  remove_seams_by_gradients(scalar_field, nx, ny, nz, spacing,
                            basin_index, dist_2_valley);
  // remove_seams(scalar_field, nx, ny, nz, basin_index, dist_2_valley);

  shuffle_labels(nx, ny, nz, num_labels, basin_index, valley_height);

  delete [] data_array;
  delete [] queue_x;
  delete [] queue_y;
  delete [] queue_z;

  delete_3d_array(dist);
}

int find_root(int node, int *father) {
  /// DEBUG ///
  if (node < 0) {
    report_error("Invalid node: %d\n", node);
  }

  if (father[node] == node) {
    return node;
  }

  return father[node] = find_root(father[node], father);
}

void merge(int node_1, int node_2, int *father) {
  father[find_root(node_1, father)] = find_root(node_2, father);
}

void output_basin(double ***scalar_field,
                  int ***basin_index,
                  int ***dist_2_valley,
                  int start_x, int start_y, int start_z,
                  int nx, int ny, int nz) {
  int label = basin_index[start_x][start_y][start_z];
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        if (basin_index[x][y][z] == label) {
          printf("%d, %d, %d: index %d, dist %d, scalar %lf\n",
                 x, y, z, basin_index[x][y][z], dist_2_valley[x][y][z],
                 scalar_field[x][y][z]);
        }
      }
    }
  } 
}

void merge_watershed(double ***scalar_field,
                     int ***basin_index,
                     int ***dist_2_valley,
                     int nx, int ny, int nz,
                     const std::vector<double> &valley_height,
                     double height_threshold,
                     double quotient_threshold,
                     int *father) {
  /// DEBUG ///
  printf("nx, ny, nz = %d, %d, %d\n", nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        /// DEBUG ///
        if (basin_index[x][y][z] < 0) {
          report_error("Invalid index (%d, %d, %d): %d\n",
                       x, y, z, basin_index[x][y][z]);
        }

        for (int k = 0; k < kConnectivity; k++) {
          int next_x = x + kDirections[k][0];
          int next_y = y + kDirections[k][1];
          int next_z = z + kDirections[k][2];
          if (outside(next_x, next_y, next_z, nx, ny, nz)) {
            continue;
          }

          /// DEBUG ///
          if (basin_index[next_x][next_y][next_z] < 0) {
            report_error("Invalid index2 (%d, %d, %d): %d\n",
                         next_x, next_y, next_z,
                         basin_index[next_x][next_y][next_z]);
          }

          if (find_root(basin_index[x][y][z], father)
              == find_root(basin_index[next_x][next_y][next_z], father)) {
            continue;
          }

          /// DEBUG ///
          // printf("x, y, z = %d, %d, %d, next_xyz = %d, %d, %d\n", x, y, z, next_x, next_y, next_z);

          // if (dist_2_valley[x][y][z] == 0) {
          //   output_basin(scalar_field, basin_index, dist_2_valley, x, y, z,
          //                nx, ny, nz);
          //   report_error("dist_2_valley[%d][%d][%d] == 0\n", x, y, z);
          // }

          // if (dist_2_valley[next_x][next_y][next_z] == 0) {
          //   output_basin(scalar_field, basin_index, dist_2_valley,
          //                next_x, next_y, next_z, nx, ny, nz);
          //   report_error("dist_2_valley[%d][%d][%d] == 0\n",
          //                next_x, next_y, next_z);
          // }

          // double quotient_1 =
          //     (scalar_field[x][y][z] - valley_height[basin_index[next_x][next_y][next_z]])
          //     / (dist_2_valley[next_x][next_y][next_z] + 1);

          // double quotient_2 =
          //     (scalar_field[next_x][next_y][next_z] - valley_height[basin_index[x][y][z]])
          //     / (dist_2_valley[x][y][z] + 1);

          // if (quotient_1 < quotient_threshold
          //     || quotient_2 < quotient_threshold) {
          //   merge(basin_index[x][y][z], basin_index[next_x][next_y][next_z],
          //         father);
          //   continue;
          // }

          double ridge_height = std::max(scalar_field[x][y][z],
                                         scalar_field[next_x][next_y][next_z]);
          double higher_valley_height = std::max(
              valley_height[basin_index[x][y][z]],
              valley_height[basin_index[next_x][next_y][next_z]]);

          if (ridge_height - higher_valley_height < height_threshold) {
            printf("ridge_height = %lf, higher_valley_height = %lf\n", ridge_height, higher_valley_height);
            printf("scalar_field[x][y][z] = %lf, valley_1 = %lf\n",
                   scalar_field[x][y][z], valley_height[basin_index[x][y][z]]);
            printf("scalar_field[next_x][next_y][next_z] = %lf, valley_2 = %lf\n",
                   scalar_field[next_x][next_y][next_z], valley_height[basin_index[next_x][next_y][next_z]]);
            merge(basin_index[x][y][z], basin_index[next_x][next_y][next_z],
                  father);
          }
        }
      }
    }
  }

  /// DEBUG ///
  printf("Finished merging.\n");
}

}

void WatershedExtractor::laplacian_smoothing(
    double ***scalar_field, int nx, int ny, int nz) {
  double ***work_field = create_3d_array<double>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        double sum = 0.0;
        int count = 0;
        for (int k = 0; k < kConnectivity; k++) {
          int next_x = x + kDirections[k][0];
          int next_y = y + kDirections[k][1];
          int next_z = z + kDirections[k][2];
          if (outside(next_x, next_y, next_z, nx, ny, nz)) {
            continue;
          }

          count++;
          sum += scalar_field[next_x][next_y][next_z];
        }

        work_field[x][y][z] = sum / count;
      }
    }
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        scalar_field[x][y][z] = work_field[x][y][z];
      }
    }
  }

  delete_3d_array(work_field);
}

void WatershedExtractor::extract_watershed(
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints **basin_index,
    vtkStructuredPoints **dist_2_valley,
    std::vector<double> *valley_height) {
  int dimensions[3];
  double spacing[3], origin[3];
  scalar_field->GetDimensions(dimensions);
  scalar_field->GetSpacing(spacing);
  scalar_field->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  *basin_index = vtkStructuredPoints::New();
  *dist_2_valley = vtkStructuredPoints::New();

  (*basin_index)->SetDimensions(dimensions);
  (*dist_2_valley)->SetDimensions(dimensions);

  (*basin_index)->SetSpacing(spacing);
  (*dist_2_valley)->SetSpacing(spacing);

  (*basin_index)->SetOrigin(origin);
  (*dist_2_valley)->SetOrigin(origin);

  int ***basin_index_array = create_3d_array<int>(nx, ny, nz);
  int ***dist_2_valley_array = create_3d_array<int>(nx, ny, nz);
  double ***scalar_field_array = create_3d_array<double>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        scalar_field_array[x][y][z] = scalar_field->GetPointData()
                                                  ->GetScalars()
                                                  ->GetTuple1(index);
      }
    }
  }

  watershed_process(scalar_field_array, nx, ny, nz, spacing,
                    basin_index_array, dist_2_valley_array,
                    valley_height);

  vtkSmartPointer<vtkIntArray> label_array =
      vtkSmartPointer<vtkIntArray>::New();

  vtkSmartPointer<vtkIntArray> dist_array =
      vtkSmartPointer<vtkIntArray>::New();

  label_array->SetNumberOfComponents(1);
  dist_array->SetNumberOfComponents(1);

  label_array->SetNumberOfTuples(nx * ny * nz);
  dist_array->SetNumberOfTuples(nx * ny * nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        label_array->SetTuple1(index, basin_index_array[x][y][z]);
        dist_array->SetTuple1(index, dist_2_valley_array[x][y][z]);
      }
    }
  }

  (*basin_index)->GetPointData()->SetScalars(label_array);
  (*dist_2_valley)->GetPointData()->SetScalars(dist_array);

  delete_3d_array(basin_index_array);
  delete_3d_array(dist_2_valley_array);
  delete_3d_array(scalar_field_array);
}

void WatershedExtractor::filter_watershed(
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints *basin_index,
    vtkStructuredPoints *dist_2_valley,
    const std::vector<double> &valley_height,
    double height_threshold,
    double quotient_threshold,
    vtkStructuredPoints **filtered_index) {
  int dimensions[3];
  scalar_field->GetDimensions(dimensions);
  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  double ***scalar_field_array = create_3d_array<double>(nx, ny, nz);
  int ***basin_index_array = create_3d_array<int>(nx, ny, nz);
  int ***dist_2_valley_array = create_3d_array<int>(nx, ny, nz);

  int num_labels = -1;

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        scalar_field_array[x][y][z] = scalar_field->GetPointData()
                                                  ->GetScalars()
                                                  ->GetTuple1(index);
        basin_index_array[x][y][z] = basin_index->GetPointData()
                                                ->GetScalars()
                                                ->GetTuple1(index);
        dist_2_valley_array[x][y][z] = dist_2_valley->GetPointData()
                                                    ->GetScalars()
                                                    ->GetTuple1(index);
        if (basin_index_array[x][y][z] > num_labels) {
          num_labels = basin_index_array[x][y][z];
        }
      }
    }
  }

  num_labels++;

  int *father = new int[num_labels];
  for (int i = 0; i < num_labels; i++) {
    father[i] = i;
  }

  /// DEBUG ///
  printf("num_labels = %d\n", num_labels);
  printf("before merge_watershed\n");

  merge_watershed(scalar_field_array, basin_index_array, dist_2_valley_array,
                  nx, ny, nz, valley_height, height_threshold, quotient_threshold, father);

  /// DEBUG ///
  printf("after merge_watershed\n");

  int *color = new int[num_labels];

  int new_num_labels = 0;
  for (int i = 0; i < num_labels; i++) {
    if (father[i] != i) {
      color[i] = -1;
    } else {
      color[i] = new_num_labels++;
    }
  }

  printf("new_num_labels = %d\n", new_num_labels);

  *filtered_index = vtkStructuredPoints::New();
  (*filtered_index)->CopyStructure(scalar_field);
  vtkSmartPointer<vtkDoubleArray> labels = vtkSmartPointer<vtkDoubleArray>::New();
  labels->SetNumberOfComponents(1);
  labels->SetNumberOfTuples(nx * ny * nz);
  for (int i = 0; i < nx * ny * nz; i++) {
    int curr_label = basin_index->GetPointData()->GetScalars()->GetTuple1(i);
    int curr_color = color[find_root(curr_label, father)];
    labels->SetTuple1(i, static_cast<double>(curr_color));
  }

  (*filtered_index)->GetPointData()->SetScalars(labels);

  delete [] father;
  delete [] color;

  delete_3d_array(scalar_field_array);
  delete_3d_array(basin_index_array);
  delete_3d_array(dist_2_valley_array);
}

vtkStructuredPoints *WatershedExtractor::image_data_2_structured_points(
    vtkImageData *image) {
  int dimensions[3];
  double origin[3];
  double spacing[3];

  image->GetDimensions(dimensions);
  image->GetOrigin(origin);
  image->GetSpacing(spacing);

  vtkStructuredPoints *structured = vtkStructuredPoints::New();

  structured->SetDimensions(dimensions);
  structured->SetOrigin(origin);
  structured->SetSpacing(spacing);

  structured->GetPointData()->SetScalars(image->GetPointData()->GetScalars());

  return structured;
}

