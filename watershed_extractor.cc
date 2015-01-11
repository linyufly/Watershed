// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_extractor.h"

#include "util.h"

#include <algorithm>
#include <map>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

const int kConnectivity = 6;

const int kDirections[6][3] = {
    {0, 0, 1}, {0, 0, -1},
    {0, 1, 0}, {0, -1, 0},
    {1, 0, 0}, {-1, 0, 0}};

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

void remove_seams(double ***scalar_field, int nx, int ny, int nz,
                  int ***basin_index, int ***dist_2_valley) {
  while (true) {
    bool flag = false;

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          if (basin_index[x][y][z] == -1) {
            std::map<int, int> labels;
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

              if (basin_index[next_x][next_y][next_z] >= 0) {
                labels[basin_index[next_x][next_y][next_z]]++;
              }
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
                    int ***basin_index) {
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

  delete [] permutations;
}

void watershed_process(double ***scalar_field, int nx, int ny, int nz,
                       int ***basin_index, int ***dist_2_valley) {
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

  for (int i = 0; i < nx * ny * nz; i++) {
    int x = data_array[i].x_;
    int y = data_array[i].y_;
    int z = data_array[i].z_;

    if (basin_index[x][y][z] != -1) {  // a local minima previously covered
      continue;
    }

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

  remove_seams(scalar_field, nx, ny, nz, basin_index, dist_2_valley);

  shuffle_labels(nx, ny, nz, num_labels, basin_index);

  delete [] data_array;
  delete [] queue_x;
  delete [] queue_y;
  delete [] queue_z;

  delete_3d_array(dist);
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
    vtkStructuredPoints **dist_2_valley) {
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

  watershed_process(scalar_field_array, nx, ny, nz,
                    basin_index_array, dist_2_valley_array);

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
        if (scalar_field_array[x][y][z] > num_labels) {
          num_labels = scalar_field_array[x][y][z];
        }
      }
    }
  }

  num_labels++;

   

  delete_3d_array(scalar_field_array);
  delete_3d_array(basin_index_array);
  delete_3d_array(dist_2_valley_array);
}

