// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_extractor.h"

#include "util.h"

#include <set>
#include <algorithm>

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
  int x_, y_, z_;
  double height_;

  GridPointData(int x, int y, int z, height)
      : x_(x), y_(y), z_(z), height_(height);

  bool operator < (const GridPointData &other) {
    return height_ < other.height_;
  }
};

bool outside(int x, int y, int z, int nx, int ny, int nz) {
  return x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz;
}

void watershed_process(double ***scalar_field, int nx, int ny, int nz,
                       int ***basin_index, int ***dist_2_valley) {
  GridPointData *data_array = new GridPointData[nx * ny * nz];
  int ***dist = create_3d_array<int>(nx, ny, nz);  // plateau distance
  int *queue_x = new int[nx * ny * nz * 2];  // double size for fictitous
  int *queue_y = new int[nx * ny * nz * 2];
  int *queue_z = new int[nx * ny * nz * 2];

  int count = 0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        data_array[count++] = GridPointData(x, y, z, scalar_field[x][y][z]);
      }
    }
  }

  std::sort(data_array, data_array + nx * ny * nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        basin_index[x][y][z] = -1;
        dist_2_valley[x][y][z] = -1;
        dist[x][y][z] = 0;
      }
    }
  }

  int num_labels = 0, head = 0, tail = -1;

  for (int curr = 0; curr < nx * ny * nz; ) {
    double height = data_array[curr].height_;

    int succ = curr + 1;
    while (succ < nx * ny * nz && data_array[succ].height_ == height) {
      succ++;
    }

    for (int i = curr; i < succ; i++) {
      int x = data_array[i].x_;
      int y = data_array[i].y_;
      int z = data_array[i].z_;

      basin_index[x][y][z] = -2;

      bool flag = false;
      for (int k = 0; k < kConnectivity; k++) {
        int next_x = x + kDirections[k][0];
        int next_y = y + kDirections[k][1];
        int next_z = z + kDirections[k][2];

        if (outside(next_x, next_y, next_z, nx, ny, nz)) {
          continue;
        }

        if (basin_index[next_x][next_y][next_z] >= 0) {  // basin or watershed
          flag = true;
          break;
        }
      }

      if (flag) {
        dist[x][y][z] = 1;
        tail++;
        queue_x[tail] = x;
        queue_y[tail] = y;
        queue_z[tail] = z;
      }
    }

    int curr_dist = 1;
    tail++;
    queue_x[tail] = -1;

    while (true) {
      int x = queue_x[head];
      int y = queue_y[head];
      int z = queue_z[head];
      head++;

      if (x == -1) {
        if (head > tail) {
          break;
        } else {
          tail++;
          queue_x[tail] = -1;
          curr_dist++;
        }
        x = queue_x[head];
        y = queue_y[head];
        z = queue_z[head];
        head++;
      }

      for (int k = 0; k < kConnectivity; k++) {
        int next_x = x + kDirections[k][0];
        int next_y = y + kDirections[k][1];
        int next_z = z + kDirections[k][2];

        if (outside(next_x, next_y, next_z, nx, ny, nz)) {
          continue;
        }

        if (dist[next_x][next_y][next_z] < curr_dist
            && basin_index[next_x][next_y][next_z] >= 0) {
          if (dist[next_x][next_y][next_z] != curr_dist - 1) {
            report_error("Found unexpected dist\n");
          }

          if (basin_index[next_x][next_y][next_z] > 0) {  // non-ridge
            if (basin_index[x][y][z] <= 0) {  // ridge or not assigned
              basin_index[x][y][z] = basin_index[next_x][next_y][next_z];
            }
          }
        }
      }
    }
  }

  

  delete [] data_array;
  delete [] queue_x;
  delete [] queue_y;
  delete [] queue_z;

  delete_3d_array(dist);
}

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

  delete_3d_array(basin_index_array);
  delete_3d_array(dist_2_valley_array);
  delete_3d_array(scalar_field_array);
}
