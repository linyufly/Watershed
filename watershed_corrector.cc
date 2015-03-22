// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_corrector.h"

#include "util.h"

#include <set>

#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

const int kConnectivity = 6;

const int kDirections[26][3] = {
    {0, 0, 1}, {0, 0, -1},
    {0, 1, 0}, {0, -1, 0},
    {1, 0, 0}, {-1, 0, 0},

    {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1},
    {1, 0, 1}, {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
    {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},

    {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
    {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};

int get_index(int *coord, int *dimensions) {
  return (coord[2] * dimensions[1] + coord[1]) *
         dimensions[0] + coord[0];
}

bool outside(int *coord, int *dimensions) {
  for (int d = 0; d < 3; d++) {
    if (coord[d] < 0 || coord[d] >= dimensions[d]) {
      return true;
    }
  }

  return false;
}

bool update(bool ***used, int x, int y, int z, int nx, int ny, int nz,
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints *label_field) {
  used[x][y][z] = true;

  int dimensions[] = {nx, ny, nz};
  int curr_coord[] = {x, y, z};
  int curr_index = get_index(curr_coord, dimensions);

  double curr_scalar =
      scalar_field->GetPointData()->GetScalars()->GetTuple1(curr_index);
  double curr_label =
      label_field->GetPointData()->GetScalars()->GetTuple1(curr_index);

  double min_scalar = curr_scalar;
  std::set<double> min_labels;

  for (int d = 0; d < kConnectivity; d++) {
    int next_coord[3];
    for (int c = 0; c < 3; c++) {
      next_coord[c] = curr_coord[c] + kDirections[d][c];
    }

    if (outside(next_coord, dimensions)) {
      continue;
    }

    int next_index = get_index(next_coord, dimensions);
    double next_scalar =
        scalar_field->GetPointData()->GetScalars()->GetTuple1(next_index);
    double next_label =
        label_field->GetPointData()->GetScalars()->GetTuple1(next_index);

    if (next_scalar < min_scalar) {
      min_scalar = next_scalar;
      min_labels.clear();
    }

    if (next_scalar == min_scalar) {
      min_labels.insert(next_label);
    }
  }

  if (min_labels.size() == 0) {
    return false;
  }

  if (min_labels.find(curr_label) == min_labels.end()) {
    double new_label = *min_labels.begin();
    label_field->GetPointData()->GetScalars()->SetTuple1(curr_index, new_label);

    for (int d = 0; d < kConnectivity; d++) {
      int next_coord[3];
      for (int c = 0; c < 3; c++) {
        next_coord[c] = curr_coord[c] + kDirections[d][c];
      }

      if (outside(next_coord, dimensions)) {
        continue;
      }

      if (used[next_coord[0]][next_coord[1]][next_coord[2]]) {
        continue;
      }

      int next_index = get_index(next_coord, dimensions);
      double next_scalar = scalar_field->GetPointData()->GetScalars()->GetTuple1(next_index);
      if (next_scalar <= curr_label) {
        continue;
      }

      update(used, next_coord[0], next_coord[1], next_coord[2], nx, ny, nz, scalar_field, label_field);
    }

    return true;
  }

  return false;
}

}

vtkStructuredPoints *WatershedCorrector::correct(
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints *old_label_field) {
  vtkStructuredPoints *label_field = vtkStructuredPoints::New();
  label_field->DeepCopy(old_label_field);

  int dimensions[3];
  scalar_field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  bool ***used = create_3d_array<bool>(nx, ny, nz);

  /// DEBUG ///
  int num_updates = 0;

  while (true) {
    bool updated = false;

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          used[x][y][z] = false;
        }
      }
    }

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
          if (update(used, x, y, z, nx, ny, nz, scalar_field, label_field)) {
            updated = true;
          }
        }
      }
    }

    if (!updated) {
      break;
    }

    /// DEBUG ///
    num_updates++;
    printf("num_updates = %d\n", num_updates);
  }

  delete_3d_array(used);

  return label_field;
}
