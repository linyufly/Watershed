// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_extractor.h"

#include "util.h"

#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkIndent.h>

#include <cstdio>
#include <cstdlib>

#include <iostream>

// const char *kScalarFile = "data/one_sphere_ftle.vtk";
// const char *kScalarFile = "data/sphere_ftle.vtk";
const char *kScalarFile = "data/gyre_half.vtk";
// const char *kScalarFile = "smoothed_scalar.vtk";
const char *kBasinFile = "basin_index.vtk";
const char *kDistFile = "dist_2_valley.vtk";
const char *kScalarToSmoothFile = "data/gyre_half.vtk";
const char *kSmoothedScalarFile = "smoothed_scalar.vtk";

const int kNumberOfSmoothing = 80;

void extract_watershed_test() {
  printf("extract_watershed_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kScalarFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> scalar_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  scalar_field->ShallowCopy(reader->GetOutput());

  WatershedExtractor extractor;
  vtkStructuredPoints *basin_index = NULL, *dist_2_valley = NULL;
  extractor.extract_watershed(scalar_field, &basin_index, &dist_2_valley);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
    vtkSmartPointer<vtkStructuredPointsWriter>::New();

  writer->SetFileName(kBasinFile);
  writer->SetInputData(basin_index);
  writer->Write();

  writer->SetFileName(kDistFile);
  writer->SetInputData(dist_2_valley);
  writer->Write();

  if (basin_index) {
    basin_index->Delete();
  }

  if (dist_2_valley) {
    dist_2_valley->Delete();
  }

  printf("} extract_watershed_test\n\n");
}

void laplacian_smoothing_test() {
  printf("laplacian_smoothing_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kScalarToSmoothFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> scalar_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  scalar_field->ShallowCopy(reader->GetOutput());

  int dimensions[3];
  scalar_field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  double ***scalar_array = create_3d_array<double>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        scalar_array[x][y][z] = scalar_field->GetPointData()
                                            ->GetScalars()
                                            ->GetTuple1(index);
      }
    }
  }

  WatershedExtractor extractor;
  for (int i = 0; i < kNumberOfSmoothing; i++) {
    extractor.laplacian_smoothing(scalar_array, nx, ny, nz);
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        scalar_field->GetPointData()
                    ->GetScalars()
                    ->SetTuple1(index, scalar_array[x][y][z]);
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();

  writer->SetFileName(kSmoothedScalarFile);
  writer->SetInputData(scalar_field);
  writer->Write();

  delete_3d_array(scalar_array);

  printf("} laplacian_smoothing_test\n\n");
}

int main() {
  extract_watershed_test();
  laplacian_smoothing_test();

  return 0;
}
