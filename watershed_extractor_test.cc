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
// const char *kScalarFile = "data/gyre_half.vtk";
const char *kScalarFile = "smoothed_scalar.vtk";
const char *kBasinFile = "basin_index.vtk";
const char *kDistFile = "dist_2_valley.vtk";
// const char *kScalarToSmoothFile = "data/gyre_half.vtk";
const char *kScalarToSmoothFile = "data/output_200.vtk";
// const char *kScalarToSmoothFile = "data/output.vtk";
const char *kSmoothedScalarFile = "smoothed_scalar.vtk";
const char *kFilteredBasinFile = "filtered_basin.vtk";

const int kNumberOfSmoothing = 10;  // 1 for gyre_half.vtk
                                   // 10 for output_200.vtk

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
  std::vector<double> valley_height;
  extractor.extract_watershed(
      scalar_field, &basin_index, &dist_2_valley, &valley_height);

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

void filter_watershed_test() {
  printf("filter_watershed_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kScalarFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> scalar_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  scalar_field->ShallowCopy(reader->GetOutput());

  WatershedExtractor extractor;
  vtkStructuredPoints *basin_index = NULL, *dist_2_valley = NULL;
  std::vector<double> valley_height;
  extractor.extract_watershed(
      scalar_field, &basin_index, &dist_2_valley, &valley_height);

  vtkStructuredPoints *filtered_index = NULL;
  extractor.filter_watershed(scalar_field, basin_index, dist_2_valley,
                             valley_height, 0.01, 0.000015, &filtered_index);
  // gyre_half
  // good height_threshold for 10 / 20 : 0.001, 0.002, 0.003, 0.004, 0.005
  // good height_threshold for 1: 0.0005, 0.001

  // output_200
  // good height_threshold for 10: 0.006, 0.005, 0.01?

  // output
  // good height_threshold for 10: 

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
    vtkSmartPointer<vtkStructuredPointsWriter>::New();

  writer->SetFileName(kFilteredBasinFile);
  writer->SetInputData(filtered_index);
  writer->Write();

  if (basin_index) {
    basin_index->Delete();
  }

  if (dist_2_valley) {
    dist_2_valley->Delete();
  }

  if (filtered_index) {
    filtered_index->Delete();
  }

  printf("} filter_watershed_test\n\n");
}

int main() {
  // extract_watershed_test();
  // laplacian_smoothing_test();
  filter_watershed_test();

  return 0;
}
