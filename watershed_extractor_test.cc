// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_extractor.h"

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

const char *kScalarFile = "data/sphere_ftle.vtk";
const char *kBasinFile = "basin_index.vtk";
const char *kDistFile = "dist_2_valley.vtk";

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

int main() {
  extract_watershed_test();

  return 0;
}
