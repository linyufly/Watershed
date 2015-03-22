// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_corrector.h"

#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

namespace {

const char *kScalarFieldFile = "/home/linyufly/Data/symmetrical_convective_half.vtk";
const char *kLabelFieldFile = "label_field.vtk";
const char *kCorrectOutputFile = "correct.vtk";

}

void correct_test() {
  printf("correct_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> scalar_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  scalar_field->DeepCopy(reader->GetOutput());

  reader->SetFileName(kLabelFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> label_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  label_field->DeepCopy(reader->GetOutput());

  vtkSmartPointer<vtkStructuredPoints> corrected =
      vtkSmartPointer<vtkStructuredPoints>(WatershedCorrector::correct(scalar_field, label_field));

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetInputData(corrected);
  writer->SetFileName(kCorrectOutputFile);
  writer->Write();

  printf("} correct_test\n");
}

int main() {
  correct_test();

  return 0;
}

