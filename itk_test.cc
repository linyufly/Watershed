// Author: Mingcheng Chen (linyufly@gmail.com)

#include <cstdio>

#include <algorithm>
#include <vector>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkSize.h>
#include <itkWatershedImageFilter.h>

typedef itk::Image<double, 3> ImageType;
typedef itk::WatershedImageFilter<ImageType> WatershedType;

// const char *kScalarFile = "smoothed_scalar.vtk";
// const char *kScalarFile = "/home/linyufly/Data/P96_bFTLE.vtk";
// const char *kScalarFile = "/home/linyufly/Data/abcflow_200.vtk";
// const char *kScalarFile = "/home/linyufly/GitHub/Playground/masked.vtk";
// const char *kScalarFile = "/home/linyufly/GitHub/Playground/cropped_convective_half.vtk";
// const char *kScalarFile = "/home/linyufly/GitHub/Playground/symmetrical_convective_half.vtk";
const char *kScalarFile = "/home/linyufly/Data/symmetrical_convective_half.vtk";

const char *kITKWatershedImageFilterResultFile =
    "itk_watershed_image_filter_result.vtk";

// const char *kScalarFile = "test_scalar.vtk";
// const char *kITKWatershedImageFilterResultFile = "test_label.vtk";

void watershed_image_filter_test() {
  printf("watershed_image_filter_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kScalarFile);
  reader->Update();

  vtkStructuredPoints *scalar_field = reader->GetOutput();

  int dimensions[3];
  double spacing[3], origin[3];

  scalar_field->GetDimensions(dimensions);
  scalar_field->GetSpacing(spacing);
  scalar_field->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  ImageType::Pointer scalar_image = ImageType::New();

  ImageType::SizeType image_size;
  for (int dimension = 0; dimension < 3; dimension++) {
    image_size[dimension] = dimensions[dimension];
  }

  scalar_image->SetRegions(image_size);
  scalar_image->Allocate();

  ImageType::IndexType image_index;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        image_index[0] = x;
        image_index[1] = y;
        image_index[2] = z;

        double value = scalar_field->GetPointData()
                                     ->GetScalars()
                                     ->GetTuple1((z * ny + y) * nx + x);

        scalar_image->SetPixel(
            image_index, value);
      }
    }
  }

  printf("Finished creating scalar_image.\n");

  WatershedType::Pointer watershed = WatershedType::New();
  watershed->SetInput(scalar_image);
  watershed->SetThreshold(0.0);
  watershed->SetLevel(0.1);  // 0.2
  watershed->Update();

  int last_region = 0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        image_index[0] = x;
        image_index[1] = y;
        image_index[2] = z;

        int region = static_cast<int>(
            watershed->GetOutput()->GetPixel(image_index));
        if (region > last_region) {
          last_region = region;
        }
      }
    }
  }

  std::vector<int> renumber(last_region + 1);
  for (int code = 0; code <= last_region; code++) {
    renumber[code] = code;
  }

  random_shuffle(renumber.begin() + 1, renumber.end());

  vtkSmartPointer<vtkStructuredPoints> watershed_result =
      vtkSmartPointer<vtkStructuredPoints>::New();

  watershed_result->CopyStructure(scalar_field);

  vtkSmartPointer<vtkDoubleArray> label_array =
      vtkSmartPointer<vtkDoubleArray>::New();

  label_array->SetNumberOfComponents(1);
  label_array->SetNumberOfTuples(nx * ny * nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        image_index[0] = x;
        image_index[1] = y;
        image_index[2] = z;

        label_array->SetTuple1(
            (z * ny + y) * nx + x,
            renumber[static_cast<int>(
                watershed->GetOutput()->GetPixel(image_index))]);
      }
    }
  }

  watershed_result->GetPointData()->SetScalars(label_array);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetInputData(watershed_result);
  writer->SetFileName(kITKWatershedImageFilterResultFile);
  writer->Write();

  printf("} watershed_image_filter_test\n");
  printf("\n");
}

int main() {
  watershed_image_filter_test();

  return 0;
}

