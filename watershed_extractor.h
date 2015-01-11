// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef WATERSHED_EXTRACTOR_H_
#define WATERSHED_EXTRACTOR_H_

#include <vector>

class vtkStructuredPoints;

class WatershedExtractor {
 public:
  // This method conducts one iteration of Laplacian smoothing.
  void laplacian_smoothing(double ***scalar_field, int nx, int ny, int nz);

  // basin_index stores the index of regions, starting from 0.
  // dist_2_valley stores the steps from the valley (local minima).
  void extract_watershed(vtkStructuredPoints *scalar_field,
                         vtkStructuredPoints **basin_index,
                         vtkStructuredPoints **dist_2_valley,
                         std::vector<double> *valley_height);

  // quotient_threshold is for the quotient of height difference over distance
  // to valley.
  void filter_watershed(vtkStructuredPoints *scalar_field,
                        vtkStructuredPoints *basin_index,
                        vtkStructuredPoints *dist_2_valley,
                        const std::vector<double> &valley_height,
                        double quotient_threshold,
                        vtkStructuredPoints **filtered_index);
};

#endif  // WATERSHED_EXTRACTOR_H_
