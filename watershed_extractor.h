// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef WATERSHED_EXTRACTOR_H_
#define WATERSHED_EXTRACTOR_H_

class vtkStructuredPoints;

class WatershedExtractor {
 public:
  // basin_index stores the index of regions, starting from 0.
  // dist_2_valley stores the steps from the valley (local minima).
  void extract_watershed(vtkStructuredPoints *scalar_field,
                         vtkStructuredPoints **basin_index,
                         vtkStructuredPoints **dist_2_valley);
};

#endif  // WATERSHED_EXTRACTOR_H_