// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef WATERSHED_CORRECTOR_H_
#define WATERSHED_CORRECTOR_H_

class vtkStructuredPoints;

class WatershedCorrector {
 public:
  static vtkStructuredPoints *correct(
      vtkStructuredPoints *scalar_field,
      vtkStructuredPoints *label_field);
};

#endif  // WATERSHED_CORRECTOR_H_
