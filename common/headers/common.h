#include <math.h>

inline double subdivide(double x0, double x1, unsigned N){
  return(x1-x0)/N;
}

inline double abs_diff(double y0, double y1){
  return fabs(y0-y1);
}

inline double rel_abs_diff(double y0, double y1){
  return fabs(y0-y1)/y1;
}
