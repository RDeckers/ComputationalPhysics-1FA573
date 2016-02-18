#include <common.h>
#include <integrators.h>
#include <stdio.h>
#include <string.h>

/*NOTE:
  tables don't match exactly but that seems to be mostly attributable to different values
  for the exact answer.
*/

double exact(double r){
  return 1.0-(0.5*r+1)*exp(-r);
}

double S(double r){
  return -r*exp(-r)/2;
}

double r_sq(double r){
  return 0;
}

void linear_correction(double *y, unsigned x0, unsigned x1, double h, unsigned N){
  double y0 = y[x0], y1 = y[x1];
  double slope = (y1-y0)/((x1-x0)*h);
  for(unsigned n = 0; n < N; n++){
    y[n] -= slope*n*h;
  }
}

int main(int argc, char **argv){
  unsigned n_steps = 200;
  double h = subdivide(0,20, n_steps);
  double y0 = 0;
  double y1 = exact(h);
  printf("#r exact_y1 perturbed_y1");
  double y[201] = {y0, y1}, y_perturbed[201] = {y0, 0.95*y1}, y_linear[201];
  Numerov_into_array(&S, &r_sq, 0, 20, y, 200);
  Numerov_into_array(&S, &r_sq, 0, 20, y_perturbed, 200);
  memcpy(y_linear, y, 201*sizeof(double));
  linear_correction(y_linear, 190, 200, h, 201);
  for(int r_end = 1; r_end <= 20; r_end += 1){
    double y_exact = exact(r_end);
    printf(
      "\n%d %e %e %e",
      r_end,
      abs_diff(y[10*r_end], y_exact),
      abs_diff(y_perturbed[10*r_end], y_exact),
      abs_diff(y_linear[10*r_end], y_exact)
    );
  }
}
