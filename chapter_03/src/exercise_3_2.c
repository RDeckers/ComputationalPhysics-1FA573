#include <common.h>
#include <integrators.h>
#include <stdio.h>

/*NOTE
  Numerov performs worse in terms of accuracy than rk_4 for a given h
  However, it is significantlly faster as it only involves 2 function evaluations per
  iteration as opposed to 4. In this specific case the evaluations are also much cheaper.
  Taking the average speed-up to be a factor 2 makes Numerov better in terms of accuracy per
  computation time.
*/

double k_sq(double x){
  return 0;
}

double S(double r){
  return -0.5*r*exp(-r);
}

double exact(double r){
  return 1-0.5*(r+2)*exp(-r);
}

int main(int argc, char **argv){
  double x0 = 1e3;
  double y0 = 1+x0;
  double x_end = 0;
  double exact_end = 0;
  printf("#h Numerov");
  for(unsigned n = 8; n <= (1 << 20); n*=2){
    double h = subdivide(x0, x_end, n);
    double y1 = y0+h;
    printf(
      "\n%e %e",
      fabs(h),
      Numerov(&S, &k_sq, x0, y0, y1, x_end, n)
    );
  }
  return 0;
}
