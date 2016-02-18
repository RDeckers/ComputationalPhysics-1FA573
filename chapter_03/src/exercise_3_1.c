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

double drv(double x, double y){
  return -2*M_PI*sin(2*M_PI*x);
}

double k_sq(double x){
  return 4*M_PI*M_PI;
}

double S(double x){
  return 0;
}

double exact(double x){
  return cos(2*M_PI*x);
}

int main(int argc, char **argv){
  double y0 = 1;
  double x0 = 0;
  double x_end = 1.25;
  double exact_end = exact(x_end);
  printf("#h Numerov rk_4");
  for(unsigned n = 8; n <= (1 << 15); n*=2){
    double h = subdivide(x0, x_end, n);
    double y1 = rk_4(&drv, x0, y0, x0+h,1);
    printf(
      "\n%e %e %e",
      h,
      abs_diff(Numerov(&S, &k_sq, x0, y0, y1, x_end, n), exact_end),
      abs_diff(rk_4(&drv, x0, y0, x_end, n), exact_end)
    );
  }
  return 0;
}
