#include <stdio.h>
#include <math.h>
#include <common.h>
#include <integrators.h>

double test_function_drv(double x, double y){
  return -x*y;
}

double test_function_g(double x){
  return -x;
}

double exact(double x){
  return exp(-x*x*0.5);
}

int main(int argc, char **argv){
  double y0 = 1;
  double x0 = 0;
  double exact_1 = exact(1);
  printf("#h implicit rk2 rk3 rk4");
  for(unsigned n = 2; n <= (1 << 15); n*=2){
    double h = subdivide(x0, 1, n);
    printf(
      "\n%e %e %e %e %e",
      h,
      rel_abs_diff(implicit_method(&test_function_g, x0, y0, 1, n),exact_1),
      rel_abs_diff(rk_2(&test_function_drv, x0, y0, 1, n),exact_1),
      rel_abs_diff(rk_3(&test_function_drv, x0, y0, 1, n),exact_1),
      rel_abs_diff(rk_4(&test_function_drv, x0, y0, 1, n),exact_1)
    );
  }
  return 0;
}
