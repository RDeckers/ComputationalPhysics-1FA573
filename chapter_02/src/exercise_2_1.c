#include <stdio.h>
#include <common.h>
#include <integrators.h>


double test_function_drv(double x, double y){
  return -x*y;
}

double test_function_drv2_x(double x, double y){
  return -y;
}

double test_function_drv2_y(double x, double y){
  return -x;
}

double test_function_g(double x){
  return -x;
}

int main(int argc, char** argv){
  double y0 = 1;
  double x0 = 0;
  double exact = 1;
  printf("#h euler_1 taylor_1 implicit_1 euler_3 taylor_3 implicit_3");
  for(unsigned n = 2; n <= (1 << 23); n*=2){
    double h = subdivide(x0, 1, n);
    printf(
      "\n%e %e %e %e %e %e %e",
      h,
      rel_abs_diff(
        Eulers_method(&test_function_drv, 1, Eulers_method(&test_function_drv, x0, y0, 1, n), 0, n),
        exact
      ),
      rel_abs_diff(
        Taylor_series(
          &test_function_drv,
          &test_function_drv2_x,
          &test_function_drv2_y,
          1,
          Taylor_series(&test_function_drv, &test_function_drv2_x, test_function_drv2_y, x0, y0, 1, n),
          0,
          n
        ),
        exact
      ),
      rel_abs_diff(
        implicit_method(&test_function_g, 1, implicit_method(&test_function_g, x0, y0, 1, n), 0, n),
        exact
      ),
      rel_abs_diff(
        Eulers_method(&test_function_drv, 3, Eulers_method(&test_function_drv, x0, y0, 3, 3*n), 0, 3*n),
        exact
      ),
      rel_abs_diff(
        Taylor_series(
          &test_function_drv,
          &test_function_drv2_x,
          &test_function_drv2_y,
          3,
          Taylor_series(&test_function_drv, &test_function_drv2_x, test_function_drv2_y, x0, y0, 3, 3*n),
          0,
          3*n
        ),
        exact
      ),
      rel_abs_diff(
        implicit_method(&test_function_g, 3, implicit_method(&test_function_g, x0, y0, 3, 3*n), 0, 3*n),
        exact
      )
    );
  }

  return 0;
}
