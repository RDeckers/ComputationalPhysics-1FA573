#include <stdio.h>
#include <math.h>
#include <common.h>
#include <integrators.h>


double dydt(double t, double p, double y){
    return p;
}

double dpdt(double t, double p, double y){
    return -4*M_PI*M_PI*y;
}

int main(int argc, char** argv){
  unsigned iteration_count = 1;
  printf("#end location after %u iterations\n#dt p y", iteration_count);
  for(unsigned n = 16; n <= (1<<18); n *= 2){
    double p = 1, y = 0;
    rk_4_2d(
      &dpdt, &dydt,
      &p, &y,
      0, iteration_count,
      n
    );
    double dt = subdivide(0, iteration_count, n);
    printf("\n%e %e %e", dt, abs_diff(p,1), abs_diff(y,0));
  }
  return 0;
}
