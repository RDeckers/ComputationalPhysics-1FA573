#include <stdio.h>
#include <math.h>

double test_function(double x){
  return sin(x);
}

double ddx_sym_3(double (*func)(double), double x0, double h){
  return (func(x0+h)-func(x0-h))/(2*h);
}

double ddx_forward_2(double (*func)(double), double x0, double h){
  return (func(x0+h)-func(x0))/h;
}

double ddx_backward_2(double (*func)(double), double x0, double h){
  return (func(x0)-func(x0-h))/h;
}

int main(){
  double x0 = 1;
  double exact = cos(x0);
  printf("#h sym_3 forward_2 backward_2\n");
  for(double h = 0.5; h > 1e-9; h /= 2){
    printf("%e %e %e %e\n",
     h,
     fabs(exact-ddx_sym_3(&test_function, x0, h)),
     fabs(exact-ddx_forward_2(&test_function, x0, h)),
     fabs(exact-ddx_backward_2(&test_function, x0, h))
   );
  }
  return 0;
}
