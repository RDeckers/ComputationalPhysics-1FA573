#include <stdio.h>
#include <math.h>

double test_function(double x){
  return exp(x);
}

double trapezoidal(double (*func)(double), double x0, double x1, unsigned N){
  double h = (x1-x0)/N;
  double x = x0+h;
  double sum = 0;
  for(unsigned n = 0; n < N; n+=2){
    sum += (h/2)*(func(x-h)+2*func(x)+func(x+h));
    x += 2*h;
  }
  return sum;
}

double simpsons(double (*func)(double), double x0, double x1, unsigned N){
  double h = (x1-x0)/N;
  double x = x0+h;
  double sum = 0;
  for(unsigned n = 0; n < N; n+=2){
    sum += (h/3)*(func(x+h)+4*func(x)+func(x-h));
    x += 2*h;
  }
  return sum;
}

int main(){
  double x1 = 1;
  double x0 = 0;
  double exact = M_E - 1;
  printf("#h Trapezoidal Simpsons");
  for(unsigned n = 2; n <= (1 << 23); n*=2 ){
    double h = (x1-x0)/n;
    printf("\n%e %e %e",
     h,
     fabs(exact-trapezoidal(&test_function, x0, x1, n)),
     fabs(exact-simpsons(&test_function, x0, x1, n))
   );
  }
  return 0;
}
