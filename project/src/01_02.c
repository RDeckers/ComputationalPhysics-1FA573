#include <stdio.h>
#include <math.h>

double exact_integrand(double r, double a, double b){
  return -(1.0/b)*asin(b/(r*sqrt(1-a)));
}

double exact_integral(double a, double b, double r_max){
  double result = 0;
  double r_min = b/sqrt(1-a);
  if(b < r_max){//only effected if we approach close enough
    result += exact_integrand(r_max, 0, b) - exact_integrand(b, 0, b);
    if((r_min < r_max) && (1-a > 0)){ //
      result -= exact_integrand(r_max, a, b) +1.0/b*asin(1);
    }
  }
  result *= 2*b;
  return result;
}

int main(int argc, char**argv){
  double r_max = 1.0;
  for(double ve = -1.5; ve <= 1.5; ve += 0.125){
    printf("\"V/E = %+1.2f\"\n", ve);
    for(double b = 0.001; b <= r_max; b+=0.01){
      printf("%f %f\n", b, exact_integral(ve, b, r_max));
    }
    puts("\n");
  }
  return 0;
}
