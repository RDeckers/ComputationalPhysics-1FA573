#include <math.h>
#include <stdio.h>

double test_function(double x){
  return x*x-5;
}

double test_function_deriv(double x){
  return 2*x;
}

void search(double (*func)(double), double *X, double x0, double dx0, unsigned N){
  double dx = dx0;
  double x = x0;
  X[0] = x0;
  double direction = copysign(1.0, func(x));
  for(unsigned u = 1; u < N; u++){
    x += dx;
    if(direction*func(x) < 0){
      dx *= -0.5;
      direction *= -1;
    }
    X[u] = x;
  }
  return;
}

void Newton_Rhapson(double (*func)(double), double (*deriv)(double), double *X, double x0, unsigned N){
  double x = x0;
  X[0] = x0;
  for(unsigned u = 1; u < N; u++){
    x -= func(x)/deriv(x);
    X[u] = x;
  }
}

void secant(double (*func)(double), double *X, double x0, unsigned N){
  double x = x0;
  X[0] = x;
  x -= func(x);
  X[1] = x;
  for(unsigned u = 2; u < N; u++){
    double x_next = x - func(x)*(x-X[u-2])/(func(x)-func(X[u-2]));
    if(isnan(x_next)){
      for(; u < N; u++)
        X[u] = x;
      return;
    }
    x = x_next;
    X[u] = x;
  }
}


int main(int argc, char **argv){
  double x0 = 1;
  double exact = sqrt(5);
  const unsigned N = 16;
  double X_search[N], X_NR[N], X_secant[N];
  search(&test_function, X_search, x0, 1.0, N);
  Newton_Rhapson(&test_function, &test_function_deriv, X_NR, x0, N);
  secant(&test_function, X_secant, x0, N);
  printf("#iteration search Newton_Rhapson secant");
  for(unsigned n = 0;n < N; n++){
    printf("\n%u %e %e %e",
      n,
      fabs(exact-X_search[n]),
      fabs(exact-X_NR[n]),
      fabs(exact-X_secant[n])
    );
  }
  return 0;
}
