#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <poly.h>
#include <gk.h>
#include <utilities/logging.h>

double second_deriv(double gamma, double c, double r){
  return gamma*(156*pow(c,12)*pow(r, -14)-42*pow(c,6)*pow(r,-8))+6*pow(r,-4);
}

int main(int argc, char** argv){
  double c = 1;
  double gamma = 10;
  double step_size = 0.1;
  while(step_size >= 1e-13){
    gamma -= step_size;
    poly_t *deriv = poly_new(10, (double[]){-12*gamma, 0, 0, 0, 0, 0, 6*gamma, 0, 0, 0, -2});
    double d2udr2, root, remainder;
    for(unsigned order = 11; order != 0; order--){
      root = poly_hm(deriv, 2, 1e-10, 1 << 14);
      remainder = poly_factor_out(deriv, root);
      if(fabs(remainder) > 1e-3){
        gamma += step_size;
        step_size /= 10;
        //report(INFO, "step size reduced to %e, reset to %e", step_size, gamma);
        break;
      }
      d2udr2 = second_deriv(gamma, c, root);
      if((root > 0) && d2udr2 > 0){
        printf("%.12f %f %e %+f\n", gamma, root, remainder, d2udr2);
        break;
      }
    }
    free(deriv);
  }
  return 0;
}
