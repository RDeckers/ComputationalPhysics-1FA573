#include <poly.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <utilities/logging.h>

poly_t* poly_new(unsigned n, double *coef){
  poly_t *poly = malloc(sizeof(poly_t)+(n+1)*sizeof(double));
  poly->n = n;
  if(coef){
    for(unsigned u = 0; u < n+1; u++){
      //report(INFO,"[%u, %f], ",u, coef[u]);
      poly->coef[u] = coef[u];
    }
    //puts("");
  }
  return poly;
}

poly_t* poly_new_deriv_from_poly(poly_t *base){
  if(base->n == 0){
    double zero = 0.0;
    return poly_new(0, &zero);
  }
  poly_t *deriv = poly_new(base->n-1, NULL);
  for(unsigned u = 0; u < base->n; u++){
    deriv->coef[u] = base->coef[u+1]*(u+1);
    //printf("[%u, %f], ",u, deriv->coef[u]);
  }
  //puts("");
  return deriv;
}

double poly_eval(poly_t *poly, double x){
  double res = poly->coef[poly->n];//last coef
  for(unsigned u = poly->n; u != 0; u--){
    res = res*x+poly->coef[u-1];
  }
  return res;
}

double poly_nr(poly_t *poly, double x0, double precision){
  poly_t *deriv = poly_new_deriv_from_poly(poly);
  for
  (
    double step = poly_eval(poly, x0)/poly_eval(deriv, x0);
    fabs(step) >= precision;
    step = poly_eval(poly, x0)/poly_eval(deriv, x0)
  ){
    //printf("%e, %e\t%+2.12f\n", poly_eval(poly, x0), poly_eval(deriv, x0), x0);
    x0 -= step;
  }
  free(deriv);
  return x0;
}

double poly_hm(poly_t *poly, double x0, double precision, unsigned max_steps){
  poly_t *deriv = poly_new_deriv_from_poly(poly);
  poly_t *second_deriv = poly_new_deriv_from_poly(deriv);
  double step;
  unsigned steps = 0;
  do{
    double f = poly_eval(poly, x0);
    double dfdx = poly_eval(deriv, x0);
    double d2fdx2 = poly_eval(second_deriv, x0);
    step = 2*f*dfdx/(2*dfdx*dfdx-f*d2fdx2);
    //report(INFO,"%e %e %e %e %e", x0, f, dfdx, d2fdx2, step);
    x0 -= step;
    steps++;
  }while((fabs(step) > precision) && ((x0+step) != x0) && (steps < max_steps));
  free(deriv);
  free(second_deriv);
  return x0;
}

double poly_factor_out(poly_t *poly, double p){
  if(poly->n == 0){
    return poly->coef[0];
  }
  //report(INFO, "factoring out %f from order %u", p, poly->n);
  double prev = poly->coef[poly->n];
  double remainder =0;
  double next;
  for(int n = poly->n-1; n >= 0; n--){
    next = poly->coef[n];
    poly->coef[n] = (remainder*p)+prev;
    //report(PASS, "updated coef %u to %f", n, poly->coef[n]);
    prev = next;
    remainder = poly->coef[n];
  }
  poly->n = poly->n - 1;
  remainder = remainder*p+next;
  //report(remainder ? WARN : PASS, "remainder = %f", remainder);
  return remainder;
}
