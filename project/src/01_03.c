#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <poly.h>
#include <gk.h>
#include <utilities/logging.h>

void approx_rinv(double at, double *base, double *slope){
  *base = 2/at;
  *slope = -1/(at*at);
}

void approx_square_potential_root(double at, double a, double b, double *base, double *slope){
  *base = 1-a-3*b*b/(at*at);
  *slope = 2*b*b/(at*at*at);
}

double f(double r, double a, double b){
  return 1.0/(r*r*sqrt(1.0-(b*b)/(r*r)-a));
}

double f_approx(double r, double at, double bt, double ab, double bb){
  return pow(at+r*bt,2)/sqrt(ab+r*bb);
}

typedef struct{
  double a;
  double b;
  double at;
  double bt;
  double ab;
  double bb;
}f_args;

f_args create_approximated_fargs(double a, double b, double root, double r_max){
  f_args args = {
    .a = a,
    .b = b
  };
  approx_rinv(root, &(args.at), &(args.bt));
  approx_square_potential_root(root, a, b, &(args.ab), &(args.bb));
  return args;
}

f_args create_fargs(double a, double b){
  f_args args = {
    .a = a,
    .b = b,
    .at = 2*sqrt(1-a)/b,
    .bt = (a-1)/(b*b),
    .ab = 2*(a-1),
    .bb = 2*pow(1-a, 1.5)/b
  };
  return args;
}

double f_corr(double x, void* args){
  f_args *c_args = args;
  //report(WARN, "a = %f\nb = %f\nat = %f\nbt = %f\nab = %f\nbb = %f", c_args->a, c_args->b, c_args->at, c_args->bt, c_args->ab, c_args->bb);
  double eval = f(x,c_args->a,c_args->b) - f_approx(x,c_args->at,c_args->bt,c_args->ab,c_args->bb);
  return eval;
}

double exact_correction_integrand(double r, double at, double bt, double ab, double bb){
  return 2*sqrt(ab+bb*r)*(15*at*at*bb*bb+10*at*bb*bt*(-2*ab+bb*r)+bt*bt*(8*ab*ab-4*ab*bb*r+3*bb*bb*r*r))/(15*bb*bb*bb);
}

double exact_correction_integral(double r_min, double r_max, f_args *args){
  //report(WARN, "a = %f\nb = %f\nat = %f\nbt = %f\nab = %f\nbb = %f", args->a, args->b, args->at, args->bt, args->ab, args->bb);
  return exact_correction_integrand(r_max, args->at, args->bt, args->ab, args->bb);// - exact_correction_integrand(r_min, args->at, args->bt, args->ab, args->bb);
}

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

double first_approximation(double a, double b, double r_max){
  double result = 0;
  if(b < r_max){//only effected if we approach close enough
    gk_output_t out;
    f_args args = create_fargs(0, b);
    g7k15_integrate(f_corr,&args,b,r_max, NULL, &out);
    result = 2*b*(out.val+exact_correction_integral(b,r_max, &args));
  }
  return result;
}

double second_approximation(double a, double b, double r_max){
  double result = 0;
  double r_min = b/sqrt(1-a);
  if(b < r_max){//only effected if we approach close enough
    if((r_min < r_max) && (1-a > 0)){
      gk_output_t out;
      f_args args = create_fargs(a, b);
      g7k15_integrate(f_corr,&args,r_min,r_max, NULL, &out);
      result = 2*b*(out.val+exact_correction_integral(r_min, r_max, &args));
    }
  }
  return result;
}

double fully_approximated(double a, double b, double r_max){
  return first_approximation(a, b, r_max) - second_approximation(a, b, r_max);
}
double partially_approximated(double a, double b, double r_max){
  double result = 0;
  if(b < r_max){//only effected if we approach close enough
    result += exact_integrand(r_max, 0, b) - exact_integrand(b, 0, b);
  }
  result *= 2*b;
  return result - second_approximation(a, b, r_max);
}

double realistically_approximated(double a, double b, double r_max){
  double result = 0;
  if(b < r_max){//only effected if we approach close enough
    result += exact_integrand(r_max, 0, b) - exact_integrand(b, 0, b);
    poly_t *poly = poly_new(1, (double[]){-b*b,1-a});
    double r_min = sqrt(poly_nr(poly, r_max, 1e-9));
    //report(INFO, "rmin = %f / %f", r_min, b/sqrt(1-a));
    f_args args = create_approximated_fargs(a, b,r_min,  r_max);
    if(b < r_max){//only effected if we approach close enough
      if((r_min < r_max) && (1-a > 0)){
        gk_output_t out;
        g7k15_integrate(f_corr,&args,r_min,r_max, NULL, &out);
        result -= (out.val+exact_correction_integral(r_min, r_max, &args));
      }
    }
    free(poly);
  }
  result *= 2*b;
  return result;// - second_approximation(a, b, r_max);
}

int main(int argc, char** argv){
  double a = -0.5;
  double r_max = 1;
  for(double b = 0.05; b < 1; b += 0.05){
    double exact = exact_integral(a, b, r_max);
    double approximation_1 = partially_approximated(a, b, r_max);
    double approximation_2 = fully_approximated(a, b, r_max);
    double approximation_3 = realistically_approximated(a, b, r_max);
    printf("%f %e %e %e %e\n", b, exact, fabs(approximation_1-exact), fabs(approximation_2-exact), fabs(approximation_3-exact));
  }
  return 0;
}
