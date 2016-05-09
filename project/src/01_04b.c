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

void approx_square_potential_root(double at, double a, double b, double c, double *base, double *slope){
  *base = 1-3*b*b/(at*at)+4*a*pow(c,6)*(-13*pow(c,6)+7*pow(at,6))/pow(at,12);
  *slope = (48*a*pow(c,12)-24*a*pow(c*at,6)+2*b*b*pow(at,10))/pow(at,13);
}

double f(double r, double a, double b, double c){
  double rsqinv = 1/(r*r);
  double cr6inv = (c*c*rsqinv)*(c*c*rsqinv)*(c*c*rsqinv);
  double cr12inv = cr6inv*cr6inv;
  return rsqinv/(sqrt(1.0-(b*b*rsqinv)-4*a*(cr12inv-cr6inv)));
}

double f_approx(double r, double at, double bt, double ab, double bb){
  return pow(at+r*bt,2)/sqrt(ab+r*bb);
}

typedef struct{
  double a;
  double b;
  double c;
  double at;
  double bt;
  double ab;
  double bb;
}f_args;

f_args create_approximated_fargs(double a, double b, double c, double root, double r_max){
  f_args args = {
    .a = a,
    .b = b,
    .c = c
  };
  approx_rinv(root, &(args.at), &(args.bt));
  approx_square_potential_root(root, a, b, c, &(args.ab), &(args.bb));
  return args;
}

double f_corr(double x, void* args){
  f_args *c_args = args;
  //report(WARN, "a = %f\nb = %f\nat = %f\nbt = %f\nab = %f\nbb = %f", c_args->a, c_args->b, c_args->at, c_args->bt, c_args->ab, c_args->bb);
  double eval = f(x,c_args->a,c_args->b, c_args->c) - f_approx(x,c_args->at,c_args->bt,c_args->ab,c_args->bb);
  return eval;
}

double exact_correction_integrand(double r, double at, double bt, double ab, double bb){
  return 2*sqrt(ab+bb*r)*(15*at*at*bb*bb+10*at*bb*bt*(-2*ab+bb*r)+bt*bt*(8*ab*ab-4*ab*bb*r+3*bb*bb*r*r))/(15*bb*bb*bb);
}

double exact_correction_integral(double r_min, double r_max, f_args *args){
  return exact_correction_integrand(r_max, args->at, args->bt, args->ab, args->bb);// - exact_correction_integrand(r_min, args->at, args->bt, args->ab, args->bb);
}

double exact_integrand(double r, double b){
  return -(1.0/b)*asin(b/r);
}


double realistically_approximated(double a, double b, double c){
  double r_max = 3*c;
  double result = 0;
  if(b < r_max){//only effected if we approach close enough
    result += exact_integrand(r_max, b) - exact_integrand(b, b);
    poly_t *poly = poly_new(12, (double[]){-4*a*pow(c,12),0,0,0,0,0,4*a*pow(c,6),0,0,0,-b*b,0,1});
    //report(INFO, "starting root-finding for %f %f %f",a ,b, c);
    double r_min = -1e20;
    for(unsigned order = 12; order != 0; order--){
      double root = poly_hm(poly, r_max, 1e-10, 1 << 14);
      double remainder = poly_factor_out(poly, root);
      if(fabs(remainder) > 1e-6){
        //report(WARN, "stopping at order %u because no root was found", order);
        break;
      }
      if(root > r_min){
        r_min = root;
      }
      //report(PASS, "Found root %f, remainder = %f", root, remainder);
    }
    // double r_min = fabs(poly_hm(poly, r_max, 1e-9));
    //return r_min;
    //report(PASS, "rmin = %f", r_min);
    f_args args = create_approximated_fargs(a, b, c, r_min,  r_max);
    gk_output_t out;
    g7k15_integrate(f_corr,&args,r_min,r_max, NULL, &out);
    result -= (out.val+exact_correction_integral(r_min, r_max, &args));
    free(poly);
  }
  result *= 2*b;
  return result;
}

int main(int argc, char** argv){
  double c = 1;
  for(double a = 1e2; a >= 1e-2; a /= 5){
    printf("\"ε/E = %.2e\"\n", a);
    for(double b = 1e-6; b <= 3*c; b += 3*c/100){
      double center = realistically_approximated(a, b, c);
      /*if(center != 0)*/{
        double h = 6*c/500;
        double left = realistically_approximated(a, b-h, c);
        double right = realistically_approximated(a, b+h, c);
        double cross_section = b/sin(center)*fabs((2*h)/(right-left));
        printf("%f %e\n", center, cross_section);
      }
      //report(PASS,"\t%e %e",a,b);
    }
    puts("\n");
    //report(PASS, "%e", a);
  }
  return 0;
}
