typedef struct{
  int n;
  double coef[0];
}poly_t;

poly_t* poly_new(unsigned n, double *coef);
poly_t* poly_new_deriv_from_poly(poly_t *base);
double poly_eval(poly_t *poly, double x);
double poly_nr(poly_t *poly, double x0, double precision);
double poly_hm(poly_t *poly, double x0, double precision,unsigned max_steps);
double poly_factor_out(poly_t *poly, double p);
