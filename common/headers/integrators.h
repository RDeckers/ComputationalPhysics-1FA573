double Eulers_method(
  double (*dydx)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps /*we pass n_steps instead of h because for small h it will infinitely loop)*/
);
double Taylor_series(
  double (*dydx)(double, double),
  double (*d2ydx2)(double, double),
  double (*d2ydxy)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
);

/*Only works for functions for which dydx is linear in y, that is g(x)*y*/
double implicit_method(
  double (*g)(double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
);

double rk_2(
  double (*dydx)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps /*we pass n_steps instead of h because for small h it will infinitely loop)*/
);
double rk_3(
  double (*dydx)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps /*we pass n_steps instead of h because for small h it will infinitely loop)*/
);
double rk_4(
  double (*dydx)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps /*we pass n_steps instead of h because for small h it will infinitely loop)*/
);
void rk_4_2d(
  double (*f)(double, double, double),
  double (*g)(double, double, double),
  double *x, double *y,
  double t0, double t_end,
  unsigned n_steps
);

double Numerov(
  double (*S)(double),
  double (*k_sq)(double),
  double x0, double y0, double y1,
  double x_end,
  unsigned n_steps
);

void Numerov_into_array(
  double (*S)(double),
  double (*k_sq)(double),
  double x0, double x_end,
  double *y,
  unsigned n_steps
);
