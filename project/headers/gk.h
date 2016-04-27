
typedef struct{
  double abs_err;
  double rel_err;
}gk_settings_t;

typedef struct{
  double val;
  double err_estimate;
  unsigned n_evals;
}gk_output_t;

void g7k15_integrate(double (*f)(double, void*), void* args, double a, double b, gk_settings_t* settings, gk_output_t* output);
