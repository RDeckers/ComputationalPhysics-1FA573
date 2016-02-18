#include <integrators.h>
#include <common.h>

double Eulers_method(
  double (*dydx)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps /*we pass n_steps instead of h because for small h it will infinitely loop)*/
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0, y = y0;
  for(unsigned u = 0; u < n_steps; u++){
    y += step_size*dydx(x, y);
    x += step_size;
  }
  return y;
}

double Taylor_series(
  double (*dydx)(double, double),
  double (*d2ydx2)(double, double),
  double (*d2ydxy)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0, y = y0;
  for(unsigned u = 0; u < n_steps; u++){
    //We have an fabs here --------V so that is step_size is negative (i.e. x_end < x0) we properly integrate backwards.
    y += step_size*dydx(x, y)+0.5*fabs(step_size)*step_size*(d2ydx2(x, y)+dydx(x, y)*d2ydxy(x, y));
    x += step_size;
  }
  return y;
}

/*Only works for functions for which dydx is linear in y, that is g(x)*y*/
double implicit_method(
  double (*g)(double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0, y = y0;
  for(unsigned u = 0; u < n_steps; u++){
    y *= (1+0.5*g(x)*step_size)/(1-0.5*g(x+step_size)*step_size);
    x += step_size;
  }
  return y;
}

double rk_2(
  double (*f)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0, y = y0;
  for(unsigned u = 0; u < n_steps; u++){
    double k = step_size*f(x,y);
    y += step_size*f(x+0.5*step_size, y+0.5*k);
    x += step_size;
  }
  return y;
}

double rk_3(
  double (*f)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0, y = y0;
  for(unsigned u = 0; u < n_steps; u++){
    double k1 = step_size*f(x,y);
    double k2 = step_size*f(x+0.5*step_size, y+0.5*k1);
    double k3 = step_size*f(x+step_size, y-k1+2*k2);
    y += 1.0/6.0*(k1+4*k2+k3);
    x += step_size;
  }
  return y;
}

double rk_4(
  double (*f)(double, double),
  double x0, double y0,
  double x_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0, y = y0;
  for(unsigned u = 0; u < n_steps; u++){
    double k1 = step_size*f(x,y);
    double k2 = step_size*f(x+0.5*step_size, y+0.5*k1);
    double k3 = step_size*f(x+0.5*step_size, y+0.5*k2);
    double k4 = step_size*f(x+step_size, y+k3);
    y += 1.0/6.0*(k1+2*k2+2*k3+k4);
    x += step_size;
  }
  return y;
}

void rk_4_2d(
  double (*f)(double, double, double),
  double (*g)(double, double, double),
  double *x, double *y,
  double t0, double t_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return;
  }
  double step_size = subdivide(t0, t_end, n_steps);
  double X = *x, Y = *y, t = t0;
  for(unsigned u = 0; u < n_steps; u++){
    double k1 = step_size*f(t,X,Y);
    double l1 = step_size*g(t,X,Y);

    double k2 = step_size*f(t+0.5*step_size, X+0.5*k1, Y+0.5*l1);
    double l2 = step_size*g(t+0.5*step_size, X+0.5*k1, Y+0.5*l1);

    double k3 = step_size*f(t+0.5*step_size, X+0.5*k2, Y+0.5*l2);
    double l3 = step_size*g(t+0.5*step_size, X+0.5*k2, Y+0.5*l2);

    double k4 = step_size*f(t+step_size, X+k3, Y+l3);
    double l4 = step_size*g(t+step_size, X+k3, Y+l3);
    X += 1.0/6.0*(k1+2*k2+2*k3+k4);
    Y += 1.0/6.0*(l1+2*l2+2*l3+l4);
    t += step_size;
  }
  *x = X;
  *y = Y;
}

void Numerov_into_array(
  double (*S)(double),
  double (*k_sq)(double),
  double x0, double x_end,
  double *y,
  unsigned n_steps
){
  if(0 == n_steps){
    return;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double x = x0+step_size;
  double S_prev = S(x0), S_cur = S(x);
  double k_prev = k_sq(x0), k_cur = k_sq(x);
  for(unsigned u = 1; u < n_steps; u++){
    double k_next = k_sq(x+step_size);
    double S_next = S(x+step_size);
    double step_square_over_12 = step_size*step_size/12.0;
    double S_term = step_square_over_12*(S_next+10*S_cur+S_prev);
    double y_term = 2*(1-5*step_square_over_12*k_cur)*y[u];
    double y_prev_term = (1+step_square_over_12*k_prev)*y[u-1];
    //update x and y
    y[u+1] = (y_term-y_prev_term+S_term)/(1+step_square_over_12*k_next);
    x += step_size;
    //move the old variables around for the next iteration
    k_prev = k_cur;
    k_cur = k_next;
    S_prev = S_cur;
    S_cur =  S_next;
  }
}
double Numerov(
  double (*S)(double),
  double (*k_sq)(double),
  double x0, double y0, double y1,
  double x_end,
  unsigned n_steps
){
  if(0 == n_steps){
    return y0;
  }
  double step_size = subdivide(x0, x_end, n_steps);
  double y_prev = y0, y = y1, x = x0+step_size;
  double S_prev = S(x0), S_cur = S(x);
  double k_prev = k_sq(x0), k_cur = k_sq(x);
  for(unsigned u = 1; u < n_steps; u++){
    double k_next = k_sq(x+step_size);
    double S_next = S(x+step_size);
    double step_square_over_12 = step_size*step_size/12.0;
    double S_term = step_square_over_12*(S_next+10*S_cur+S_prev);
    double y_term = 2*(1-5*step_square_over_12*k_cur)*y;
    double y_prev_term = (1+step_square_over_12*k_prev)*y_prev;
    //no longer need y_prev
    y_prev = y;
    //update x and y
    y = (y_term-y_prev_term+S_term)/(1+step_square_over_12*k_next);
    x += step_size;
    //move the old variables around for the next iteration
    k_prev = k_cur;
    k_cur = k_next;
    S_prev = S_cur;
    S_cur =  S_next;
  }
  return y;
}
