#include <gk.h>
#include <stdlib.h>
#include <math.h>
#include <utilities/logging.h>

const double g30k61_node[] = {
  0.999484410050490637571325895705811,
  0.996893484074649540271630050918695,
  0.991630996870404594858628366109486,
  0.983668123279747209970032581605663,
  0.973116322501126268374693868423707,
  0.960021864968307512216871025581798,
  0.944374444748559979415831324037439,
  0.926200047429274325879324277080474,
  0.905573307699907798546522558925958,
  0.882560535792052681543116462530226,
  0.857205233546061098958658510658944,
  0.829565762382768397442898119732502,
  0.799727835821839083013668942322683,
  0.767777432104826194917977340974503,
  0.733790062453226804726171131369528,
  0.697850494793315796932292388026640,
  0.660061064126626961370053668149271,
  0.620526182989242861140477556431189,
  0.579345235826361691756024932172540,
  0.536624148142019899264169793311073,
  0.492480467861778574993693061207709,
  0.447033769538089176780609900322854,
  0.400401254830394392535476211542661,
  0.352704725530878113471037207089374,
  0.304073202273625077372677107199257,
  0.254636926167889846439805129817805,
  0.204525116682309891438957671002025,
  0.153869913608583546963794672743256,
  0.102806937966737030147096751318001,
  0.051471842555317695833025213166723,
  0.000000000000000000000000000000000
};
const double k61_weights[] = {
    0.001389013698677007624551591226760,
    0.003890461127099884051267201844516,
    0.006630703915931292173319826369750,
    0.009273279659517763428441146892024,
    0.011823015253496341742232898853251,
    0.014369729507045804812451432443580,
    0.016920889189053272627572289420322,
    0.019414141193942381173408951050128,
    0.021828035821609192297167485738339,
    0.024191162078080601365686370725232,
    0.026509954882333101610601709335075,
    0.028754048765041292843978785354334,
    0.030907257562387762472884252943092,
    0.032981447057483726031814191016854,
    0.034979338028060024137499670731468,
    0.036882364651821229223911065617136,
    0.038678945624727592950348651532281,
    0.040374538951535959111995279752468,
    0.041969810215164246147147541285970,
    0.043452539701356069316831728117073,
    0.044814800133162663192355551616723,
    0.046059238271006988116271735559374,
    0.047185546569299153945261478181099,
    0.048185861757087129140779492298305,
    0.049055434555029778887528165367238,
    0.049795683427074206357811569379942,
    0.050405921402782346840893085653585,
    0.050881795898749606492297473049805,
    0.051221547849258772170656282604944,
    0.051426128537459025933862879215781,
    0.051494729429451567558340433647099
};
const double g30_weights[] = {
    0.007968192496166605615465883474674,
    0.018466468311090959142302131912047,
    0.028784707883323369349719179611292,
    0.038799192569627049596801936446348,
    0.048402672830594052902938140422808,
    0.057493156217619066481721689402056,
    0.065974229882180495128128515115962,
    0.073755974737705206268243850022191,
    0.080755895229420215354694938460530,
    0.086899787201082979802387530715126,
    0.092122522237786128717632707087619,
    0.096368737174644259639468626351810,
    0.099593420586795267062780282103569,
    0.101762389748405504596428952168554,
    0.102852652893558840341285636705415
};

const double g7k15_node[] ={
  0.991455371120812639206854697526329,
  0.949107912342758524526189684047851,
  0.864864423359769072789712788640926,
  0.741531185599394439863864773280788,
  0.586087235467691130294144838258730,
  0.405845151377397166906606412076961,
  0.207784955007898467600689403773245,
  0.000000000000000000000000000000000
};

const double k15_weights[] ={
  0.022935322010529224963732008058970,
  0.063092092629978553290700663189204,
  0.104790010322250183839876322541518,
  0.140653259715525918745189590510238,
  0.169004726639267902826583426598550,
  0.190350578064785409913256402421014,
  0.204432940075298892414161999234649,
  0.209482141084727828012999174891714
};

const double g7_weights[] ={
  0.129484966168869693270611432679082,
  0.279705391489276667901467771423780,
  0.381830050505118944950369775488975,
  0.417959183673469387755102040816327
};

typedef struct{
  double error;
  double val;
  double a;
  double b;
} gk_step_out_t;

gk_settings_t DEFAULT_GK_SETTINGS = {
  .abs_err = 1e-9, //TODO: use 10 & 21 instead?
  .rel_err = 0//1e-12
};

gk_step_out_t gk_step(double(*f)(double,void*), void* args, double a, double b, const double *nodes, const double *wg, const double *wk, const unsigned n){
  double evaluations[n];

  for(unsigned u = 0; u < n/2; u++){
    double node = nodes[u];
    evaluations[2*u+0] = f(0.5*((b-a)*node+(a+b)), args);
    evaluations[2*u+1] = f(0.5*((b-a)*(-node)+(a+b)), args);
  }
  evaluations[n-1] = f(0.5*(a+b), args);

  double konrod = 0;
  for(unsigned u = 0; u < n-1; u++){
    konrod += evaluations[u]*wk[u/2];
  }
  konrod +=evaluations[n-1]*wk[n/2];

  double mean = 0.5*konrod;
  konrod *= (b-a)/2;

  double asc = fabs(evaluations[n-1]-mean)*wk[n/2];
  for(unsigned u = 0; u < n-1; u++){
    asc += fabs(evaluations[u]-mean)*wk[u/2];
  }
  asc *= (b-a)/2;

  double gauss = 0;
  for(unsigned u = 0; u < n/4; u++){
    gauss += evaluations[4*u+2]*wg[u];//2,6,10
    gauss += evaluations[4*u+3]*wg[u];//3,7,11
  }
  gauss +=evaluations[n-1]*wg[n/4];
  gauss *= (b-a)/2;
  double scale = pow(200*fabs(gauss-konrod)/asc,1.5);
  double error;
  if(scale < 1)
    error = asc * scale;
  else
    error = asc;
  //report(INFO,"g: %e\tk: %e\t Error: %e (%e)\t @[%e, %e]", gauss, konrod, error, error/konrod, a, b);
  gk_step_out_t out = {.val = konrod, .error = error, .a = a, .b = b};
  return out;
 }

gk_step_out_t g7k15_step(double (*f)(double, void*), void* args, double a, double b){
  return gk_step(f, args, a, b, g7k15_node, g7_weights, k15_weights, 15);
}

gk_step_out_t g30k61_step(double (*f)(double, void*), void* args, double a, double b){
  return gk_step(f, args, a, b, g30k61_node, g30_weights, k61_weights, 61);
}

void g7k15_integrate(double (*f)(double, void*), void* args, double a, double b, gk_settings_t* settings, gk_output_t* output){
  if(!settings){
    settings = &DEFAULT_GK_SETTINGS;
  }
  unsigned buffers_used = 1;
  unsigned buffer_size = 32;

  gk_step_out_t *segments = malloc(buffer_size*sizeof(gk_step_out_t));//TODO: track global error estimate instead
  segments[0] = g7k15_step(f, args, a, b);
  double total_error = segments[0].error;
  double total_estimate = segments[0].val;
  double largest_error = total_error;
  unsigned largest_error_index = 0;
  while( (total_error > settings->abs_err)  && (fabs(total_error/total_estimate) > settings->rel_err)){
    for(unsigned u = 0; u < buffers_used;u++){
      double error = segments[u].error;
      if(error > largest_error){
        largest_error = error;
        largest_error_index = u;
      }
    }
    gk_step_out_t segment = segments[largest_error_index];

    double midpoint = 0.5*(segment.a+segment.b);
    if(buffer_size == buffers_used){
      buffer_size *=2;
      segments = realloc(segments,buffer_size*sizeof(gk_step_out_t));
    }
    gk_step_out_t right = g7k15_step(f, args, midpoint, segment.b);
    gk_step_out_t left = g7k15_step(f, args, segment.a, midpoint);
    if(!isfinite(left.val) || !isfinite(right.val)){
      report(WARN,"integrand is no longer finite!");//TODO: be smarter here.
      break;
    }
    segments[buffers_used++] = right;
    segments[largest_error_index] = left;
    total_error += left.error;
    total_estimate += left.val;
    total_error += right.error;
    total_estimate += right.val;
    total_error -= segment.error;
    total_estimate -= segment.val;

    largest_error = left.error;
  }

  free(segments);
  output->val = total_estimate;
  output->err_estimate = total_error;
  output->n_evals = 61*buffers_used;
  return;
}
