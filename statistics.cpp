/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "statistics.h"
#include "pw_gv.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>

using namespace std;

double getNegativeBinomial(int k, double p, double n)
{
  return gsl_ran_negative_binomial_pdf(k, p, n);
}

double getZINB(int k, double p, double n, double p0)
{
  double r;
  if(!k) {
    r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
  }else{
    r = (1 - p0) * gsl_ran_negative_binomial_pdf(k, p, n);
  }
  return r;
}

gsl_multimin_fminimizer *gsl_multimin_fminimizer_new(size_t ndim)
{
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;  /* ネルダーとミードのシンプレックス法 */
  gsl_multimin_fminimizer            *s = gsl_multimin_fminimizer_alloc(T, ndim); 
  return s;
}

gsl_vector *gsl_vector_new(int ndim, double init)
{
  gsl_vector *p = gsl_vector_alloc(ndim);
  gsl_vector_set_all(p, init);
  return p;
}

double f_zinb_const(const gsl_vector *v, void *params)
{
  double p, n, r, fxy(0), p0;
  double *par = (double *)params;
  p = gsl_vector_get(v, 0);
  if(p <= 0) p = 0.01;
  if(p >= 1) p = 0.99;
  n = gsl_vector_get(v, 1);
  p0 = gsl_vector_get(v, 2);
  if(p0 < 0) p0 = 0;
  if(p0 > 1) p0 = 1.0;
  //  LOG("zinb_const p=%f, n=%f, p0=%f\n", p, n, p0);
  
  for(int i=0; i<NUM_WIGDISTARRAY; i++){
    if(!i) r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
    else   r =      (1 - p0) * gsl_ran_negative_binomial_pdf(i, p, n);
    //    printf("%d: %f - %f\n", i, par[i] , r);
    fxy += (par[i] - r)*(par[i] - r);
  }
  return fxy;
}

void func_iteration(gsl_multimin_fminimizer *s)
{
  size_t iter(0);
  int status;
  double size;
  do {
    ++iter;
    status = gsl_multimin_fminimizer_iterate(s);  /* 繰り返し計算を1回行う。予期しない問題が発生した場合はエラーコードを返す。 */
    if(status) break;
    size = gsl_multimin_fminimizer_size(s);       /* sのその時点での最小点の最良推定値を返す */
    status = gsl_multimin_test_size(size, 1e-3);  /* sizeが閾値(1e-3)より小さければGSL_SUCCESS を、そうでなければGSL_CONTINUEを返す。 */
    if(status == GSL_SUCCESS) { cout << "converged to minimum at" << endl; }
#ifdef DEBUG   
    BPRINT("%1% p=%2% n=%3% p0=%4% f() = %5% size = %6%\n") % iter % gsl_vector_get(s->x, 0) % gsl_vector_get(s->x, 1) % gsl_vector_get(s->x, 2) % gsl_multimin_fminimizer_minimum(s) % size;
#endif
  } while(status == GSL_CONTINUE && iter < 1000);

  return;
}

void estimateZINB(Mapfile &p)
{
  size_t ndim(3);

  // initialization 
  int thre = NUM_WIGDISTARRAY;
  double par[thre];
  for(int i=0; i<thre; ++i) {
    if(!p.genome.wigDist[i]) par[i] = 0;
    else par[i] = p.genome.wigDist[i] /(double)p.genome.nbin;
#ifdef DEBUG   
    BPRINT("par[%1%]=%2%, darray_bg=%3%, num=%4%\n") % i % par[i] % p.genome.wigDist[i] % p.genome.nbin;
#endif
  }
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_new(ndim);
  gsl_vector *x = gsl_vector_alloc(ndim);
  gsl_vector_set(x, 0, p.lchr->nb_p);
  gsl_vector_set(x, 1, p.lchr->nb_n);
  gsl_vector_set(x, 2, 0.02); // p0
  gsl_vector *ss = gsl_vector_new(ndim, 0.1); // step size 

  gsl_multimin_function minex_func;
  minex_func.n = ndim;
  minex_func.f = &f_zinb_const;
  minex_func.params = (void *)&par;
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  func_iteration(s);

  p.genome.nb_p  = gsl_vector_get(s->x, 0);
  if(p.genome.nb_p <= 0) p.genome.nb_p = 0.01;
  if(p.genome.nb_p >= 1) p.genome.nb_p = 0.99;
  p.genome.nb_n  = gsl_vector_get(s->x, 1);
  p.genome.nb_p0 = gsl_vector_get(s->x, 2);
  if(p.genome.nb_p0 < 0) p.genome.nb_p0 = 0.0;
  if(p.genome.nb_p0 > 1) p.genome.nb_p0 = 1.0;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return;
}
