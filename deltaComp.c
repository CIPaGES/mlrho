/***** deltaComp.c ********************************
 * Description: Compute diesequilibrium coefficient,
 *   delta, for single diploid individual.
 * Reference: Lynch, M. (2008). Estimation of nuc-
 *   leotide diversity, disequilibrium coefficients,
 *   and mutation rates from high-coverage genome-
 *   sequencing projects. Mol. Biol. Evol. 25:
 *   2409-2419.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Mar 17 21:21:02 2009
 **************************************************/
#include <float.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "interface.h"
#include "profileTree.h"
#include "mlComp.h"
#include "eprintf.h"

Node *root;
double likelihood;
int numPara;
double globalPi, globalEpsilon;

double rhoFromDelta(double t, double d);
void lik(Node *p, double pi, double ee, double de);
double myF(const gsl_vector *v, void *params);
double confFun(double x, void *params);
void conf(Args *args, Result *result);
double iterate(Args *args, gsl_root_fsolver *s, double xLo, double xHi);

/* estimateDelta: estimate pi, delta, and epsilon using 
 * the Nelder-Mead Simplex algorithm; code adapted 
 * from Galassi, M., Davies, J., Theiler, J., Gough, B., 
 * Jungman, G., Booth, M., Rossi, F. (2005). GNU 
 * Scientific Library Reference Manual. Edition 1.6, 
 * for GSL Version 1.6, 17 March 2005, p 472f.
 */
void estimateDelta(Node *r, Args *args, Result *result, int np){
  const gsl_multimin_fminimizer_type *T =  gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter = 0;
  int status;
  double size;

  /*set up likelihood computation */
  numPara = np;
  root = r;
  iniMlComp(root,1);
  /* initialize vertex size vector */
  ss = gsl_vector_alloc(np);
  /* set all step sizes */
  gsl_vector_set_all(ss, args->s);
  /* starting point */
  x = gsl_vector_alloc(np);
  gsl_vector_set(x, 0, args->P);
  gsl_vector_set(x, 1, args->E);
  if(np == 3)
    gsl_vector_set(x, 2, args->D);
  /* initialize method and iterate */
  minex_func.f = &myF;
  minex_func.n = np;
  minex_func.params = (void *)NULL;
  s = gsl_multimin_fminimizer_alloc(T, np);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if(status)
      break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, args->t);

  }while(status == GSL_CONTINUE && iter < args->i);
  if(status)
    printf("WARNING: Estimation of \\Delta failed.\n");
  if(np == 1)
    result->de = gsl_vector_get(s->x, 0);
  else if(np == 2 || np == 3){
    result->pi = gsl_vector_get(s->x, 0);
    result->ee = gsl_vector_get(s->x, 1);
    if(np == 3)
      result->de = gsl_vector_get(s->x, 2);
  }
  result->l = s->fval;
  result->i = iter;
  if(!args->f){
    result->pi = globalPi;
    result->ee = globalEpsilon;
  }
  result->rhoFromDelta = rhoFromDelta(result->pi,result->de);
  conf(args, result);
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  freeMlComp();

}

/* myF: function called during the minimization procedure */
double myF(const gsl_vector* v, void *params){
  double pi, ee, de;

  if(numPara == 1){
    pi = globalPi;
    ee = globalEpsilon;
    de = gsl_vector_get(v, 0);
  }else{
    pi = gsl_vector_get(v, 0);
    ee = gsl_vector_get(v, 1);
    if(numPara == 3)
      de = gsl_vector_get(v, 2);
    else
      de = 0;
  }
  likelihood = 0;
  if(pi < 0 || ee < 0 || de < -1 || pi > 1 || ee > 1 || de > 1)
    return DBL_MAX;
  lik(root, pi, ee, de);
  return -likelihood;
}

/* lik: traverse profile tree and compute likelihood according to
 * equations (11) and (12) in Lynch (2008)
 */
void lik(Node *p, double pi, double ee, double de){
  double l, lOneA, lOneB, lTwoA, lTwoB;
  double h0, h2;

  if(p != NULL){
    lik(p->left, pi, ee, de);
    h0 = 1/(1+pi)/(1+pi)+de*pi/(1+pi)/(1+pi);
    h2 = pi*pi/(1+pi)/(1+pi)+de*pi/(1+pi)/(1+pi);
    lOneA = lOne(p->c1, p->profile1, ee);
    lOneB = lOne(p->c2, p->profile2, ee);
    lTwoA = lTwo(p->c1, p->profile1, ee);
    lTwoB = lTwo(p->c2, p->profile2, ee);
    l = h0*lOneA*lOneB+h2*lTwoA*lTwoB			\
      + 0.5*(1-h0-h2)*(lOneA*lTwoB+lTwoA*lOneB);
    if(l>0)
      likelihood += log(l) * p->n;
    else
      likelihood += log(DBL_MIN) * p->n;
    lik(p->right, pi, ee, de);
  }
}

/* conFun: called for confidence interval estimation of pi */
double confFun(double x, void *params){
  Result *res;
  double l;

  res = (Result *)params;
  likelihood = 0;
  if(res->type == 0)
    lik(root, x, res->ee, res->de);
  else if(res->type == 1){
    lik(root, res->pi, x, res->de);
  }else
    lik(root, res->pi, res->ee, x);
  l =  likelihood + res->l + 2;
  return l;
}

void conf(Args *args, Result *result){
  const gsl_root_fsolver_type *solverType;
  gsl_root_fsolver *s;
  gsl_function fun;
  double xLo, xHi;
  int status;

  /* preliminaries */
  gsl_set_error_handler_off();
  fun.function = &confFun;
  fun.params = result;
  solverType = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(solverType);
  /* search for lower bound of pi */
  result->type = 0;
  xLo = 0.0;
  xHi = result->pi;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->pLo = 0.0;
    printf("WARNING: Lower confidence limit of \theta cannot be estimated; setting it to 0.\n");
  }else
    result->pLo = iterate(args, s, xLo, xHi);
  /* search for upper bound of pi */
  xLo = result->pi;
  xHi = 0.75;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->pUp = 0.75;
    printf("WARNING: Upper confidence limit of \\theta cannot be estimated; setting it to 0.75.\n");
  }else
    result->pUp = iterate(args, s, xLo, xHi);
  /* search for lower bound of epsilon */
  result->type = 1;
  xLo = 0.0;
  xHi = result->ee;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->eLo = 0;
    printf("WARNING: Lower confidence limit of \\epsilon cannot be estimated; setting it to 0.\n");
  }else
    result->eLo = iterate(args, s, xLo, xHi);
  /* search for upper bound of epsilon */
  xLo = result->ee;
  xHi = 0.5;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->eUp = 1;
    printf("WARNING: Upper confidence limit of \\epsilon cannot be estimated; setting it to 1.\n");
  }else
    result->eUp = iterate(args, s, xLo, xHi);
  /* search for lower bound of delta */
  result->type = 2;
  if(result->de < 0)
    xLo = -1;
  else
    xLo = 0;
  xHi = result->de;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    printf("WARNING: Lower confidence limit of \\Delta cannot be estimated; setting it to -1.\n");
    result->dLo = -1.0;
  }else
    result->dLo = iterate(args, s, xLo, xHi);
  /* search for upper bound of delta */
  xLo = result->de;
  xHi = 1.0;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    printf("WARNING: Upper confidence limit of \\Delta cannot be estimated; setting it to 1.\n");
    result->dUp = 1;
  }else
    result->dUp = iterate(args, s, xLo, xHi);
  gsl_root_fsolver_free(s);
}

double iterate(Args *args, gsl_root_fsolver *s, double xLo, double xHi){
  int iter;
  double r;
  int status;

  iter = 0;
  do{
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    xLo = gsl_root_fsolver_x_lower(s);
    xHi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(xLo, xHi, 0, args->t);
  }while(status == GSL_CONTINUE && iter < args->i);
  
  return r;
}

void setPi(double pi){
  globalPi = pi;
}

void setEpsilon(double ee){
  globalEpsilon = ee;
}


double rhoFromDelta(double t, double d){
  double r;

  r = (t + pow(t,2) - 
       d*(13 + 19*t + 
	  6*pow(t,2)) + 
       sqrt(1 + t)*
       sqrt(pow(t,2)*(1 + t) + 
	    2*d*t*
	    (23 + 17*t + 
	     2*pow(t,2)) + 
	    pow(d,2)*
	    (97 + 109*t + 
	     32*pow(t,2) + 
	     4*pow(t,3))))/
    (2.*d*(1 + t));

  return r;
}
