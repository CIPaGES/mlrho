/***** rhoComp.c **********************************
 * Description: Compute recombination parameter,
 *   rho, for single diploid individual.
 * Reference: Haubold et al. (2009).
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

void rhoLik(Node *p, double pi, double ee, double de);
double rhoMyF(const gsl_vector *v, void *params);
double rhoConfFun(double x, void *params);
void rhoConf(Args *args, Result *result);
double rhoIterate(Args *args, gsl_root_fsolver *s, double xLo, double xHi);

/* estimateRho: estimate pi, rho, and epsilon using 
 * the Nelder-Mead Simplex algorithm; code adapted 
 * from Galassi, M., Davies, J., Theiler, J., Gough, B., 
 * Jungman, G., Booth, M., Rossi, F. (2005). GNU 
 * Scientific Library Reference Manual. Edition 1.6, 
 * for GSL Version 1.6, 17 March 2005, p 472f.
 */
int estimateRho(Node *r, Args *args, Result *result, int np){
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
    gsl_vector_set(x, 2, args->R);
  /* initialize method and iterate */
  minex_func.f = &rhoMyF;
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
  /* if(status) */
  /*   printf("WARNING: \\rho estimation failed.\n"); */
  if(np == 1)
    result->rh = gsl_vector_get(s->x, 0);
  else if(np == 2 || np == 3){
    result->pi = gsl_vector_get(s->x, 0);
    result->ee = gsl_vector_get(s->x, 1);
    if(np == 3)
      result->rh = gsl_vector_get(s->x, 2);
  }
  result->l = s->fval;
  result->i = iter;
  if(!args->f){
    result->pi = globalPi;
    result->ee = globalEpsilon;
  }
  if(status == GSL_SUCCESS)
    rhoConf(args, result);
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  freeMlComp();
  return status;
}

/* myF: function called during the minimization procedure */
double rhoMyF(const gsl_vector* v, void *params){
  double pi, ee, rh;

  if(numPara == 1){
    pi = globalPi;
    ee = globalEpsilon;
    rh = gsl_vector_get(v, 0);
  }else{
    pi = gsl_vector_get(v, 0);
    ee = gsl_vector_get(v, 1);
    if(numPara == 3)
      rh = gsl_vector_get(v, 2);
    else
      rh = 0;
  }
  likelihood = 0;
  if(pi < 0 || ee < 0 || pi > 1 || ee > 1)
    return DBL_MAX;
  rhoLik(root, pi, ee, rh);
  return -likelihood;
}

/* lik: traverse profile tree and compute likelihood according to
 * equations (11) and (12) in Lynch (2008)
 */
void rhoLik(Node *p, double t, double ee, double r){
  double l, lOneA, lOneB, lTwoA, lTwoB;
  double h0, h2;
  double r2, t2, t3, de;
  int profile1[4], profile2[4];
  int i, c1, c2;

  r2 = r*r;
  t2 = t*t;
  t3 = t2*t;
  
  de = t*(18.+r+18.*t+r*t+4.*t2)/ \
    (18.+13.*r+r2+54.*t+40.*t2+8.*t3+r*(r*t+19.*t+6.*t2));

  if(p != NULL){
    rhoLik(p->left, t, ee, r);
    sscanf(p->key,"%d %d %d %d %d %d %d %d",&profile1[0],&profile1[1],&profile1[2],&profile1[3], \
	   &profile2[0],&profile2[1],&profile2[2],&profile2[3]);
    c1 = 0; 
    c2 = 0;
    for(i=0;i<4;i++){
      c1 += profile1[i];
      c2 += profile2[i];
    }
    h0 = 1/(1+t)/(1+t)+de*t/(1+t)/(1+t);
    h2 = t2/(1+t)/(1+t)+de*t/(1+t)/(1+t);
    lOneA = lOne(c1,profile1, ee);
    lOneB = lOne(c2, profile2, ee);
    lTwoA = lTwo(c1, profile1, ee);
    lTwoB = lTwo(c2, profile2, ee);
    l = h0*lOneA*lOneB+h2*lTwoA*lTwoB			\
      + 0.5*(1-h0-h2)*(lOneA*lTwoB+lTwoA*lOneB);
    if(l>0)
      likelihood += log(l) * p->n;
    
    rhoLik(p->right, t, ee, r);
  }
}

/* conFun: called for confidence interval estimation of pi */
double rhoConfFun(double x, void *params){
  Result *res;
  double l;

  res = (Result *)params;
  likelihood = 0;
  if(res->type == 0)
    rhoLik(root, x, res->ee, res->rh);
  else if(res->type == 1){
    rhoLik(root, res->pi, x, res->rh);
  }else
    rhoLik(root, res->pi, res->ee, x);
  l =  likelihood + res->l + 2;
  return l;
}

void rhoConf(Args *args, Result *result){
  const gsl_root_fsolver_type *solverType;
  gsl_root_fsolver *s;
  gsl_function fun;
  double xLo, xHi;
  int status;

  /* preliminaries */
  gsl_set_error_handler_off();
  fun.function = &rhoConfFun;
  fun.params = result;
  solverType = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(solverType);
  /* search for lower bound of pi */
  result->type = 0;
  xLo = 0.0;
  xHi = result->pi;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    printf("WARNING: Estimation of lower confidence limit of \\theta failed; setting it to 0.\n");
    result->pLo = 0.0;
  }else
    result->pLo = rhoIterate(args, s, xLo, xHi);
  /* search for upper bound of pi */
  xLo = result->pi;
  xHi = 0.75;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->pUp = 0.75;
    printf("WARNING: Estimation of upper confidence limit of \\theta failed; setting it to 0.75.\n");
  }else
    result->pUp = rhoIterate(args, s, xLo, xHi);
  /* search for lower bound of epsilon */
  result->type = 1;
  xLo = 0.0;
  xHi = result->ee;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->eLo = 0;
    printf("WARNING: Estimation of lower confidence limit of \\epsilon failed; setting it to 0.\n");
  }else
    result->eLo = rhoIterate(args, s, xLo, xHi);
  /* search for upper bound of epsilon */
  xLo = result->ee;
  xHi = 0.5;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->eUp = 1;
    printf("WARNING: Estimation of upper confidence limit of \\epsilon failed; setting it to 1.\n");
  }else
    result->eUp = rhoIterate(args, s, xLo, xHi);
  /* search for lower bound of rho */
  result->type = 2;
  xLo = 0;
  xHi = result->rh;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->rLo = 0;
    printf("WARNING: Estimation of lower confidence limit of \\rho failed; setting it to 0.\n");
  }else
    result->rLo = rhoIterate(args, s, xLo, xHi);
  /* search for upper bound of rho */
  xLo = result->rh;
  xHi = FLT_MAX;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status){
    result->rUp = 1;
    printf("WARNING: Estimation of upper confidence limit of \\rho failed; setting it to 1.\n");
  }else{
    result->rUp = rhoIterate(args, s, xLo, xHi);
  }
  gsl_root_fsolver_free(s);
}

double rhoIterate(Args *args, gsl_root_fsolver *s, double xLo, double xHi){
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

void rhoSetPi(double pi){
  globalPi = pi;
}

void rhoSetEpsilon(double ee){
  globalEpsilon = ee;
}
