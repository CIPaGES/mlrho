/***** piComp.c ************************************
 * Description: Compute heterozygosity, pi, for
 *   single diploid individual.
 * Reference: Lynch, M. (2008). Estimation of nuc-
 *   leotide diversity, disequilibrium coefficients,
 *   and mutation rates from high-coverage genome-
 *   sequencing projects. Mol. Biol. Evol. 25:
 *   2409-2419.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Mar 17 21:20:57 2009
 ***************************************************/
#include <float.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <assert.h>
#include "interface.h"
#include "profileTree.h"
#include "mlComp.h"
#include "eprintf.h"

Node *root;
double likelihood;

void likP(Node *p, double pi, double ee);
double myP(const gsl_vector *v, void *params);
double myPconf(double x, void *params);
double myEconf(double x, void *params);
void confP(Args *args, Result *result);
void confE(Args *args, Result *result);


/* estimatePi: estimate pi and epsilon using the Nelder-Mead
 * Simplex algorithm; code adapted from Galassi, M., Davies, 
 * J., Theiler, J., Gough, B., Jungman, G., Booth, M., 
 * Rossi, F. (2005). GNU Scientific Library Reference Manual. 
 * Edition 1.6, for GSL Version 1.6, 17 March 2005, p 472f.
 */
void estimatePi(Node *r, Args *args, Result *result){
  size_t np = 2;
  const gsl_multimin_fminimizer_type *T =  gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter = 0;
  int status;
  double size;

  /*set up likelihood computation */
  root = r;
  iniMlComp(root,0);
  /* initialize vertex size vector */
  ss = gsl_vector_alloc(np);
  /* set all step sizes */
  gsl_vector_set_all(ss, args->s);
  /* starting point */
  x = gsl_vector_alloc(np);
  gsl_vector_set(x, 0, args->P);
  gsl_vector_set(x, 1, args->E);
  /* initialize method and iterate */
  minex_func.f = &myP;
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
  if(GSL_CONTINUE)
    assert("ERROR[mlRho]: increase the number of iterations using option -i.\n");

  result->pi = gsl_vector_get(s->x, 0);
  result->ee = gsl_vector_get(s->x, 1);
  result->l = s->fval;
  result->i = iter;
  gsl_set_error_handler_off();
  confP(args, result);
  confE(args, result);
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  freeMlComp();

}
/* myF: function called during the minimization procedure */
double myP(const gsl_vector* v, void *params){
  double pi, ee;

  pi = gsl_vector_get(v, 0);
  ee = gsl_vector_get(v, 1);
  likelihood = 0;
  if(pi < 0 || ee < 0 || pi > 1 || ee > 1){
    return DBL_MAX;
  }
  likP(root, pi, ee);
  return -likelihood;
}

/* lik: traverse profile tree and compute likelihood */
void likP(Node *p, double pi, double ee){
  double l;

  /* printf("likP - p->c1: %d; pi: %f\n",p->c1, pi); */
  
  if(p != NULL){
    likP(p->left, pi, ee);

    l = lOne(p->c1, p->profile, ee) * (1.0 - pi)	\
      + lTwo(p->c1, p->profile, ee) * pi;
    if(l>0){
      likelihood += (log(l) * p->n);
    }

    likP(p->right, pi, ee);
  }
}

/* myPconf: called for confidence interval estimation of pi */
double myPconf(double x, void *params){
  Result *res;

  res = (Result *)params;
  likelihood = 0;
  likP(root, x, res->ee);

  return likelihood + res->l + 2;
}

/* myEconf: called for confidence interval estimation of epsilon */
double myEconf(double x, void *params){
  Result *res;
  double l;

  res = (Result *)params;
  likelihood = 0;
  likP(root, res->pi, x);

  l =  likelihood + res->l + 2;
  return l;
}

void confE(Args *args, Result *result){
  int status;
  int iter = 0;
  int max_iter = args->i;
  const gsl_root_fsolver_type *solverType;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function fun;
  double xLo, xHi;

  fun.function = &myEconf;
  fun.params = result;

  solverType = gsl_root_fsolver_brent;

  s = gsl_root_fsolver_alloc(solverType);
  /* search for upper bound */
  xLo = result->ee;
  xHi = 0.5;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status)
    result->eUp = 0.5;
  else{
    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r =gsl_root_fsolver_root(s);
      xLo = gsl_root_fsolver_x_lower(s);
      xHi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(xLo, xHi, 0, args->t);
    }while(status == GSL_CONTINUE && iter < max_iter);
    result->eUp = r;
  }
  /* search for lower bound */
  xLo = 0.0;
  xHi = result->ee;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status)
    result->eLo = 0;
  else{
    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r =gsl_root_fsolver_root(s);
      xLo = gsl_root_fsolver_x_lower(s);
      xHi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(xLo, xHi, 0, args->t);
    }while(status == GSL_CONTINUE && iter < max_iter);
    result->eLo = r;
  }
  gsl_root_fsolver_free(s);
}


void confP(Args *args, Result *result){
  int status;
  int iter = 0;
  int max_iter = args->i;
  const gsl_root_fsolver_type *solverType;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function fun;
  double xLo, xHi;

  fun.function = &myPconf;
  fun.params = result;

  solverType = gsl_root_fsolver_brent;

  s = gsl_root_fsolver_alloc(solverType);
  /* search for upper bound */
  xLo = result->pi;
  xHi = 0.75;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status)
    result->pUp = 0.75;
  else{
    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r =gsl_root_fsolver_root(s);
      xLo = gsl_root_fsolver_x_lower(s);
      xHi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(xLo, xHi, 0, 0.001);
    }while(status == GSL_CONTINUE && iter < max_iter);
    result->pUp = r;
  }
  /* search for lower bound */
  xLo = 0.0;
  xHi = result->pi;
  status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
  if(status)
    result->pLo = 0;
  else{
    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r =gsl_root_fsolver_root(s);
      xLo = gsl_root_fsolver_x_lower(s);
      xHi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(xLo, xHi, 0, 0.001);
    }while(status == GSL_CONTINUE && iter < max_iter);
    result->pLo = r;
  }
  gsl_root_fsolver_free(s);
}
