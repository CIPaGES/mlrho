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
#include <assert.h>
#include "interface.h"
#include "ld.h"
#include "mlComp.h"
#include "eprintf.h"
#include "profile.h"

double globalPi, globalEpsilon;
int globalDist;
Node **globalProfilePairs;
int globalNumProfiles;
double likelihood;

double rhoFromDelta(double t, double d);
void lik(double de);
double myF(const gsl_vector *v, void *params);
double confFun(double x, void *params);
void conf(Args *args, Result *result);
double iterate(Args *args, gsl_root_fsolver *s, double xLo, double xHi);
void traverse(int a, Node *np, double h0, double h2, double complementHalf);

/* estimateDelta: estimate pi, delta, and epsilon using 
 * the Nelder-Mead Simplex algorithm; code adapted 
 * from Galassi, M., Davies, J., Theiler, J., Gough, B., 
 * Jungman, G., Booth, M., Rossi, F. (2005). GNU 
 * Scientific Library Reference Manual. Edition 1.6, 
 * for GSL Version 1.6, 17 March 2005, p 472f.
 */
Result *estimateDelta(Node **profilePairs, int numProfiles, Args *args, Result *result, int dist){
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter;
  int status, numPara;
  double size;

  globalProfilePairs = profilePairs;
  globalNumProfiles = numProfiles;
  T =  gsl_multimin_fminimizer_nmsimplex;
  s = NULL;
  iter = 0;
  numPara = 1; /* one parameter estimation */
  globalPi = result->pi;
  globalEpsilon = result->ee;
  globalDist = dist;
  /* initialize vertex size vector */
  ss = gsl_vector_alloc(numPara);
  /* set all step sizes */
  gsl_vector_set_all(ss, args->s);
  /* starting point */
  x = gsl_vector_alloc(numPara);
  gsl_vector_set(x, 0, args->D);
  /* initialize method and iterate */
  minex_func.f = &myF;
  minex_func.n = numPara;
  minex_func.params = (void *)NULL;
  s = gsl_multimin_fminimizer_alloc(T, numPara);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if(status)
      break;
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, args->t);
  }while(status == GSL_CONTINUE && iter < args->i);
  if(status != GSL_SUCCESS)
    printf("WARNING: Estimation of \\Delta failed: %d\n",status);
  result->de = gsl_vector_get(s->x, 0);
  result->l = s->fval;
  result->i = iter;
  result->rh = rhoFromDelta(result->pi,result->de)/dist;
  conf(args, result);
  result->rLo = rhoFromDelta(result->pi,result->dLo)/dist;
  result->rUp = rhoFromDelta(result->pi,result->dUp)/dist;
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  /* freeMlComp(); */
  return result;
}


double myF(const gsl_vector* v, void *params){
  double de;

  /* get delta */  
  de = gsl_vector_get(v, 0);
  if(de < -1 || de > 1)
    return DBL_MAX;
  else
    lik(de);

  return -likelihood;
}

void lik(double de){
  int i;
  double pi, h0, h2;
  double complementHalf;

  pi = globalPi;
  h0 = 1./(1.+pi)/(1.+pi) + de*pi/(1.+pi)/(1.+pi);
  h2 = pi*pi/(1.+pi)/(1.+pi) + de*pi/(1.+pi)/(1.+pi);
  complementHalf = (1.-h0-h2)/2.;
  likelihood = 0.;
  for(i=0;i<globalNumProfiles;i++)
    traverse(i,globalProfilePairs[i],h0,h2,complementHalf);
  
}


void traverse(int a, Node *np, double h0, double h2, double complementHalf){
  double li;
  int b;
  double *lOnes, *lTwos;

  if(np != NULL){
    traverse(a,np->left,h0,h2,complementHalf);

    lOnes = getLones();
    lTwos = getLtwos();
    b = np->key;
    li = h0*lOnes[a]*lOnes[b]
       + h2*lTwos[a]*lTwos[b]
       + complementHalf*(lOnes[a]*lTwos[b]+lTwos[a]*lOnes[b]);
    if(li>0){
      /* printf("li: %f\n",li); */
      /* printf("np->n: %f\n",log(li)*np->n); */
      likelihood += log(li) * np->n;
    }else
      likelihood += log(DBL_MIN) * np->n;

    traverse(a,np->right,h0,h2,complementHalf);
  }
}

/* conFun: called for confidence interval estimation of pi */
double confFun(double x, void *params){
  Result *res;
  double l;

  res = (Result *)params;
  lik(x);
  l = likelihood + res->l + 2.;

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

double rhoFromDelta(double t, double d){
  double r;

  r = (t + pow(t,2.) - 
       d*(13. + 19.*t + 
	  6.*pow(t,2.)) + 
       sqrt(1. + t)*
       sqrt(pow(t,2.)*(1. + t) + 
	    2.*d*t*
	    (23. + 17*t + 
	     2.*pow(t,2.)) + 
	    pow(d,2.)*
	    (97. + 109.*t + 
	     32.*pow(t,2.) + 
	     4*pow(t,3.))))/
    (2.*d*(1. + t));

  return r;
}
