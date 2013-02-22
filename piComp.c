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
#include <string.h>
#include <sys/stat.h>
#include "interface.h"
#include "profile.h"
#include "mlComp.h"
#include "eprintf.h"

double *lOnes = NULL;
double *lTwos = NULL;
double numPos;

double likP(Profile *profiles, int numProfiles, double pi, double ee);
double myP(const gsl_vector *v, void *params);
double myPconf(double x, void *params);
double myEconf(double x, void *params);
void confP(Args *args, Result *result);
void confE(Args *args, Result *result);
void compSiteLik(Profile *profiles, int numProfiles, double ee);
FILE *openLikFile(char *baseName);
Result *readLik(FILE *fp, int numProfiles, Result *result);

/* estimatePi: estimate pi and epsilon using the Nelder-Mead
 * Simplex algorithm; code adapted from Galassi, M., Davies, 
 * J., Theiler, J., Gough, B., Jungman, G., Booth, M., 
 * Rossi, F. (2005). GNU Scientific Library Reference Manual. 
 * Edition 1.6, for GSL Version 1.6, 17 March 2005, p 472f.
 */
Result *estimatePi(Profile *profiles, int numProfiles, Args *args, Result *result){
  size_t np;
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter;
  int status;
  double size;
  FILE *fp;

  if((fp = openLikFile(args->n)) != NULL && !args->I){
    result = readLik(fp,numProfiles,result);
    fclose(fp);
    return result;
  }

  np = 2;
  T = gsl_multimin_fminimizer_nmsimplex;
  s = NULL;
  iter = 0;

  /*set up likelihood computation */
  iniMlComp(profiles, numProfiles);
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
  compSiteLik(getProfiles(),getNumProfiles(),result->ee);
  return result;
}
/* myF: function called during the minimization procedure */
double myP(const gsl_vector* v, void *params){
  double pi, ee;
  double likelihood;

  pi = gsl_vector_get(v, 0);
  ee = gsl_vector_get(v, 1);
  likelihood = 0;
  if(pi < 0 || ee < 0 || pi > 1 || ee > 1){
    return DBL_MAX;
  }
  likelihood = likP(getProfiles(), getNumProfiles(), pi, ee);
  return -likelihood;
}

double likP(Profile *profiles, int numProfiles, double pi, double ee){
  int i;
  int *coverages;
  double l, likelihood, lOneLoc, lTwoLoc;

  coverages = getCoverages();
  likelihood = 0.;
  for(i=0;i<numProfiles;i++){
    lOneLoc = lOne(coverages[i],profiles[i].profile,ee);
    lTwoLoc = lTwo(coverages[i],profiles[i].profile,ee);
    l = lOneLoc * (1.0 - pi) + lTwoLoc * pi;
    if(l>0)
      likelihood += log(l) * profiles[i].n;
  }
  return likelihood;
}

void compSiteLik(Profile *profiles, int numProfiles, double ee){
  int i;
  int *coverages;

  lOnes = (double *)emalloc(numProfiles*sizeof(double));
  lTwos = (double *)emalloc(numProfiles*sizeof(double));
  coverages = getCoverages();
  for(i=0;i<numProfiles;i++){
    lOnes[i] = lOne(coverages[i],profiles[i].profile,ee);
    lTwos[i] = lTwo(coverages[i],profiles[i].profile,ee);
  }
}

double piComp_getNumPos(Profile *profiles, int numProfiles){
  int i;
  
  if(profiles){
    numPos = 0;
    for(i=0;i<numProfiles;i++)
      numPos += profiles[i].n;
  }
  return numPos;
}

/* myPconf: called for confidence interval estimation of pi */
double myPconf(double x, void *params){
  Result *res;
  double likelihood;

  res = (Result *)params;

  likelihood = likP(getProfiles(), getNumProfiles(), x, res->ee);

  return likelihood + res->l + 2;
}

/* myEconf: called for confidence interval estimation of epsilon */
double myEconf(double x, void *params){
  Result *res;
  double l;
  double likelihood;

  res = (Result *)params;
  likelihood = likP(getProfiles(), getNumProfiles(), res->pi, x);

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

double *getLones(){
  return lOnes;
}

double *getLtwos(){
  return lTwos;
}


void writeLik(char *baseName, Result *result){
  char *fileName;
  int n;
  double np;
  FILE *fp;

  fileName = (char *)emalloc(256*sizeof(char));
  fileName = strcpy(fileName,baseName);
  fileName = strcat(fileName,".lik");

  fp = efopen(fileName,"wb");
  n = fwrite("lik",sizeof(char),3,fp);
  assert(n == 3);
  n = fwrite(result,sizeof(Result),1,fp);
  assert(n == 1);
  np = getNumProfiles();
  n = fwrite(&np,sizeof(double),1,fp);
  assert(n == 1);
  n = fwrite(getLones(),sizeof(double),getNumProfiles(),fp);
  assert(n == getNumProfiles());
  n = fwrite(getLtwos(),sizeof(double),getNumProfiles(),fp);
  assert(n == getNumProfiles());
  printf("#Likelihoods written to %s\n",fileName);
  free(fileName);
  fclose(fp);
}

FILE *openLikFile(char *baseName){
  char *fileName;
  FILE *fp;
  struct stat stbuf;
  int n;
  char *tag;

  fileName = (char *)emalloc(256*sizeof(char));
  tag = (char *)emalloc(4*sizeof(char));
  fileName = strcpy(fileName,baseName);
  fileName = strcat(fileName,".lik");
  if(stat(fileName,&stbuf) != -1){
    fp = efopen(fileName,"rb");
    n = fread(tag,sizeof(char),3,fp);
    tag[3] = '\0';
    assert(n == 3);
    if(strcmp(tag,"lik") != 0){
      assert(0);
    }
  }else
    fp = NULL;
  free(tag);
  free(fileName);
  return fp;
}

Result *readLik(FILE *fp, int numProfiles, Result *result){
  int n, i, sod;
  double np;

  assert(lOnes == NULL);
  assert(lTwos == NULL);
  fseek(fp,3,SEEK_SET);
  n = fread(result,sizeof(Result),1,fp);
  assert(n == 1);
  sod = sizeof(double);
  n = fread(&np,sod,1,fp);
  assert(n == 1 && np == numProfiles);
  lOnes = (double *)emalloc(numProfiles*sod);
  lTwos = (double *)emalloc(numProfiles*sod);
  for(i=0;i<numProfiles;i++){
    n = fread(&lOnes[i],sod,1,fp);
    assert(n == 1);
  }
  for(i=0;i<numProfiles;i++){
    n = fread(&lTwos[i],sod,1,fp);
    assert(n == 1);
  }
  return result;
}
