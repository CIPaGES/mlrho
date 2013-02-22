/***** mlComp.h ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Feb 19 08:49:00 2009
 * License: GNU General Public
 **************************************************/
#ifndef MLCOMP
#define MLCOMP

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include "profile.h"
#include "profileTree.h"

#define MAX_ITER 1000

typedef struct result{
  double pi;   
  double ee;
  double de;   /* delta */
  double rh;   /* rho */
  double rhoFromDelta;
  double l;    /* likelihood */		 
  double pLo;  /* lower bound of pi */
  double pUp;  /* upper bound of pi */
  double eLo;  /* lower bound of epsilon */
  double eUp;  /* upper bound of epsilon */
  double dLo;  /* lower bound of delta */
  double dUp;  /* upper bound of delta */
  double rLo;  /* lower bound of rho */
  double rUp;  /* upper bound of rho */
  int i;       /* number of iterations */
  char type;   /* parameter type */
}Result;

Result *estimatePi(Profile *profiles, int numProfiles, Args *args, Result *result);
Result *estimateDelta(Node **profilePairs, int numProfiles, Args *args, Result *result, int dist);
double piComp_getNumPos(Profile *profiles, int numProfiles);
double deltaComp_getNumPos();
/* void estimateDelta(Node *r, Args *args, Result *res, int np); */
/* int estimateRho(Node *r, Args *args, Result *res, int np); */
inline double lOne(int cov, int *profile, double ee);
inline double lTwo(int cov, int *profile, double ee);
void iniMlComp(Profile *profiles, int numProfile);
void setPi(double pi);
void setEpsilon(double ee);
void rhoSetPi(double pi);
void rhoSetEpsilon(double ee);
void freeMlComp();
void setFreqNuc(double *f);
void setMaxCov(int m);
int *getCoverages();
double *getLones();
double *getLtwos();

inline double lOneDelta(int cov, int *profile, double ee);
void writeLik(char *baseName, Result *result);
Result *newResult();

#endif
