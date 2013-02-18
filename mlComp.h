/***** mlComp.h ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Feb 19 08:49:00 2009
 * License: GNU General Public
 **************************************************/
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

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

void estimatePi(Node *r, Args *args, Result *res);
void estimateDelta(Node *r, Args *args, Result *res, int np);
int estimateRho(Node *r, Args *args, Result *res, int np);
double lOne(double cov, int *profile, double ee);
double lTwo(double cov, int *profile, double ee);
void iniMlComp(Node *node, int d);
void setPi(double pi);
void setEpsilon(double ee);
void rhoSetPi(double pi);
void rhoSetEpsilon(double ee);
void freeMlComp();
void setFreqNuc(double *f);
void setMaxCov(int m);
void compS();
