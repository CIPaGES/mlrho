/***** mlComp.c ***********************************
 * Description: ML computations.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Feb 19 09:56:43 2009
 **************************************************/
#include <stdio.h>
#include <assert.h>
#include "eprintf.h"
#include "interface.h"
#include "profileTree.h"
#include "mlComp.h"

double *freqNuc, *logFreqNuc;
double totalNuc;
double S, logS;
Node *root;
double likelihood;
double lo2;
double *lngamma;
gsl_sf_result *result;
int maxCov;

void countNuc(Node *node, int d);
void compS();
inline double powInt(double x, int y);

void iniMlComp(Node *node, int d){
  int i;

  freqNuc = (double *)emalloc(4*sizeof(double));
  logFreqNuc = (double *)emalloc(4*sizeof(double));
  result = (gsl_sf_result *)emalloc(sizeof(gsl_sf_result));
  root = node;
  for(i=0;i<4;i++)
    freqNuc[i] = 0.0;
  totalNuc = 0;
  maxCov = 0;
  countNuc(root, d);
  lngamma = (double *)emalloc((maxCov+1)*sizeof(double));
  for(i=1;i<=maxCov;i++)
    lngamma[i] = gsl_sf_lngamma(i);
  for(i=0;i<4;i++){
    freqNuc[i] /= totalNuc;
    logFreqNuc[i] = log(freqNuc[i]);
  }
  compS();
  lo2 = log(2);
}

/* lOne: equation (4a) of Lynch (2008) as revised by Stephen Bates on June 24, 2012 */
inline double lOne(int cov, int *profile, double ee){
  int i;
  double s, compEe, eeThird, g;

  s = 0.0;
  compEe = 1.0 - ee;
  eeThird = ee / 3.0;
  g = 0.0;
  for(i=0;i<4;i++){
    s += freqNuc[i] * powInt(compEe, profile[i]) * powInt(eeThird,cov-profile[i]);
    g += lngamma[profile[i]+1];
  }
  g = exp(lngamma[cov+1]-g);
  s *= g;

  return s;
}

/* lTwo: equation (4b) of Lynch (2008) as revised by Stephen Bates on June 24, 2012 */
inline double lTwo(int cov, int *profile, double ee){
  int i, j;
  double s, g, x, eeThird;

  /* compute multinomial coefficient */
  g = 0.0;
  for(i=0;i<4;i++)
    g += lngamma[profile[i]+1];
  g = exp(lngamma[cov+1]-g);
  /* compute likelihood */
  s = 0.0;
  eeThird = ee/3.0;
  x = (1.0-2.0*eeThird)/2.0;
  for(i=0;i<4;i++)
    for(j=i+1;j<4;j++){
      s += freqNuc[i]*freqNuc[j]/S * powInt(x,profile[i]+profile[j]) * powInt(eeThird,cov-profile[i]-profile[j]);
    }
  s *= g;

  return s;
}

inline double powInt(double x, int y){
  int i;
  double p;

  p = 1.;
  for(i=0;i<y;i++)
    p *= x;

  return p;
}

/* compS: compute global variable S */
void compS(){
  int i;

  S = 0.0;
  for(i=0;i<4;i++)
    S += freqNuc[i]*freqNuc[i];
  S = 1.0 - S;
  logS = log(S);
}

/* countNuc: traverse profile tree and count individual 
 *   nucleotides in freqNuc, total nucleotides in
 *   totalNuc, and determine the maximum coverage
 */
void countNuc(Node *node, int d){
  int i;

  if(node != NULL){
      countNuc(node->left, d);

      for(i=0;i<4;i++){
	freqNuc[i] += node->profile1[i]*node->n;
	if(d)
	  freqNuc[i] += node->profile2[i]*node->n;
      }
      if(maxCov < node->c1)
	maxCov = node->c1;
      if(maxCov < node->c2)
	maxCov = node->c2;
      totalNuc += (node->c1+node->c2)*node->n;

      countNuc(node->right, d);
  }
}

/* setFreqNuc: set nucleotide frequencies for testing purposes */
void setFreqNuc(double *f)
{
  freqNuc = f;
}

void freeMlComp()
{
  /* int i; */

  free(freqNuc);
  free(logFreqNuc);
  free(result);
  free(lngamma);
}
