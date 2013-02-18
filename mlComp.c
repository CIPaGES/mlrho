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
gsl_sf_result *result;

void countNuc(Node *node, int d);
double logBinProb(int n, int k, double p);
void compS();


void iniMlComp(Node *node, int d){
  int i;

  freqNuc = (double *)emalloc(4*sizeof(double));
  logFreqNuc = (double *)emalloc(4*sizeof(double));
  result = (gsl_sf_result *)emalloc(sizeof(gsl_sf_result));

  root = node;
  for(i=0;i<4;i++)
    freqNuc[i] = 0.0;
  totalNuc = 0;
  countNuc(root, d);
  for(i=0;i<4;i++){
    freqNuc[i] /= totalNuc;
    logFreqNuc[i] = log(freqNuc[i]);
  }
  compS();
  lo2 = log(2);
}

/* lOne: equation (4a) of Lynch (2008) */
double lOne(double cov, int *profile, double ee){
  int i;
  double s;

  s = 0.0;
  for(i=0;i<4;i++)
    s += exp(logFreqNuc[i] + logBinProb(cov, cov-profile[i], ee));

  return s;
}

/* lTwo: equation (4b) of Lynch (2008) */
double lTwo(double cov, int *profile, double ee){
  int i, j;
  double s;

  s = 0.0;
  for(i=0;i<4;i++)
    for(j=i+1;j<4;j++){

      s += exp(lo2 + logFreqNuc[i] + logFreqNuc[j]    \
	+ logBinProb(cov,cov-profile[i]-profile[j],2.0*ee/3.0) \
	       + logBinProb(profile[i]+profile[j],profile[i],0.5) - logS);

    }
  
  return s;
}

double logBinProb(int n, int k, double p){
  int status;

  status = gsl_sf_lnchoose_e(n,k,result);
  if(status != GSL_SUCCESS)
    assert("ERROR[logBinPro]: Error in binomial coefficient.\n");

  return result->val + log(pow(p,k)) + log(pow(1.0-p,n-k));

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

      for(i=0;i<4;i++)
	freqNuc[i] += (node->profile[i] + node->profile[i+4])*node->n;
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
  free(freqNuc);
  free(logFreqNuc);
  free(result);
}
