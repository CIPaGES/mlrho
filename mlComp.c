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

#define MAXNK 100

double *freqNuc, *logFreqNuc;
double totalNuc;
double S, logS;
Node *root;
double likelihood;
double lo2;
double **lnchoose;
gsl_sf_result *result;

void countNuc(Node *node, int d);
double logBinProb(int n, int k, double p);
void compS();


void iniMlComp(Node *node, int d){
  int i, j, status;

  freqNuc = (double *)emalloc(4*sizeof(double));
  logFreqNuc = (double *)emalloc(4*sizeof(double));
  result = (gsl_sf_result *)emalloc(sizeof(gsl_sf_result));
  lnchoose = (double **)emalloc((MAXNK+1)*sizeof(double *));
  for(i=0;i<=MAXNK;i++){
    lnchoose[i] = (double *)emalloc((i+1)*sizeof(double));
    for(j=0;j<=i;j++){
      status = gsl_sf_lnchoose_e(i,j,result);
      if(status != GSL_SUCCESS)
	assert("ERROR[iniMlComp]: Error in binomial coefficient.\n");
      lnchoose[i][j] = result->val;
    }
  }
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

/* lOne: equation (4a) of Lynch (2008) as revised by Stephen Bates on June 24, 2012 */
double lOne(double cov, int *profile, double ee){
  int i;
  double s, compEe, eeThird, g;

  s = 0.0;
  compEe = 1.0 - ee;
  eeThird = ee / 3.0;
  g = 0.0;
  for(i=0;i<4;i++){
    s += freqNuc[i] * pow(compEe, profile[i]) * pow(eeThird,cov-profile[i]);
    g += gsl_sf_lngamma(profile[i]+1);
  }
  g = exp(gsl_sf_lngamma(cov+1) - g);
  s *= g;

  return s;
}

/* lTwo: equation (4b) of Lynch (2008) as revised by Stephen Bates on June 24, 2012 */
double lTwo(double cov, int *profile, double ee){
  int i, j;
  double s, g, x, eeThird;

  /* compute multinomial coefficient */
  g = 0.0;
  for(i=0;i<4;i++)
    g += gsl_sf_lngamma(profile[i] + 1);
  g = exp(gsl_sf_lngamma(cov+1) - g);

  /* compute likelihood */
  s = 0.0;
  x = (1.0-2.0*ee/3.0)/2.0;
  eeThird = ee/3.0;
  for(i=0;i<4;i++)
    for(j=i+1;j<4;j++){
      s += freqNuc[i]*freqNuc[j]/S * pow(x,profile[i]+profile[j]) * pow(eeThird,cov-profile[i]-profile[j]);
    }
  s *= g;

  return s;
}

double logBinProb(int n, int k, double p){
  int status;

  if(n<=MAXNK && k<=MAXNK){
    result->val = lnchoose[n][k];
    status = GSL_SUCCESS;
  }else
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
  int profile[8];
  int c1, c2;

  if(node != NULL){
    c1 = 0;
    c2 = 0;
      countNuc(node->left, d);
      if(d){
	  sscanf(node->key,"%d %d %d %d %d %d %d %d",&profile[0],&profile[1],&profile[2],&profile[3], \
		 &profile[4],&profile[5],&profile[6],&profile[7]);
	  for(i=0;i<4;i++){
	    c1 += profile[i];
	    c2 += profile[i+4];
	  }
      }else{
	  sscanf(node->key,"%d %d %d %d",&profile[0],&profile[1],&profile[2],&profile[3]);
	  for(i=0;i<4;i++)
	    c1 += profile[i];
	  c2 = 0;
	  for(i=4;i<8;i++)
	      profile[i] = 0;
      }

      for(i=0;i<4;i++)
	freqNuc[i] += (profile[i] + profile[i+4])*node->n;
      totalNuc += (c1+c2)*node->n;
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
  int i;

  free(freqNuc);
  free(logFreqNuc);
  free(result);
  for(i=0;i<=MAXNK;i++)
    free(lnchoose[i]);
  free(lnchoose);
}
