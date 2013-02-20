/***** mlRho.c ************************************
 * Description: Maximum-likelihood estimation of 
 *   mutation and recombination rates from re-
 *   sequencing data.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Feb 18 16:35:16 2009.
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include "eprintf.h"
#include "interface.h"
#include "tab.h"
#include "profileTree.h"
#include "mlComp.h"

void runAnalysis(Args *args);

int main(int argc, char *argv[]){
  Args *args;
  char *version;

  version = "1.21";
  setprogname2("mlRho");
  args = getArgs(argc, argv);
  if(args->p)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  runAnalysis(args);
  free(args);
  freeProfileTree();
  free(progname());
  return 0;
}

void runAnalysis(Args *args){
  Node *pairwTree, *indivTree;
  Result *r;
  int i, status;
  char *headerPi, *headerDelta, *headerRho, *outStrPi; 
  char *outStrDelta1, *outStrDelta2;
  char *outStrRho1, *outStrRho2, *outStrRho3, *outStrRho4;	
  double numPos;
  FILE *fp;

  headerPi = "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\n";
  headerDelta = "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\t\tdelta\t\t\t\trho=f(delta)\n";
  headerRho =  "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\t\trho\n";
  outStrPi = "%d\t%.0f\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\n";
  outStrDelta1 = "%d\t%.0f\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\n";
  outStrDelta2 = "%d\t%.0f\t\t\t\t\t\t\t\t\t%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\n";
  outStrRho1 = "%d\t%.0f\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\t%8.2e<%8.2e<%8.2e\n";
  outStrRho2 = "%d\t%.0f\t\t\t\t\t\t\t\t\t%8.2e\t%8.2e<%8.2e<%8.2e\n";
  outStrRho3 = "%d\t%.0f\t\t\t\t\t\t\t\t\t%8.2e\t%s\n";
  outStrRho4 = "%d\t%0f\t%s\t%s\t%s\t%s\n";
  r = (Result *)emalloc(sizeof(Result));
  /* heterozygosity analysis */
  indivTree = getSummarizedProfiles(args);
  numPos = getNumPos();
  if(args->M == 0 && numPos)
    printf("%s", headerPi);
  else{
    if(args->l && numPos > 1)
      printf("%s", headerDelta);
    else if(!args->l && numPos > 1)
      printf("%s", headerRho);
  }
  if(numPos){
    estimatePi(indivTree,args,r);
    printf(outStrPi,0,numPos,r->pLo,r->pi,r->pUp,r->eLo,r->ee,r->eUp,r->l);
  }
  fflush(NULL);
  /* linkage analysis */
  if(args->T)
    setTestMode();
  /* make pi & epsilon available for one-param version of delta and rho estimation? */
  if(!args->f){
    setPi(r->pi);      
    setEpsilon(r->ee);
    rhoSetPi(r->pi);
    rhoSetEpsilon(r->ee);
  }
  if(getNumPos() > 1 && args->M - args->m > 0){
    fp = iniLinkAna(args);
    for(i=args->m;i<=args->M;i+=args->S){
      pairwTree = getProfileTree(fp,args, i);
      if(args->f){
	if(args->l){
	  estimateDelta(pairwTree, args, r, 3);
	  printf(outStrDelta1,i,getNumPos(),r->pLo,r->pi,r->pUp,r->eLo,r->ee,r->eUp,r->l, \
		 r->dLo,r->de,r->dUp,r->rhoFromDelta/(double)i);
	}else{
	  status = estimateRho(pairwTree, args, r, 3);
	  if(!status){
	    r->rh  /= (double)i;
	    r->rLo /= (double)i;
	    if(r->rUp != 1.0)
	      r->rUp /= (double)i;
	    printf(outStrRho1,i,getNumPos(),r->pLo,r->pi,r->pUp,r->eLo,r->ee,r->eUp,r->l, \
		   r->rLo,r->rh,r->rUp);
	  }else
	    printf(outStrRho4,i,getNumPos(),"n/a","n/a","n/a","n/a");
	}
      }else{
	if(args->l){
	  estimateDelta(pairwTree, args, r, 1);
	  printf(outStrDelta2,i,getNumPos(),r->l,r->dLo,r->de,r->dUp,r->rhoFromDelta/(double)i);
	}else{
	  status = estimateRho(pairwTree, args, r, 1);
	  if(status == GSL_SUCCESS){
	    r->rh  /= (double)i;
	    r->rLo /= (double)i;
	    if(r->rUp != 1.0)
	      r->rUp /= (double)i;
	    printf(outStrRho2,i,getNumPos(),r->l,r->rLo,r->rh,r->rUp);
	  }else{
	    printf(outStrRho3,i,getNumPos(),r->l,"n/a");
	  }
	}
      }
      fflush(NULL);
      freeTree(pairwTree);
    }
  }
  freeTree(indivTree);
  free(r);
}

