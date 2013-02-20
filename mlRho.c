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
#include "eprintf.h"
#include "interface.h"
#include "tab.h"
#include "profileTree.h"
#include "mlComp.h"

void runAnalysis(int fd, Args *args);

int main(int argc, char *argv[]){
  Args *args;
  char *version;
  int i, fd;

  version = "1.11";
  setprogname2("mlRho");
  args = getArgs(argc, argv);
  if(args->p)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  for(i=0;i<args->numInputFiles;i++){
    fd = open(args->inputFiles[i],O_RDONLY,0);
    runAnalysis(fd, args);
    close(fd);
  }
  free(args);
  free(progname());
  return 0;
}

void runAnalysis(int fd, Args *args){
  Node *root;
  Result *r;
  int i, status;
  char *headerPi, *headerDelta, *headerRho, *outStrPi; 
  char *outStrDelta1, *outStrDelta2;
  char *outStrRho1, *outStrRho2, *outStrRho3, *outStrRho4;

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
  root = getProfileTree(fd,args, 0);
  if(args->M == 0 && getNumPos())
    printf("%s", headerPi);
  else{
    if(args->l && getNumPos() > 1)
      printf("%s", headerDelta);
    else if(!args->l && getNumPos() > 1)
      printf("%s", headerRho);
  }
  if(getNumPos()){
    estimatePi(root,args,r);
    printf(outStrPi,0,getNumPos(),r->pLo,r->pi,r->pUp,r->eLo,r->ee,r->eUp,r->l);
  }
  fflush(NULL);
  freeTree(root);
  /* linkage analysis */
  if(args->T)
    setTestMode();
  /* make pi & epsilon available for one-param version of delta estimation? */
  if(!args->f){
    setPi(r->pi);      
    setEpsilon(r->ee);
    rhoSetPi(r->pi);
    rhoSetEpsilon(r->ee);
  }
  if(getNumPos() > 1){
    for(i=args->m;i<=args->M;i+=args->S){
      lseek(fd, 0L, 0);
      root = getProfileTree(fd,args, i);
      if(args->f){
	if(args->l){
	  estimateDelta(root, args, r, 3);
	  printf(outStrDelta1,i,getNumPos(),r->pLo,r->pi,r->pUp,r->eLo,r->ee,r->eUp,r->l, \
		 r->dLo,r->de,r->dUp,r->rhoFromDelta/(double)i);
	}else{
	  status = estimateRho(root, args, r, 3);
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
	  estimateDelta(root, args, r, 1);
	  printf(outStrDelta2,i,getNumPos(),r->l,r->dLo,r->de,r->dUp,r->rhoFromDelta/(double)i);
	}else{
	  status = estimateRho(root, args, r, 1);
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
      freeTree(root);
    }
  }
  free(r);
}

