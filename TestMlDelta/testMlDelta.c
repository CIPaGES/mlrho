/* 
   Program to:
   1) generate site-specific sequence data for random genomic sequence, with arbitrary heterozygosity
   and read-error rates;
   2) do so with pairs of sites assumed to have correlated homozygosity/heterozygosity;
   3) use maximum likelihood to estimate average genome-wide nucleotide diversity and across-locus correlation.

   NOTE: This program can also be used to get the statistics for single-site analyses by setting the initial
   disequilibrium coefficient = 0.0, and not interating over it.
   The number of observed sites is then equal to twice the number of pairs entered (2 * nsites).

   Data are for a single diploid individual.
   Coverage is assumed to be the same at all sites.
   Assumes nucleotide frequencies are known without error.
*/
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <time.h>
#include "ranNum.h"
#include "eprintf.h"
#include "interface.h"

/* point to the output file */
FILE *stream;
int main(int argc, char *argv[]){
  int iters;					/* number of individuals to simulate */
  double pnuc[5];				/* nucleotide frequencies */
  int ig, jg, kg, mg, mhg, rg, sg, meg, dhg, lg;/* indicators for loops */
  int truenuc, truenuc1, truenuc2;		/* true nucleotide at a homozygous locus */
  double testnuc;				/* random number draw fro truenuc */
  int errornum;					/* number of sequence errors at the site */
  int errornuc;					/* nucleotide at error site */
  int numn[5], numn2[5];			/* number of observed nucleotides of the four types at a site */
  int numfirst, numsecond;			/* number of sequences from first and second alleles in heterozygotes */
  double facto = 0;      			/* term for Poissson */
  double poiterm[11][51], poiterm2[11][51];
  double test1, test2, test3;
  double maxll;					/* maximum log likelihood for best estimates */
  double lhood;
  double llhood[21][51][21];			/* likelihood, log likelihood */
  double besthet=0;                             /* ML estimate */
  double besterror=0;                           /* ML estimate */
  double bestdis=0;		                /* ML estimate */
  double ll1a, ll1b, ll2a, ll2b;		/* terms for the site-specific likelihoods */
  double ll1ab, ll2ab, ll3ab;
  double front1[51][51], front2[51][51], front3[51][51];
  double tothet;
  int nalleles;
  double coef = 0;
  double binopar[11][11];
  double sumhet, sumerr, sumdis, sumsqhet, sumsqerr, sumsqdis, totests;
  double meanhet, meanerr, meandis, msqhet, msqerr, msqdis, varhet, varerr, vardis, sdhet, sderr, sddis;
  double errfit[51], disfit[51], hetfit[51];
  int nhetiters, ndisiters, nerriters;	/* number of spaces on the grid search for hetero, delta, and phi */
  double hetero, phi, delta;
  int nsites, niters, coverage;
  Args *args;
  char *version;
  FILE *fp;
  int i;
  int ij = 1802;
  int kl = 9373;
  int pos = 1;

  ij = clock() % 31329; 
  rmarin(ij,kl);

  stream = stdout;

  args = getArgs(argc, argv);
  version = "0.7";
  if(args->h || args->e)
    printUsage(version);
  if(args->d)
    fp = efopen(args->d,"w");
  else
    fp = NULL;
  hetero = args->H;
  phi = args->P;
  delta = args->D;
  coverage = args->c;
  nsites = args->s;
  niters = args->i;

  /* SET THE INITIAL NUCEOTIDE FREQUENCIES. */
  pnuc[1] = 0.25;			/* frequency of A */
  pnuc[2] = 0.25;			/* frequency of C */
  pnuc[3]	= 0.25;			/* frequency of G */
  pnuc[4]	= 0.25;			/* frequency of T */
  test1 = pnuc[1];
  test2 = pnuc[1] + pnuc[2];
  test3 = pnuc[1] + pnuc[2] + pnuc[3];
  tothet = 1.0 - pow(pnuc[1],2.0) - pow(pnuc[2],2.0) - pow(pnuc[2],2.0)- pow(pnuc[3],2.0) ; 
  /* Initialize counters to zero. */
  sumhet = 0.0;
  sumerr = 0.0;
  sumdis = 0.0;
  sumsqhet = 0.0;
  sumsqerr = 0.0;
  sumsqdis = 0.0;
  totests = 0.0;
  /* Precalculate binomial probabilities of allelic draws. */
  for (ig = 0; ig <= coverage; ++ig)								
    for (jg = 0; jg <= ig; ++jg){
      if (jg == 0 || jg == ig)
	coef = 1.0;
      else  
	coef = coef * ((double)(ig) - (double)jg + 1.0) / (double)(jg);
      binopar[jg][ig] = coef * pow(0.5,(double)(ig)); 
    }
  for (rg = 0; rg <= coverage; ++rg)
    for (sg = 0; sg <= coverage; ++sg)
      if (sg > rg)
	binopar[sg][rg] = 0.0;
  binopar[0][1] = 0.5;
  binopar[1][1] = 0.5;
  binopar[0][0] = 0.0;
  /* SET UP THE GRID FOR PARAMETER EXPLORATION. */
  nhetiters = 15;
  ndisiters = 50;
  nerriters = 15;

  hetfit[1] = 0.0069;
  for (mhg = 1; mhg <= nhetiters; ++mhg)   /* loop over possible heterozygosities */
    hetfit[mhg+1] = hetfit[mhg] + 0.0001;
  disfit[1] = 0.06;
  for (dhg = 1; dhg <= ndisiters; ++dhg)   /* loop over possible disequilibria */
    disfit[dhg+1] = disfit[dhg] + 0.002;
  errfit[1] = 0.0014;
  for (meg = 1; meg <= nerriters; ++meg)   /* loop over read-error rates */
    errfit[meg+1] = errfit[meg] + 0.000015;

  for (mhg = 1; mhg <= (nhetiters + 1); ++mhg)			
    for (dhg = 1; dhg <= (ndisiters + 1); ++dhg) {			
      front1[mhg][dhg] = pow((1.0 - hetfit[mhg]),2.0) + (disfit[dhg] * hetfit[mhg] * (1.0 - hetfit[mhg]));
      front2[mhg][dhg] = pow(hetfit[mhg],2.0) + (disfit[dhg] * hetfit[mhg] * (1.0 - hetfit[mhg]));
      front3[mhg][dhg] = (1.0 - disfit[dhg]) * hetfit[mhg] * (1.0 - hetfit[mhg]);
    }
  /* Get binomial error probabilities for different read-error rates. */
  for (ig = 1; ig <= (nerriters+1); ++ig)
    for (jg = 0; jg <= coverage; ++jg) {
      if (jg == 0 || jg == coverage)
	facto = 1.0;
      else
	facto = facto * ((double)(coverage) - (double)(jg) + 1.0) / ((double)jg);
      poiterm[jg][ig] = facto * pow(errfit[ig],(double)(jg)) * pow((1.0 - errfit[ig]),(double)(coverage - jg));
    }
  for (ig = 1; ig <= (nerriters+1); ++ig)
    for (jg = 0; jg <= coverage; ++jg) {
      if (jg == 0 || jg == coverage)
	facto = 1.0;
      else
	facto = facto * ((double)(coverage) - (double)(jg) + 1.0) / (double)(jg);
      poiterm2[jg][ig] = facto * pow((0.6667 * errfit[ig]),(double)(jg)) * pow((1.0 - (0.6667 * errfit[ig])),(double)(coverage - jg));
    }
  /* Loop over the number of different simulations. */
  for (iters = 1; iters <= niters; ++iters) {
    fprintf(fp,">DataSet_%d\n",iters); 
    /* Zero the likelihood. */
    for (rg = 1; rg <= (nhetiters + 1); ++rg) 				
      for (sg = 1; sg <= (ndisiters + 1); ++sg) 
	for (ig = 1; ig <= (nerriters + 1); ++ig)
	  llhood[rg][sg][ig] = 0.0;
    /* Start generating the data. */
    for (ig = coverage; ig <= coverage; ++ig) {			
      for (lg = 1; lg <= nsites; ++lg) {		
	/* Generate data for the first site. */
	for (mg = 1; mg < 5; ++mg)
	  numn[mg] = 0;
	if (ranmar() > hetero) {					/* generate data for homozygotes */
	  nalleles = 1;
	  testnuc = ranmar();					
	  if (testnuc < test1)
	    truenuc = 1;
	  else if (testnuc < test2)
	    truenuc = 2;
	  else if (testnuc < test3)
	    truenuc = 3;
	  else 
	    truenuc = 4;
	  errornum = ignbin(ig,phi);			
	  numn[truenuc] = ig - errornum;
	  if (errornum > 0) {
	    for (kg = 1; kg <= errornum; ++kg) {
	      testnuc = ranmar();
	      if (testnuc < 0.3333333)
		errornuc = 1;
	      else if (testnuc < 0.6666667)
		errornuc = 2;
	      else
		errornuc = 3;
	      if (errornuc >= truenuc)
		errornuc++;
	      numn[errornuc]++;
	    }
	  }
	}else {									/* generate data for heterozygotes */
	  nalleles = 2;
	  testnuc = ranmar();					
	  if (testnuc < test1)
	    truenuc1 = 1;
	  else if (testnuc < test2)
	    truenuc1 = 2;
	  else if (testnuc < test3)
	    truenuc1 = 3;
	  else 
	    truenuc1 = 4;
	  testnuc = ranmar();
	  if (testnuc < 0.3333333)
	    truenuc2 = 1;
	  else if (testnuc < 0.6666667)
	    truenuc2 = 2;
	  else
	    truenuc2 = 3;
	  if (truenuc2 >= truenuc1)
	    truenuc2++;
	  numfirst = ignbin(ig,0.5);					
	  numsecond = ig - numfirst;
	  if (numfirst >= 1) {
	    errornum = ignbin(numfirst,phi);				
	    numn[truenuc1] = numfirst - errornum;
	    if (errornum > 0) {
	      for (kg = 1; kg <= errornum; ++kg) {
		testnuc = ranmar();
		if (testnuc < 0.3333333)
		  errornuc = 1;
		else if (testnuc < 0.6666667)
		  errornuc = 2;
		else
		  errornuc = 3;
		if (errornuc >= truenuc1)
		  errornuc++;
		numn[errornuc]++;
	      }
	    }
	  }
	  if (numsecond >= 1) {
	    errornum = ignbin(numsecond,phi);				
	    numn[truenuc2] = numsecond - errornum;
	    if (errornum > 0) {
	      for (kg = 1; kg <= errornum; ++kg) {
		testnuc = ranmar();
		if (testnuc < 0.3333333)
		  errornuc = 1;
		else if (testnuc < 0.6666667)
		  errornuc = 2;
		else
		  errornuc = 3;
		if (errornuc >= truenuc2)
		  errornuc++;
		numn[errornuc]++;
	      }
	    }
	  }
	}
	/* Generate data for the second site. */
	for (mg = 1; mg < 5; ++mg)
	  numn2[mg] = 0;
	/* Conditional on first site being homozygous. */
	if (nalleles == 1) {
	  if (ranmar() > (hetero - (hetero*delta)) ) {		/* generate data for homozygotes */
	    testnuc = ranmar();					
	    if (testnuc < test1)
	      truenuc = 1;
	    else if (testnuc < test2)
	      truenuc = 2;
	    else if (testnuc < test3)
	      truenuc = 3;
	    else 
	      truenuc = 4;
	    errornum = ignbin(ig,phi);			
	    numn2[truenuc] = ig - errornum;
	    if (errornum > 0) {
	      for (kg = 1; kg <= errornum; ++kg) {
		testnuc = ranmar();
		if (testnuc < 0.3333333)
		  errornuc = 1;
		else if (testnuc < 0.6666667)
		  errornuc = 2;
		else
		  errornuc = 3;
		if (errornuc >= truenuc)
		  errornuc++;
		numn2[errornuc]++;
	      }
	    }
	  }else {										/* generate data for heterozygotes */
	    testnuc = ranmar();					
	    if (testnuc < test1)
	      truenuc1 = 1;
	    else if (testnuc < test2)
	      truenuc1 = 2;
	    else if (testnuc < test3)
	      truenuc1 = 3;
	    else 
	      truenuc1 = 4;
	    testnuc = ranmar();
	    if (testnuc < 0.3333333)
	      truenuc2 = 1;
	    else if (testnuc < 0.6666667)
	      truenuc2 = 2;
	    else
	      truenuc2 = 3;
	    if (truenuc2 >= truenuc1)
	      truenuc2++;
	    numfirst = ignbin(ig,0.5);					
	    numsecond = ig - numfirst;
				
	    if (numfirst >= 1) {
	      errornum = ignbin(numfirst,phi);				
	      numn2[truenuc1] = numfirst - errornum;
	      if (errornum > 0) {
		for (kg = 1; kg <= errornum; ++kg) {
		  testnuc = ranmar();
		  if (testnuc < 0.3333333)
		    errornuc = 1;
		  else if (testnuc < 0.6666667)
		    errornuc = 2;
		  else
		    errornuc = 3;
		  if (errornuc >= truenuc1)
		    errornuc++;
		  numn2[errornuc]++;
		}
	      }
	    }

	    if (numsecond >= 1) {
	      errornum = ignbin(numsecond,phi);				
	      numn2[truenuc2] = numsecond - errornum;
	      if (errornum > 0) {
		for (kg = 1; kg <= errornum; ++kg) {
		  testnuc = ranmar();
		  if (testnuc < 0.3333333)
		    errornuc = 1;
		  else if (testnuc < 0.6666667)
		    errornuc = 2;
		  else
		    errornuc = 3;
		  if (errornuc >= truenuc2)
		    errornuc++;
		  numn2[errornuc]++;
		}
	      }
	    }
	  }
	}else { 	/* Conditional on first site being heterozygous. */
	  if (ranmar() > (hetero + (delta * (1.0 - hetero)))) {  /* generate data for homozygotes */
	    testnuc = ranmar();					
	    if (testnuc < test1)
	      truenuc = 1;
	    else if (testnuc < test2)
	      truenuc = 2;
	    else if (testnuc < test3)
	      truenuc = 3;
	    else 
	      truenuc = 4;
	    errornum = ignbin(ig,phi);			
	    numn2[truenuc] = ig - errornum;
	    if (errornum > 0) {
	      for (kg = 1; kg <= errornum; ++kg) {
		testnuc = ranmar();
		if (testnuc < 0.3333333)
		  errornuc = 1;
		else if (testnuc < 0.6666667)
		  errornuc = 2;
		else
		  errornuc = 3;
		if (errornuc >= truenuc)
		  errornuc++;
		numn2[errornuc]++;
	      }
	    }
	  } else{  /* generate data for heterozygotes */
	    testnuc = ranmar();					
	    if (testnuc < test1)
	      truenuc1 = 1;
	    else if (testnuc < test2)
	      truenuc1 = 2;
	    else if (testnuc < test3)
	      truenuc1 = 3;
	    else 
	      truenuc1 = 4;
	    testnuc = ranmar();
	    if (testnuc < 0.3333333)
	      truenuc2 = 1;
	    else if (testnuc < 0.6666667)
	      truenuc2 = 2;
	    else
	      truenuc2 = 3;
	    if (truenuc2 >= truenuc1)
	      truenuc2++;
	    numfirst = ignbin(ig,0.5);					
	    numsecond = ig - numfirst;
	    if (numfirst >= 1) {
	      errornum = ignbin(numfirst,phi);				
	      numn2[truenuc1] = numfirst - errornum;
	      if (errornum > 0) {
		for (kg = 1; kg <= errornum; ++kg) {
		  testnuc = ranmar();
		  if (testnuc < 0.3333333)
		    errornuc = 1;
		  else if (testnuc < 0.6666667)
		    errornuc = 2;
		  else
		    errornuc = 3;
		  if (errornuc >= truenuc1)
		    errornuc++;
		  numn2[errornuc]++;
		}
	      }
	    }
	    if (numsecond >= 1) {
	      errornum = ignbin(numsecond,phi);				
	      numn2[truenuc2] = numsecond - errornum;
	      if (errornum > 0) {
		for (kg = 1; kg <= errornum; ++kg) {
		  testnuc = ranmar();
		  if (testnuc < 0.3333333)
		    errornuc = 1;
		  else if (testnuc < 0.6666667)
		    errornuc = 2;
		  else
		    errornuc = 3;
		  if (errornuc >= truenuc2)
		    errornuc++;
		  numn2[errornuc]++;
		}
	      }
	    }
	  }
	}
	if(args->d){
	  fprintf(fp,"%d",pos++);
	  for(i=1;i<5;i++)
	    fprintf(fp,"\t%d",numn[i]); 
	  fprintf(fp,"\n");
	  fprintf(fp,"%d",pos++);
	  for(i=1;i<5;i++)
	    fprintf(fp,"\t%d",numn2[i]); 
	  fprintf(fp,"\n");
	}
	if(args->l){   /* likelihood computation */
	  /* Add to the likelihood values for all combinations of values along the grid. */
	  for (mhg = 1; mhg <= (nhetiters + 1); ++mhg) {   /* loop over possible heterozygosities */
	    for (dhg = 1; dhg <= (ndisiters + 1); ++dhg) { /* loop over possible disequilibria */
	      for (meg = 1; meg <= (nerriters + 1); ++meg){/* loop over read-error rates */
		ll1a = 0.0;
		ll1b = 0.0;
		ll2a = 0.0;
		ll2b = 0.0;
		for (jg = 1; jg <= 4; ++jg) {
		  ll1a = ll1a + ( pnuc[jg] * poiterm[coverage - numn[jg]][meg] );
		  ll1b = ll1b + ( pnuc[jg] * poiterm[coverage - numn2[jg]][meg] );
		  for (kg = 1; kg <= 4; ++kg) {
		    if (kg >= jg) {		
		      ll2a = ll2a + ( 2.0 * pnuc[jg] * pnuc[kg] * poiterm2[coverage - numn[jg] - numn[kg]][meg] * binopar[numn[jg]][numn[jg]+numn[kg]] / tothet);
		      ll2b = ll2b + ( 2.0 * pnuc[jg] * pnuc[kg] * poiterm2[coverage - numn2[jg] - numn2[kg]][meg] * binopar[numn2[jg]][numn2[jg]+numn2[kg]] / tothet);

		    }
		  }
		}
		ll1ab = ll1a * ll1b;
		ll2ab = ll2a * ll2b;
		ll3ab = (ll1a * ll2b) + (ll1b * ll2a);
		lhood = front1[mhg][dhg] * ll1ab;
		lhood = lhood + (front2[mhg][dhg] * ll2ab);
		lhood = lhood + (front3[mhg][dhg] * ll3ab);
		llhood[mhg][dhg][meg] = llhood[mhg][dhg][meg] + log(lhood);
	      }
	    }
	  }
	}
      }
      maxll = -10000000000.0;
      for (mhg = 1; mhg <= (nhetiters + 1); ++mhg) 			
	for (dhg = 1; dhg <= (ndisiters + 1); ++dhg) 			
	  for (meg = 1; meg <= (nerriters + 1); ++meg) 
	    if (llhood[mhg][dhg][meg] > maxll) {
	      maxll = llhood[mhg][dhg][meg];
	      besthet = hetfit[mhg];
	      bestdis = disfit[dhg];
	      besterror = errfit[meg];
	    }
      sumhet = sumhet + besthet;
      sumerr = sumerr + besterror;
      sumdis = sumdis + bestdis;
      sumsqhet = sumsqhet + pow(besthet,2.0);
      sumsqerr = sumsqerr + pow(besterror,2.0);
      sumsqdis = sumsqdis + pow(bestdis,2.0);
      totests = totests + 1.0;
    }
    if(args->l)
      printf("pi:             %8.2e             epsilon:             %8.2e             delta:             %8.2e             -log(L): %8.2e\n",besthet,besterror,bestdis,-maxll);
  }
  meanhet = sumhet / totests;
  meanerr = sumerr / totests;
  meandis = sumdis / totests;
  msqhet = sumsqhet / totests;
  msqerr = sumsqerr / totests;
  msqdis = sumsqdis /totests;
  varhet = (totests/(totests - 1.0)) * (msqhet - pow(meanhet,2.0));
  varerr = (totests/(totests - 1.0)) * (msqerr - pow(meanerr,2.0));
  vardis = (totests/(totests - 1.0)) * (msqdis - pow(meandis,2.0));

  sdhet = pow(varhet,0.5);
  sderr = pow(varerr,0.5);
  sddis = pow(vardis,0.5);

/*   printf("\n"); */
/*   printf("      mean heterozygosity = %8.7f \n",meanhet);  */
/*   printf("          mean read-error = %8.7f \n",meanerr);  */
/*   printf("mean disequilibrium coef. = %8.7f \n",meandis);  */
/*   printf("\n"); */

/*   fprintf(stream, "%8.7f , %8.7f , %d, %d ,, %8.7f , %8.7f ,, %8.7f , %8.7f ,,, %8.7f ,, %8.7f , %8.7f ,, %d\n", hetero , phi , coverage , nsites , meanhet, sdhet, meanerr , sderr, delta, meandis, sddis, niters); */
  fclose(fp);
  exit(0);
}




