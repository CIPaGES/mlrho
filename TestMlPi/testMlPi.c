/* 

   Program testMlPi.c to:

   1) generate site-specific sequence data for random genomic sequence, with arbitrary heterozygosity
   and read-error rates;
   2) use ML to estimate these parameters in retrospect.

   Data are for a single diploid individual.
   Coverage is assumed to be the same at all sites.

   NOTE THAT SOME SETTINGS NEED TO BE MADE IN THE MAIN PROGRAM WHERE HEADINGS ARE IN CAPS:
   1) THE NUCLEOTIDE USAGE FREQUENCIES.
   2) THE GRID PARAMETERS FOR THE ML SEARCH.

   Author: Michael Lynch
   Modifications: 
   Author            Date                Comment
   Bernhard Haubold  January  23, 2009  Reorganization of code
   Bernhard Haubold  February  6, 2009  Reorganization of code
   Bernhard Haubold  February 25, 2009  Included option for printing data to file

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "eprintf.h"
#include "interface.h"
#include "ranNum.h"

int main(int argc, char *argv[]){
  int ****nobs;		                    /* number of sites observed with specific counts for nucleotides A, C, G, and T */
  int iters;          			    /* number of individuals to simulate */
  double pnuc[5];                           /* nucleotide frequencies */
  int ig, jg, kg, lg, mg, rg, sg;           /* indicators for loops */
  int truenuc, truenuc1, truenuc2;	    /* true nucleotide at a homozygous locus */
  double testnuc;		            /* random number draw fro truenuc */
  int errornum;				    /* number of sequence errors at the site */
  int errornuc;				    /* nucleotide at error site */
  int numn[5];				    /* number of observed nucleotides of the four types at a site */
  int numfirst, numsecond;		    /* number of sequences from first and second alleles in heterozygotes */
  double test1, test2, test3;
  int totnuc1, totnuc2, totnuc3, totnuc4, totnucs;
  double tothet;			    /* estimated total heterozygosity */
  double estpnuc[5];			    /* estimated genome-wide nucleotide frequencies */
  double mlerror, mlhet;
  double maxll;				    /* maximum log likelihood for best estimates */
  double lhood,  llhood;		    /* likelihood, log likelihood */
  double besthet, besterror;		    /* ML estimates */
  double ll1, ll2;			    /* terms for the site-specific likelihoods */

  double hetero;                            /* average nucleotide heterozygosity */
  double phi;                               /* read-error rate */
  double coverage;                          /* average sequencing coverage per site */
  int nsites;                               /* number of informative sites */
  int niters;                               /* number of individuals iterated */

  double coef, mlstep;
  double binopar[25][25];
  double sumhet, sumerr, sumsqhet, sumsqerr, totests;
  double meanhet, meanerr, msqhet, msqerr, varhet, varerr, sdhet, sderr;
  double poiterm[100], poiterm2[100], facto;
  int jgg, kgg;
  int nobst[5];
  int ij, kl;
  Args *args;
  char *version = "0.7";
  int i,j,k,l,m;
  FILE *fp;
  int pos;

  /* placate compiler */
  jg = 0;
  coef = 0;
  besterror = 0;
  besthet = 0;
  facto = 0;

  args = getArgs(argc, argv);

  if(args->h || args->e)
    printUsage(version);

  hetero = args->H;
  phi = args->P;
  coverage = args->c;
  nsites = args->s;
  niters = args->i;

  nobs = (int ****)emalloc((args->c+1)*sizeof(int ***));
  for(i=0;i<args->c+1;i++)
    nobs[i] = (int ***)emalloc((args->c+1)*sizeof(int **));
  for(i=0;i<args->c+1;i++)
    for(j=0;j<args->c+1;j++)
      nobs[i][j] = (int **)emalloc((args->c+1)*sizeof(int *));
  for(i=0;i<args->c+1;i++)
    for(j=0;j<args->c+1;j++)
      for(k=0;k<args->c+1;k++)
	nobs[i][j][k] = (int *)emalloc((args->c+1)*sizeof(int));

  if(args->d != NULL)
    fp = efopen(args->d,"w");
  else
    fp = NULL;

  /* NOTE: The seed variables can have values between:    0 <= IJ <= 31328 */
  /*                                                      0 <= KL <= 30081 */
  ij = 1802;
  kl = 9373;
/*   ij = time((long *) 0) % 31329; */
  ij = clock() % 31329;
  /* Initialization of ranmar() */
  rmarin(ij,kl);

  /* SET THE INITIAL NUCLEOTIDE FREQUENCIES. */
  pnuc[1] = 0.25;			/* frequency of A */
  pnuc[2] = 0.25;			/* frequency of C */
  pnuc[3] = 0.25;			/* frequency of G */
  pnuc[4] = 0.25;			/* frequency of T */

  test1 = pnuc[1];
  test2 = pnuc[1] + pnuc[2];
  test3 = pnuc[1] + pnuc[2] + pnuc[3];
  /* Initialize counters to zero. */
  sumhet = 0.0;
  sumerr = 0.0;
  sumsqhet = 0.0;
  sumsqerr = 0.0;
  totests = 0.0;
  /* Binomial probabilities of allelic draws. */
  for (ig = 0; ig <= coverage; ++ig)
    for (jg = 0; jg <= ig; ++jg){
      if (jg == 0 || jg == ig)
	coef = 1.0;
      else 
	coef = coef * ((double)(ig) - (double)(jg) + 1.0) / (double)(jg);
      binopar[jg][ig] = coef * pow(0.5,(double)(ig));
    }
  for (rg = 0; rg <= coverage; ++rg)
    for (sg = 0; sg <= coverage; ++sg)
      if (sg > rg)
	binopar[sg][rg] = 0.0;
  binopar[0][1] = 0.5;
  binopar[1][1] = 0.5;
  binopar[0][0] = 0.0;

  /* Loop over the number of different simulations. */
  for (iters = 1; iters <= niters; ++iters) {
    totnuc1 = 0;
    totnuc2 = 0;
    totnuc3 = 0;
    totnuc4 = 0;
    totnucs = 0;
    /* Zero the counters for various sequence configurations. */
    for (ig = 0; ig <= coverage; ++ig) 
      for (jg = 0; jg <= coverage; ++jg) 
	for (kg = 0; kg <= coverage; ++kg) 
	  for (lg = 0; lg <= coverage; ++lg) 
	    nobs[ig][jg][kg][lg] = 0;
    /* Start generating the data. */
    for (ig = coverage; ig <= coverage; ++ig){
      for (jg = 1; jg <= nsites; ++jg){
	for (mg = 1; mg < 5; ++mg)
	  numn[mg] = 0;
	if (ranmar() > hetero) {   /* generate data for homozygotes */
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
		errornuc = errornuc + 1;
	      numn[errornuc] = numn[errornuc] + 1;
	    }
	  }
	  nobs[numn[1]][numn[2]][numn[3]][numn[4]] = nobs[numn[1]][numn[2]][numn[3]][numn[4]] + 1;
	}else {			/* generate data for heterozygotes */
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
	    truenuc2 = truenuc2 + 1;
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
		  errornuc = errornuc + 1;
		numn[errornuc] = numn[errornuc] + 1;
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
		  errornuc = errornuc + 1;
		numn[errornuc] = numn[errornuc] + 1;
	      }
	    }
	  }

	  nobs[numn[1]][numn[2]][numn[3]][numn[4]] = nobs[numn[1]][numn[2]][numn[3]][numn[4]] + 1;
	}
      }
    }
    if(fp != NULL){
      pos = 0;
      fprintf(fp,">DataSet_%d\n",iters);
      for(i=0;i<=coverage;i++)
	for(j=0;j<=coverage;j++)
	  for(k=0;k<=coverage;k++)
	    for(l=0;l<=coverage;l++)
	      for(m=0;m<nobs[i][j][k][l];m++)
		fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",++pos,i,j,k,l);
    }

    /* Estimate the nucleotide frequencies. */
    for (ig = 0; ig <= coverage; ++ig) 
      for (jg = 0; jg <= coverage; ++jg) 
	for (kg = 0; kg <= coverage; ++kg) 
	  for (lg = 0; lg <= coverage; ++lg) 
	    if (nobs[ig][jg][kg][lg] > 0) {
	      totnuc1 = totnuc1 + (nobs[ig][jg][kg][lg] * ig);
	      totnuc2 = totnuc2 + (nobs[ig][jg][kg][lg] * jg);
	      totnuc3 = totnuc3 + (nobs[ig][jg][kg][lg] * kg);
	      totnuc4 = totnuc4 + (nobs[ig][jg][kg][lg] * lg);
	    }
    totnucs = totnuc1 + totnuc2 + totnuc3 + totnuc4;
    estpnuc[1] = (double)(totnuc1)/(double)(totnucs);						
    estpnuc[2] = (double)(totnuc2)/(double)(totnucs);
    estpnuc[3] = (double)(totnuc3)/(double)(totnucs);
    estpnuc[4] = (double)(totnuc4)/(double)(totnucs);
    tothet = 1.0 - pow(estpnuc[1],2.0) - pow(estpnuc[2],2.0) - pow(estpnuc[3],2.0) - pow(estpnuc[4],2.0) ; 

    /* Obtain the ML estimates. */
    maxll = -10000000000.0;
    mlstep = 0.00001;
    for (mlhet=0;mlhet<=0.01;mlhet+=mlstep) {	/* LOOP OVER POSSIBLE HETEROZYGOSITIES. */
      for (mlerror=0;mlerror<=0.01;mlerror+=mlstep) {	/* LOOP OVER ERROR RATES TO EXAMINE */
	llhood = 0.0;
	/* get binomial error probabilities for different read-error rates to use in the likelihood function*/
	for (jg = 0; jg <= coverage; ++jg) {
	  if (jg == 0 || jg == coverage)
	    facto = 1.0;
	  else 
	    facto = facto * ((double)(coverage) - (double)(jg) + 1.0) / (double)(jg);
	  poiterm[jg] = facto * pow(mlerror,(double)(jg)) * pow((1.0 - mlerror),(double)(coverage - jg)) ;
	}

	for (jg = 0; jg <= coverage; ++jg) {
	  if (jg == 0 || jg == coverage)
	    facto = 1.0;
	  else 
	    facto = facto * ((double)(coverage) - (double)(jg) + 1.0) / (double)(jg);
	  poiterm2[jg] = facto * pow((2.0/3.0 * mlerror),(double)(jg)) * pow((1.0 - (2.0/3.0 * mlerror)),(double)(coverage - jg)) ;
	}
	/* sum over the likelihood function for each pair of pi and phi estimates */
	for (ig = 0; ig <= coverage; ++ig) {
	  for (jg = 0; jg <= coverage; ++jg) {
	    for (kg = 0; kg <= coverage; ++kg) {
	      for (lg = 0; lg <= coverage; ++lg) {
		if (nobs[ig][jg][kg][lg] > 0) {
		  ll1 = 0.0;
		  ll2 = 0.0;
		  nobst[1] = ig;
		  nobst[2] = jg;
		  nobst[3] = kg;
		  nobst[4] = lg;
		  for (jgg = 1; jgg <= 4; ++jgg) {
		    ll1 = ll1 + ( estpnuc[jgg] * poiterm[(int)(coverage - nobst[jgg])] );
		    for (kgg = 1; kgg <= 4; ++kgg) {
		      if (kgg >= jgg) {
			ll2 = ll2 + ( 2.0 * estpnuc[jgg] * estpnuc[kgg] * poiterm2[(int)(coverage - nobst[jgg] - nobst[kgg])] \
				      * binopar[nobst[jgg]][nobst[jgg] + nobst[kgg]] / tothet);
		      }}}
		  lhood = ((1.0 - mlhet) * ll1) + (mlhet * ll2);
		  llhood = llhood + ((double)(nobs[ig][jg][kg][lg]) * log(lhood));
		}
	      }}}}
	if (llhood > maxll) {
	  maxll = llhood;
	  besthet = mlhet;
	  besterror = mlerror;
	}
      }	 /* end the error loop for the likelihood surface */
    }	 /* end the heterozygosity loop for the likelihood surface */
    printf("pi:             %8.2e             epsilon:             %8.2e             -log(L): %8.2e\n", besthet, besterror,-maxll);
    sumhet = sumhet + besthet;
    sumerr = sumerr + besterror;
    sumsqhet = sumsqhet + pow(besthet,2.0);
    sumsqerr = sumsqerr + pow(besterror,2.0);
    totests = totests + 1.0;
  } /* end the loop for individuals */
  fclose(fp);

  /* calculate the means and variances of estimates over all individual runs */
  meanhet = sumhet / totests;
  meanerr = sumerr / totests;
  msqhet = sumsqhet / totests;
  msqerr = sumsqerr / totests;
  varhet = (totests/(totests - 1.0)) * (msqhet - pow(meanhet,2.0));
  varerr = (totests/(totests - 1.0)) * (msqerr - pow(meanerr,2.0));
  sdhet = pow(varhet,0.5);
  sderr = pow(varerr,0.5);

/*   printf("Het: %8.7f +/- %8.7f\n",meanhet , sdhet);  */
/*   printf("Err: %8.7f +/- %8.7f\n",meanerr , sderr);  */

/*   printf("Data summary: hetero, phi, coverage, nsites, meanhet, sdhet, meanerr, sderr, niters\n"); */
/*   printf("%7.6f, %7.6f, %d, %d, %8.7f, %8.7f, %8.7f, %8.7f, %d\n", hetero , phi , (int)coverage , nsites , meanhet, sdhet, meanerr , sderr , niters); */

  return 0;
}
