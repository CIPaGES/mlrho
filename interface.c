/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]){
  int c;
  char *optString = "P:E:D:R:t:s:i:c:rThfpM:lm:S:b:";

  args = (Args *)emalloc(sizeof(Args));
  args->P = INI_PI;
  args->E = INI_EPSILON;
  args->D = INI_DELTA;
  args->R = INI_RHO;
  args->t = THRESHOLD;
  args->c = MIN_COV;
  args->T = 0;
  args->s = STEP_SIZE;
  args->S = DEFAULT_S;
  args->i = MAX_IT;
  args->M = INT_MAX;
  args->b = DEFAULT_B;
  args->f = 0;
  args->m = 1;
  args->l = 0;
  args->r = 0;
  args->h = 0;
  args->e = 0;
  args->p = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'P':                           /* initial value of pi */
      args->P = atof(optarg);
      break;
    case 'E':                           /* initial error rate, epsilon */
      args->E = atof(optarg);
      break;
    case 'D':                           /* initial disequilibrium coefficient, delta */
      args->D = atof(optarg);
      break;
    case 'R':                           /* initial recombination parameter, Rho */
      args->R = atof(optarg);
      break;
    case 't':                           /* simplex size threshold */
      args->t = atof(optarg);
      break;
    case 's':                           /* step size in ML computation */
      args->s = atof(optarg);
      break;
    case 'S':                           /* step size in LD computation */
      args->S = atoi(optarg);
      break;
    case 'i':                           /* maximum number of iterations */
      args->i = atoi(optarg);
      break;
    case 'd':                           /* distance between pairs of profiles for H0 and H2 computation */
      args->d = atoi(optarg);
      break;
    case 'c':                           /* minimum coverage */
      args->c = atoi(optarg);
      break;
    case 'M':                           /* maximum distance investigated in disequilibrium analysis */
      args->M = atoi(optarg);
      break;
    case 'm':                           /* minimum distance investigated in disequilibrium analysis */
      args->m = atoi(optarg);
      break;
    case 'b':                           /* size of buffer for reading input data */
      args->b = atoi(optarg);
      break;
    case 'T':                           /* test mode for linkage analysis */
      args->T = 1;
      break;
    case 'l':                           /* estimate delta */
      args->l = 1;
      break;
    case 'r':                           /* print profiles */
      args->r = 1;
      break;
    case 'f':                           /* full likelihood analysis for delta computation */
      args->f = 1;
      break;
    case 'p':                           /* print program information */
      args->p = 1;
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    default:
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  if(args->numInputFiles == 0 && !args->h && !args->e && !args->p){
    printf("ERROR[mlRho]: input must come from file rather than stdin.\n");
    args->e = 1;
  }
  return args;
}


void printUsage(char *version){
  printf("mlRho version %s\n", version);
  printf("purpose: maximum likelihood estimation of theta and rho\n");
  printf("usage: mlRho [options] [inputFile(s)]\n");
  printf("options:\n");
  printf("\t[-c <NUM> minimum coverage; default: %d]\n",MIN_COV);
  printf("\t[-m <NUM> minimum distance analyzed in rho computation; default: 1]\n");
  printf("\t[-M <NUM> maximum distance analyzed in rho computation; default: all]\n");
  printf("\t[-S <NUM> step size in rho computation; default: %d]\n",DEFAULT_S);
  printf("\t[-P <NUM> initial theta value; default: %10.3e]\n",INI_PI);
  printf("\t[-E <NUM> initial epsilon value; default: %10.3e]\n",INI_EPSILON);
  printf("\t[-R <NUM> initial rho value; default: %10.3e]\n",INI_RHO);
  printf("\t[-D <NUM> initial delta value; default: %10.3e]\n",INI_DELTA);
  printf("\t[-t <NUM> simplex size threshold; default: %10.3e]\n",THRESHOLD);
  printf("\t[-s <NUM> size of first step in ML estimation; default: %10.3e]\n",STEP_SIZE);
  printf("\t[-b <NUM> size of buffer for reading input; default: %d]\n",DEFAULT_B);
  printf("\t[-l estimate delta; default: estimate rho]\n");
  printf("\t[-r print profiles and exit]\n");
  printf("\t[-T test mode for linkage analysis]\n");
  printf("\t[-f full likelihood computation over varying distances (slow);\n");
  printf("\t\tdefault: use initial estimates of \\epsilon and \\theta]\n");
  printf("\t[-p print information about program and exit]\n");			     
  printf("\t[-h print this help message and exit]\n");
  exit(0);
}

void printSplash(char *version){
  printf("********************************************************\n");
  printf("*                 mlRho, version %s                   *\n", version);
  printf("*   ML Estimation of Mutation and Recombination Rate   *\n");
  printf("* Bernhard Haubold, Peter Pfaffelhuber, Michael Lynch  *\n");
  printf("*------------------------------------------------------*\n");
  printf("* REFERENCES                                           *\n");
  printf("* 1) Lynch, M. (2008). Estimation of nucleotide        *\n");
  printf("* diversity, disequilibrium coefficients, and mutation *\n");
  printf("* rates from high-coverage genome-sequencing projects. *\n");
  printf("* Mol. Biol. Evol., 25:2409--2419.                     *\n");
  printf("* 2) Haubold, B., Pfaffelhuber, P. and Lynch, M.       *\n");
  printf("* (2009). mlRho - A program for estimating the pop-    *\n");
  printf("* ulation mutation and recombination rates from shot-  *\n");
  printf("* gun-sequenced diploid genomes. Mol. Ecol., sub-      *\n");
  printf("* mitted.                                              *\n");
  printf("*                                                      *\n");  
  printf("* CONTACT                                              *\n");
  printf("* Code maintained by Bernhard Haubold,                 *\n");
  printf("* haubold@evolbio.mpg.de                               *\n");
  printf("*                                                      *\n");
  printf("* LICENSE                                              *\n");
  printf("* This software is distributed under the GNU General   *\n");
  printf("* Public License. You should have received a copy      *\n");
  printf("* of the licence together with this software. If       *\n");
  printf("* not, see http://www.gnu.org/licenses/.               *\n");
  printf("********************************************************\n");
  exit(0);
}
