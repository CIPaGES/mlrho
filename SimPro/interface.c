/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 * License: GNU General Public
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
  char c;
  char *optString = "D:t:P:c:s:C:lhv";

  args = (Args *)emalloc(sizeof(Args));

  args->s = DEFAULT_S;
  args->h = 0;
  args->e = 0;
  args->v = 0;
  args->D = DEFAULT_D;
  args->E = DEFAULT_E;
  args->t = DEFAULT_T;
  args->c = DEFAULT_C;
  args->C = DEFAULT_CC;
  args->l = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 's':                           /* number of sites */
      args->s = atoi(optarg);
      break;
    case 't':                           /* theta */
      args->t = atof(optarg);
      break;
    case 'E':                           /* epsilon */
      args->E = atof(optarg);
      break;
    case 'D':                           /* delta */
      args->D = atof(optarg);
      break;
    case 'c':                           /* coverage */
      args->c = atoi(optarg);
      break;
    case 'C':                           /* number of contigs */
      args->C = atoi(optarg);
      break;
    case 'l':                           /* compute likelihood */
      args->l = 1;
      break;
    case 'v':                           /* print program information */
      args->v = 1;
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
  return args;
}


void printUsage(char *version){
  printf("simPro version %s written by Michael Lynch, adaped by Bernhard Haubold\n", version);
  printf("Usage: simPro [options]\n");
  printf("Simulate data for maximum likelihood estimation of linkage disequilibrium from resequencing data\n");
  printf("Example: simPro > test.dat\n");
  printf("options:\n");
  printf("\t[-t <NUM> mutation rate, theta; default: NUM=%f]\n",DEFAULT_T);
  printf("\t[-E <NUM> error probability, epsilon; default: NUM=%f]\n",DEFAULT_E);
  printf("\t[-D <NUM> linkage disequilibrium for homozygosity, delta; default: NUM=%f]\n",DEFAULT_D);
  printf("\t[-c <NUM> coverage; default: %d]\n",DEFAULT_C);
  printf("\t[-C <NUM> number of contigs; default: %d]\n",DEFAULT_CC);
  printf("\t[-s <NUM> pairs of sites per contig; default: %d]\n",DEFAULT_S);
  printf("\t[-v print information about program and exit]\n");			     
  printf("\t[-h print this help message and exit]\n");
  exit(0);
}

void printSplash(char *version){
  printf("%s %s\n",progname(),version);
  printf("Written by Michael Lynch, adapted by Bernhard Haubold.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}
