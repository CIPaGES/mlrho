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
  char *optString = "D:H:P:c:s:i:lhd:";

  args = (Args *)emalloc(sizeof(Args));

  args->s = DEFAULT_S;
  args->h = 0;
  args->e = 0;
  args->p = 0;
  args->D = DEFAULT_D;
  args->P = DEFAULT_P;
  args->H = DEFAULT_H;
  args->c = DEFAULT_C;
  args->i = DEFAULT_I;
  args->l = 0;
  args->d = NULL;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'd':                           /* data file */
      args->d = optarg;
      break;
    case 's':                           /* number of sites */
      args->s = atoi(optarg);
      break;
    case 'H':                           /* heterozygosity */
      args->H = atof(optarg);
      break;
    case 'P':                           /* phi */
      args->P = atof(optarg);
      break;
    case 'D':                           /* delta */
      args->D = atof(optarg);
      break;
    case 'c':                           /* coverage */
      args->c = atoi(optarg);
      break;
    case 'i':                           /* number of iterations */
      args->i = atoi(optarg);
      break;
    case 'l':                           /* compute likelihood */
      args->l = 1;
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
  return args;
}


void printUsage(char *version){
  printf("testMlDiseq version %s written by Michael Lynch, adaped by Bernhard Haubold\n", version);
  printf("purpose: maximum likelihood estimation of linkage disequilibrium from resequencing data\n");
  printf("usage: testMlDiseq [options]\n");
  printf("options:\n");
  printf("\t[-H <NUM> heterozygosity; default: NUM=%f]\n",DEFAULT_H);
  printf("\t[-P <NUM> error probability, phi; default: NUM=%f]\n",DEFAULT_P);
  printf("\t[-D <NUM> linkage disequilibrium for homozygosity, delta; default: NUM=%f]\n",DEFAULT_D);
  printf("\t[-c <NUM> coverage; default: %d]\n",DEFAULT_C);
  printf("\t[-s <NUM> number of sites; default: %d]\n",DEFAULT_S);
  printf("\t[-i <NUM> number of iterations; default: %d]\n",DEFAULT_I);
  printf("\t[-d <FILE> print data to FILE; default: don't print data]\n");
  printf("\t[-l compute likelihood]\n");
  printf("\t[-p print information about program and exit]\n");			     
  printf("\t[-h print this help message and exit]\n");
  exit(0);
}
