/***** interface.h ********************************
 * Description: Header file for user interface.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:07:28 2004.
 * Licence: GNU General Public
 ***************************************************/ 
#ifndef INTERFACE
#define INTERFACE

#define INI_PI 1e-3
#define INI_EPSILON 1e-3
#define INI_DELTA 1e-3
#define INI_RHO 1.
#define THRESHOLD 1e-8
#define STEP_SIZE 1e-4
#define MAX_IT 1000
#define DEFAULT_S 1
#define DEFAULT_N "profileDb"

/* define argument container */
typedef struct args{
  double P; /* initial pi */
  double E; /* initial epsilon */
  double D; /* initial delta */
  double R; /* initial rho */
  double t; /* threshold of simplex size */
  double s; /* step size in ML analysis */
  int S;    /* step size in LD analysis */
  int d;    /* distance for H0 and H2 computation */
  int i;    /* maximum number of iterations */
  int m;    /* minimum distance in LD analysis */
  int M;    /* maximum distance in LD analysis */
  char l;   /* compute delta */
  char I;   /* print likelihood values */
  char L;   /* lump the number of distance classes indicated by "step"? */
  char r;   /* print profiles and exit */
  char p;   /* print program information */
  char T;   /* test mode */
  char h;   /* help message? */
  char e;   /* error message? */
  char *n;  /* name of database */
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);

#endif
