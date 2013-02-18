/***** interface.h ********************************
 * Description: Header file for user interface.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:07:28 2004.
 * Licence: GNU General Public
 ***************************************************/ 

#define INI_PI 1e-3
#define INI_EPSILON 1e-3
#define INI_DELTA 1e-3
#define INI_RHO 1.
#define THRESHOLD 1e-8
#define STEP_SIZE 1e-4
#define MAX_IT 100
#define MIN_COV 4
#define DEFAULT_S 1
#define DEFAULT_B 4096

/* define argument container */
typedef struct args{
  double P; /* initial pi */
  double E; /* initial epsilon */
  double D; /* initial delta */
  double R; /* initial rho */
  double t; /* threshold of simplex size */
  double s; /* step size in ML analysis */
  int b;    /* buffer size */
  int S;    /* step size in LD analysis */
  int d;    /* distance for H0 and H2 computation */
  int i;    /* maximum number of iterations */
  int c;    /* minimum coverage */
  int m;    /* minimum distance in LD analysis */
  int M;    /* maximum distance in LD analysis */
  char l;   /* compute delta */
  char r;   /* print profiles and exit */
  char p;   /* print program information */
  char f;   /* full likelihood analysis over varying distances */
  char T;   /* test mode */
  char h;   /* help message? */
  char e;   /* error message? */
  char **inputFiles;
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);
