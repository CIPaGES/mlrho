/***** profileTree.h ******************************
 * Description: Header file for profileTree.c.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Feb 18 17:06:39 2009.
 * Distributed under the GNU General Public License
 **************************************************/
#include <stdio.h>

typedef struct node{  /* the tree node: */
  char *key;          /* points to the text */
  int profile[8];     /* allele profile */
  int n;              /* number of occurrences */
  int c1, c2;         /* coverage */
  struct node *left;  /* left child */
  struct node *right; /* right child */
}Node;

typedef struct profileCollection{ /* collection of profile sets */
  int numSet;                     /* number of sets */
  int *setSize;                   /* set sizes */
  int ***intDat;                  /* integer representation of quartets */
  int **cov;                      /* coverage */
  int **pos;                      /* positions */
  char ***strDat;                 /* string representation of quartets */
  int maxD;                       /* maximum distance between pairs of positions */
  int np;                         /* number of positions surveyed */
  double hetPairFreq;             /* frequency of heterozygous pairs some distance apart */
  double homPairFreq;             /* frequency of homozygous pairs some distance apart */
}ProfileCollection;

Node *getProfileTree(Args *args, ProfileCollection *pc, int d);
void printTree(FILE *fp, Node *node);
int getNumFields();
ProfileCollection *readProfiles(FILE *fp);
void freeProfiles(ProfileCollection *pc);
void printProfiles(ProfileCollection *pc);
void freeTree(Node *n);
void setTestMode();
void pairFreq(ProfileCollection *pc, int d);
