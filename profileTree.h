/***** profileTree.h ******************************
 * Description: Header file for profileTree.c.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Feb 18 17:06:39 2009.
 * Distributed under the GNU General Public License
 **************************************************/
#ifndef PROFILETREE
#define PROFILETREE
#include <stdio.h>

#define TMPFILE "tempFile.pro"

typedef struct node{  /* the tree node: */
  char *key;          /* sort key */
  int profile1[4];    /* first profile */
  int profile2[4];    /* second profile */
  int c1;             /* coverage of first profile */
  int c2;             /* coverage of second profile */
  int n;              /* number of occurrences */
  int id;             /* id of node */
  struct node *left;  /* left child */
  struct node *right; /* right child */
}Node;

typedef struct profile{
  int pos;
  int nodeIndex;
}Profile;

typedef struct contigDescr{
  int n;                 /* number of contigs */
  int *len;              /* contig lengths */
  Profile *profileBuf;   /* buffer of profiles */
}ContigDescr;

Node *getProfileTree(int fd, Args *args, int d);
void printTree(FILE *fp, Node *node);
void freeTree(Node *n);
void setTestMode();
double getNumPos();
double getNumPosBam();
Node *newNode(char *key, int count, Node *n1, Node *n2);
Node *addTree(Node *node, char *key, int count, Node *n1, Node *n2);
int coverage(char *key, int d);
ContigDescr *getContigDescr();
void freeProfileTree();
#endif
