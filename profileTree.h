/***** profileTree.h ******************************
 * Description: Header file for profileTree.c.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Feb 18 17:06:39 2009.
 * Distributed under the GNU General Public License
 **************************************************/
#ifndef PROFILETREE
#define PROFILETREE
#include <stdio.h>

typedef struct node{  /* the tree node: */
  char *key;          /* points to the text */
  int *profile1;      /* first profile */
  int *profile2;      /* second profile */
  int c1;             /* coverage 1 */
  int c2;             /* coverage 2 */
  int n;              /* number of occurrences */
  struct node *left;  /* left child */
  struct node *right; /* right child */
}Node;

Node *getProfileTree(int fd, Args *args, int d);
void printTree(FILE *fp, Node *node);
void freeTree(Node *n);
void setTestMode();
double getNumPos();
double getNumPosBam();
Node *newNode(char *key, int count, int pair);
Node *addTree(Node *node, char *key, int count, int pair);
int coverage(char *key, int d);

#endif
