/***** profileTree.c ******************************
 * Description: Tree for counting allele profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Feb 18 17:06:35 2009.
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "eprintf.h"
#include "tab.h"
#include "stringUtil.h"
#include "interface.h"
#include "profileTree.h"
#include "ld.h"

int numPos;
Node **resetProfilePairs(Node **profilePairs, int numProfiles);
Node **countPairs(Node **pairs, int numProfiles, ContigDescr *contigDescr, FILE *fp, int dist);

Node **getProfilePairs(int numProfiles, ContigDescr *contigDescr, FILE *fp, Args *args, int d){
  static Node **profilePairs = NULL;
  int i;

  numPos = 0;
  profilePairs = resetProfilePairs(profilePairs,numProfiles);
  if(args->L){
    for(i=0;i<args->S;i++)
      profilePairs = countPairs(profilePairs, numProfiles, contigDescr, fp, d+i);
  }else
    profilePairs = countPairs(profilePairs, numProfiles, contigDescr, fp, d);

  return profilePairs;
}

Node **resetProfilePairs(Node **profilePairs, int numProfiles){
  int i;

  if(profilePairs == NULL){
    profilePairs = (Node **)emalloc(numProfiles * sizeof(Node *));
    for(i=0;i<numProfiles;i++)
      profilePairs[i] = NULL;
  }else{
    for(i=0;i<numProfiles;i++){
      freeTree(profilePairs[i]);
      profilePairs[i] = NULL;
    }
  }
  return profilePairs;
}

Node **countPairs(Node **pairs, int numProfiles, ContigDescr *contigDescr, FILE *fp, int dist){
  int a, b, i, l, r, tmp, numRead;
  Position *pb;
  
  pb = contigDescr->posBuf;
  fseek(fp,3,SEEK_SET);
  for(i=0;i<contigDescr->n;i++){
    numRead = fread(pb,sizeof(Position),contigDescr->len[i],fp);
    assert(numRead == contigDescr->len[i]);
    l = 0;
    r = 0;
    while(r<contigDescr->len[i]){
      while(pb[r].pos - pb[l].pos < dist && r<contigDescr->len[i])
	r++;
      while(pb[r].pos - pb[l].pos > dist && l<=r)
  	l++;
      if(pb[r].pos-pb[l].pos == dist){
	a = pb[l].pro;
	b = pb[r].pro;
	if(a > b){ /* sort indexes to halve the number of index pairs */
	  tmp = b;
	  b = a;
	  a = tmp;
	}
	pairs[b] = addTree(pairs[b],a);
	numPos++;
      }
      r++;
      l++;
    }
  }
  return pairs;
}



/* addTree: add key to tree */
Node *addTree(Node *node, int key){

  if(node == NULL)  /* new key has arrived */
    node = newNode(key);
  else if(node->key == key)
    node->n++;      /* repeated key */
  else if(node->key < key) /* descend into left subtree */
    node->left = addTree(node->left, key);
  else              /* descend into right subtree */
    node->right = addTree(node->right, key);

  return node;
}

/*newNode: generate and initialize new node */
Node *newNode(int key){
  Node *node;

  node = (Node *)emalloc(sizeof(Node));
  node->key = key;
  node->n = 1;
  node->left = NULL;
  node->right = NULL;

  return node;
}

/* freeTree: traverse tree and free the memory for each node */
void freeTree(Node *n){
  if(n != NULL){
    freeTree(n->left);
    freeTree(n->right);
    free(n);
  }
}

double getNumPos(){
  return numPos;
}

