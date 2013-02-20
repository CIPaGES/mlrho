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
#include "queue.h"

Node *getPairwProfiles(FILE *fp, Node *root, Args *args, int d);
void setContigDescr(ContigDescr *contigDescr);

char testMode = 0;
double numPos;
ContigDescr *globalContigDescr = NULL;
int foundNodeIndex;
Node **nodeArray = NULL;

FILE *iniLinkAna(Args *args){
  char *fileName;
  char *tag;
  FILE *fp;
  ContigDescr *cp;
  int i, max, n, numRead;

  /* get contig lengths from file */
  tag = (char *)emalloc(4*sizeof(char));
  fileName = (char *)emalloc(256*sizeof(char));
  fileName = strcpy(fileName,args->n);
  fileName = strcat(fileName,".con");
  fp = efopen(fileName,"rb");
  numRead = fread(tag,sizeof(char),3,fp);
  assert(numRead == 3);
  tag[3] = '\0';
  if(strcmp(tag,"con") != 0)
    assert(0);
  cp = (ContigDescr *)emalloc(sizeof(ContigDescr));
  numRead = fread(&n,sizeof(int),1,fp);
  assert(numRead == 1);
  cp->n = n;
  cp->len = (int *)emalloc(n*sizeof(int));
  numRead = fread(cp->len,sizeof(int),n,fp);
  assert(numRead == n);
  fclose(fp);
  max = 0;
  for(i=0;i<n;i++){
    if(cp->len[i] > max)
      max = cp->len[i];
  }
  cp->posBuf = (Position *)emalloc(max*sizeof(Position));
  setContigDescr(cp);

  /* get file pointer for position file */
  fileName = strcpy(fileName,args->n);
  fileName = strcat(fileName,".pos");
  fp = efopen(fileName,"rb");
  numRead = fread(tag,sizeof(char),3,fp);
  if(strcmp(tag,"pos") != 0)
    assert(0);
  free(fileName);
  free(tag);
  return fp;
}

Node *getProfileTree(FILE *fp, Args *args, int d){
  Node *root;
  int i;

  root = NULL;
  numPos = 0;
  if(args->L){
    for(i=0;i<args->S;i++)
      root = getPairwProfiles(fp, root, args, d+i);
  }else
    root = getPairwProfiles(fp, root, args, d);

  return root;
}

/* getProfileTree: get tree of profiles for specific profile distance, d, from file */
Node *getPairwProfiles(FILE *fp, Node *root, Args *args, int d){
  int i, l, r, numRead;
  char key[256], tag[4];
  ContigDescr *contigDescr;
  Node *nodeL, *nodeR;
  Position *pb;               /* profile buffer */

  /* make sure you are at the start of a position file */
  rewind(fp);
  numRead = fread(&tag,sizeof(char),3,fp);
  assert(numRead == 3);
  tag[3] = '\0';
  if(strcmp(tag,"pos") != 0)
    assert(0);
  
  contigDescr = getContigDescr();
  pb = contigDescr->posBuf;

  for(i=0;i<contigDescr->n;i++){
    numRead = fread(pb,sizeof(Position),contigDescr->len[i],fp);
    assert(numRead == contigDescr->len[i]);
    l = 0;
    r = 0;
    while(r<contigDescr->len[i]){
      while(pb[r].pos - pb[l].pos < d && r<contigDescr->len[i])
	r++;
      while(pb[r].pos - pb[l].pos > d && l<=r)
  	l++;
      if(pb[r].pos-pb[l].pos == d){
	nodeL = nodeArray[pb[l].pro];
	nodeR = nodeArray[pb[r].pro];
  	if(nodeL->c1 >= args->c && nodeR->c1 >= args->c){
  	  strcpy(key,nodeL->key);
	  strcat(key," ");
	  strcat(key,nodeR->key);
  	  if(args->T){
  	    if(pb[l].pos % 2 != 0){
	      numPos++;
  	      root = addTree(root,key,1,nodeL,nodeR);
	    }
  	  }else{
	    numPos++;
  	    root = addTree(root,key,1,nodeL,nodeR);
  	  }
  	}
	r++;
	l++;
      }
    }
  }
  if(args->r){
    printTree(stdout,root);
    exit (0);
  }
  return root;
}

Node *getSummarizedProfiles(Args *args){
  int i, j, numRead, n;
  char *fileName, *tag, *strBuf;
  Node *node;
  Profile *profileArray;
  FILE *fp;

  tag = (char *)emalloc(4*sizeof(char));
  strBuf = (char *)emalloc(256*sizeof(char));
  tag[0] = '\0';
  fileName = (char *)emalloc(256*sizeof(char));
  fileName = strcpy(fileName,args->n);
  fileName = strcat(fileName,".sum");
  fp = efopen(fileName,"rb");
  
  numRead = fread(tag,sizeof(char),3,fp);
  assert(numRead == 3);
  tag[3] = '\0';
  if(strcmp(tag,"sum") != 0)
    assert(0);
  numRead = fread(&n,sizeof(int),1,fp);
  profileArray = (Profile *)emalloc(n*sizeof(Profile));
  numRead = fread(profileArray,sizeof(Profile),n,fp);
  assert(numRead == n);
  nodeArray = (Node **)emalloc(n*sizeof(Node *));
  for(i=0;i<n;i++){
    node = (Node *)emalloc(sizeof(Node));
    node->c1 = 0;
    for(j=0;j<4;j++){
      node->profile1[j] = profileArray[i].profile[j];
      node->c1 += node->profile1[j];
    }
    node->n = profileArray[i].n;
    numPos += node->n;
    node->key = (char *)emalloc(24*sizeof(char));
    itoa(node->profile1[0],strBuf);
    node->key = strcpy(node->key,strBuf);
    for(j=1;j<4;j++){
      node->key = strcat(node->key," ");
      itoa(node->profile1[j],strBuf);
      node->key = strcat(node->key,strBuf);
    }
    if(i)
      nodeArray[i-1]->left = node;
    node->left = NULL;
    node->right = NULL;
    node->c2 = 0;
    node->id = i;
    nodeArray[i] = node;
  }
  fclose(fp);
  free(fileName);
  free(profileArray);
  free(strBuf);

  return nodeArray[0];
}

void setContigDescr(ContigDescr *contigDescr){
  globalContigDescr = contigDescr;
}

ContigDescr *getContigDescr(){
  return globalContigDescr;
}

/* addTree: add key to tree */
Node *addTree(Node *node, char *key, int count, Node *n1, Node *n2){
  int cond;

  if(node == NULL){  /* new key has arrived */
    node = newNode(key,count,n1,n2);
    foundNodeIndex = node->id;
  }else if((cond = strcmp(key, node->key)) == 0){
    node->n++;      /* repeated key */
    foundNodeIndex = node->id;
  }else if(cond < 0) /* descend into left subtree */
    node->left = addTree(node->left, key, count, n1, n2);
  else              /* descend into right subtree */
    node->right = addTree(node->right, key, count, n1, n2);

  return node;
}

/*newNode: generate and initialize new node */
Node *newNode(char *key, int count, Node *n1, Node *n2){
  Node *node;
  int i, l;
  static int id = 0;

  node = (Node *)emalloc(sizeof(Node));
  node->id = id++;
  l = strlen(key);
  node->key = emalloc(sizeof(char)*(l+1));
  strcpy(node->key,key);
  node->n = count;
  node->left = NULL;
  node->right = NULL;
  for(i=0;i<4;i++)
    node->profile1[i] = n1->profile1[i];
  node->c1 = n1->c1;
  if(n2){
    for(i=0;i<4;i++)
      node->profile2[i] = n2->profile1[i];
    node->c2 = n2->c1;
  }else{
    node->c2 = 0;
    for(i=0;i<4;i++)
      node->profile2[i] = 0;
  }
  return node;
}

/* printTree: traverse tree and print profile & count for every node */
void printTree(FILE *fp, Node *node){

  if(node != NULL){
    printTree(fp, node->left);

    fprintf(fp,"%d|%s|%d %d %d %d",node->n,node->key,node->profile1[0],node->profile1[1],node->profile1[2],node->profile1[3]);
    if(node->profile2)
      fprintf(fp," %d %d %d %d",node->profile2[0],node->profile2[1],node->profile2[2],node->profile2[3]);
    fprintf(fp,"\n");

    printTree(fp, node->right);
  }
}

/* freeTree: traverse tree and free the memory for each node */
void freeTree(Node *n){
  if(n != NULL){
    freeTree(n->left);
    freeTree(n->right);
    free(n->key);
    free(n);
  }
}

void setTestMode(){
  testMode  = 1;
}

double getNumPos(){
  return numPos;
}

void freeProfileTree(){
  if(globalContigDescr){
    free(globalContigDescr->len);
    free(globalContigDescr->posBuf);
    free(globalContigDescr);
  }
  free(nodeArray);
}
