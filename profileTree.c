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

Node *getPairwProfiles(Node *root, Args *args, int d);
Node *getIndivProfiles(Node *root, int fd, Args *args);
Node *getSummarizedProfiles(Node *root, int fd, Args *args);
void setContigDescr(ContigDescr *contigDescr);
void testBinFile();

char testMode = 0;
double numPos;
ContigDescr *globalContigDescr;
int foundNodeIndex;
int numNode = 0;
Node **nodeArray = NULL;

Node *getProfileTree(int fd, Args *args, int d){
  Node *root;
  int i;

  root = NULL;
  numPos = 0;
  if(args->u){
    root = getSummarizedProfiles(root, fd, args);
    return root;
  }
  if(d){
    if(args->L){
      for(i=0;i<args->S;i++)
	root = getPairwProfiles(root, args, d+i);
    }else
      root = getPairwProfiles(root, args, d);
  }else
    root = getIndivProfiles(root, fd, args);

  return root;
}

Node *getSummarizedProfiles(Node *root, int fd, Args *args){
  FILE *fp;
  int count, i, status;
  char **strArr;
  Node *node;

  node = (Node *)emalloc(sizeof(Node));
  node->key = (char *)emalloc(256*sizeof(char));
  strArr = (char **)emalloc(4*sizeof(char *));
  for(i=0;i<4;i++)
    strArr[i] = (char *)emalloc(10*sizeof(char));
  fp = fdopen(fd,"r");
  if(fp == NULL){
    printf("ERROR: could not open file.\n");
    exit(-1);
  }
  while((status = fscanf(fp,"%d\t%s\t%s\t%s\t%s\n",&count,strArr[0],strArr[1],strArr[2],strArr[3])) != EOF){
    node->key[0] = '\0';
    node->c1 = 0;
    for(i=0;i<3;i++){
      strcat(node->key,strArr[i]);
      node->profile1[i] = atoi(strArr[i]);
      node->c1 += node->profile1[i];
    }
    strcat(node->key,strArr[3]);
    node->profile1[3] = atoi(strArr[3]);
    node->c1 += node->profile1[3];
    if(node->c1 >= args->c){
      numPos += count;
      root = addTree(root,node->key,count,node,NULL);
    }
  }
  for(i=0;i<4;i++)
    free(strArr[i]);
  free(strArr);
  fclose(fp);
  free(node->key);
  free(node);
  return root;
}

/* getProfileTree: get tree of profiles for specific profile distance, d, from file */
Node *getPairwProfiles(Node *root, Args *args, int d){
  int i, l, r, numRead;
  FILE *tmpF;
  char key[256];
  ContigDescr *contigDescr;
  Node *nodeL, *nodeR;
  Profile *pb;               /* profile buffer */
  

  tmpF = fopen(TMPFILE,"rb");
  contigDescr = getContigDescr();
  pb = contigDescr->profileBuf;
  for(i=0;i<contigDescr->n;i++){
    numRead = fread(pb,sizeof(Profile),contigDescr->len[i],tmpF);
    if(numRead != contigDescr->len[i]){
      printf("ERROR: getPariwProfiles\n");
      exit(-1);
    }
    l = 0;
    r = 0;
    while(r<contigDescr->len[i]){
      while(pb[r].pos - pb[l].pos < d && r<contigDescr->len[i])
	r++;
      while(pb[r].pos - pb[l].pos > d && l<=r)
  	l++;
      if(pb[r].pos-pb[l].pos == d){
	nodeL = nodeArray[pb[l].nodeIndex];
	nodeR = nodeArray[pb[r].nodeIndex];
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
  fclose(tmpF);
  if(args->r){
    printTree(stdout,root);
    exit (0);
  }
  return root;
}

Node *getIndivProfiles(Node *root, int fd, Args *args){
  char *buf, *line;
  QueueItem *qi;
  int n, i, j, l, max, maxNumProf, numProf;
  short headerOpen;
  FILE *tmpF;
  ContigDescr *contigDescr;
  Node *node;
  Profile profile;
  Profile *profileArray;

  tmpF = fopen(TMPFILE,"wb");
  contigDescr = (ContigDescr *)emalloc(sizeof(ContigDescr));
  node = (Node *)emalloc(sizeof(Node));
  node->key = (char *)emalloc(256*sizeof(char));
  node->left = NULL;
  node->right = NULL;
  contigDescr->n = 0;
  contigDescr->len = NULL;
  root = NULL;
  maxNumProf = 10000;
  profileArray = (Profile *)emalloc(maxNumProf*sizeof(Profile));
  buf = (char *)emalloc(args->b*sizeof(char));
  line = (char *)emalloc(100*sizeof(char));
  qi = newQueueItem();
  headerOpen = 0;
  l = 0;
  numProf = 0;
  while((n = read(fd, buf, args->b)) > 0){
    for(i=0; i<n; i++){
      if(buf[i] == '>'){
	headerOpen = 1;
	contigDescr->n++;
	contigDescr->len = (int *)erealloc(contigDescr->len,contigDescr->n*sizeof(int));
	contigDescr->len[contigDescr->n-1] = 0;
      }else if(headerOpen && buf[i] == '\n'){  /* reach end of header */
	headerOpen = 0;
      }else if(!headerOpen){
	if(buf[i] != '\n')
	  line[l++] = buf[i];
	else{                                  /* line filled */
	  line[l] = '\0';
	  l = 0;
	  fillQueueItem(line,qi);
	  node->key[0] = '\0';
	  node->c1 = 0;
	  for(j=0;j<3;j++){
	    node->profile1[j] = atoi(qi->profile[j]);
	    strcat(node->key,qi->profile[j]);
	    strcat(node->key," ");
	    node->c1 += node->profile1[j];
	  }
	  node->profile1[j] = atoi(qi->profile[j]);
	  strcat(node->key,qi->profile[j]);
	  node->c1 += node->profile1[j];
	  profile.pos = qi->pos;
	  if(node->c1 >= args->c){
	    numPos++;
	    contigDescr->len[contigDescr->n-1]++;
	    root = addTree(root,node->key,1,node,NULL);
	    profile.nodeIndex = foundNodeIndex;
	    profileArray[numProf++] = profile;
	    if(numProf == maxNumProf){
	      fwrite(profileArray,sizeof(Profile),numProf,tmpF);
	      numProf = 0;
	    }
	  }
	}
      }
    }
  }
  fwrite(profileArray,sizeof(Profile),numProf,tmpF);

  max = 0;
  for(i=0;i<contigDescr->n;i++)
    if(max<contigDescr->len[i])
      max = contigDescr->len[i];
  contigDescr->profileBuf = (Profile *)emalloc(max*sizeof(Profile));

  setContigDescr(contigDescr);
  freeQueueItem(qi);
  free(buf);
  free(line);
  fclose(tmpF);
  free(node->key);
  free(node);
  free(profileArray);
  return root;
}

void testBinFile(){
  int i, j, numRead;
  FILE *tmpF;
  ContigDescr *contigDescr;
  Profile *pb;               /* profile buffer */
  Node *np;

  tmpF = fopen(TMPFILE,"rb");
  contigDescr = getContigDescr();
  pb = contigDescr->profileBuf;
  numRead = 0;
  for(i=0;i<contigDescr->n;i++){
    numRead = fread(pb,sizeof(Profile),contigDescr->len[i],tmpF);
    printf("numRead: %d\n",numRead);
    for(j=0;j<contigDescr->len[i];j++){
      np = nodeArray[pb[j].nodeIndex];
      printf("%d\t%d\t%d\t%d\t%d\n",pb[j].pos,np->profile1[0],np->profile1[1],np->profile1[2],np->profile1[3]);
    }
  }
  fclose(tmpF);
  exit(0);
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
  numNode++;
  nodeArray = (Node **)erealloc(nodeArray,numNode*sizeof(Node *));
  nodeArray[numNode-1] = node;
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
  free(globalContigDescr->len);
  free(globalContigDescr->profileBuf);
  free(globalContigDescr);
  free(nodeArray);
}
