/***** profileTree.c ******************************
 * Description: Tree for counting allele profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Feb 18 17:06:35 2009.
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "eprintf.h"
#include "tab.h"
#include "stringUtil.h"
#include "interface.h"
#include "profileTree.h"

Node *newNode(char *key, int *profile1, int *profile2, int d);
Node *addTree(Node *node, char *key, int *profile1, int *profile2, int d);
void setCoverage(Node *p, int d);
void adjustContig(ProfileCollection *pc, int ii, int jj);
void iniContig(ProfileCollection *pc, int ii);
ProfileCollection *newProfileCollection();

char testMode = 0;

ProfileCollection *readProfiles(FILE *fp){
  ProfileCollection *pc;
  char *key, *line, *field;
  int i, j, l, kl, ii, jj;
  int numProfiles;
  int prevPos, currPos, first;

  pc = newProfileCollection();
  key = (char *)emalloc(100*sizeof(char));
  pc->numSet = 0;
  pc->maxD = -1;
  ii = 0;
  numProfiles = 0;
  prevPos = 0;
  first = 1;
  while((line = tabGetLine(fp)) != NULL){
    if(line[0] == '>'){
      if(pc->maxD < numProfiles - 1)
	pc->maxD = numProfiles - 1;
      numProfiles = 0;
      ii = pc->numSet++;
      iniContig(pc,ii);
      first = 1;
    }else{
      if(tabNfield() != 5){
        printf("ERROR[readProfiles]: Line '%s' should consist of 5 fields but in fact contains %d fields - aborting.\n",line,tabNfield());
        exit(-1);
      }
      currPos = atoi(tabField(0));
      if(first){
	first = 0;
	prevPos = currPos;
      }
      /* padding positions */
      for(i=0;i<currPos-prevPos-1;i++){
	numProfiles++;
	jj = pc->setSize[ii]++;
	adjustContig(pc,ii,jj);
	for(j=0;j<4;j++)
	  pc->intDat[ii][jj][j] = 0;
	pc->strDat[ii][jj] = NULL;
	pc->cov[ii][jj] = 0;
	pc->pos[ii][jj] = -1;
      }
      prevPos = currPos;
      /* filled position */
      numProfiles++;
      jj = pc->setSize[ii]++;
      adjustContig(pc,ii,jj);
      kl = 0;
      pc->cov[ii][jj] = 0;
      pc->pos[ii][jj] = currPos;
      for(i=1;i<5;i++){ 
	field = tabField(i);
	pc->intDat[ii][jj][i-1] = atoi(field);
	pc->cov[ii][jj] += pc->intDat[ii][jj][i-1];
	l = strlen(field);
	for(j=0;j<l;j++)
	  key[kl++] = field[j];
	key[kl++] = '|';
      }
      key[kl] = '\0';
      pc->strDat[ii][jj] = strdup2(key);
    }
  }
  if(pc->maxD < numProfiles - 1)
    pc->maxD = numProfiles - 1;
  freeTab();
  free(key);
  return pc;
}

/* getProfileTree: get tree of profiles for specific profile distance, d */
Node *getProfileTree(Args *args, ProfileCollection *pc, int d){
  Node *root;
  char *key;
  int i, j;

  key = (char *)emalloc(100*sizeof(char));
  root = NULL;
  pc->np = 0;
  for(i=0;i<pc->numSet;i++){
    j = 0;
    while(j+d < pc->setSize[i]){
      if(pc->cov[i][j] >= args->c && pc->cov[i][j+d] >= args->c){
	pc->np++;
	key[0] = '\0';
	if(d){
	  key = strcat(key,pc->strDat[i][j]);
	  key = strcat(key,pc->strDat[i][j+d]);
	}else
	  key = strcat(key,pc->strDat[i][j]);
	root = addTree(root,key,pc->intDat[i][j],pc->intDat[i][j+d],d);
      }
      if(testMode)
	j += 2;
      else
	j++;
    }
  }
  setCoverage(root,d);
  free(key);
  return root;
}

/* setCoverage: Recursion across tree to set coverage per profile */
void setCoverage(Node *p, int d){
  int i;

  if(p != NULL){
    setCoverage(p->left, d);
    p->c1 = 0;
    p->c2 = 0;
    for(i=0;i<4;i++)
      p->c1 += p->profile[i];
    if(d)
      for(i=4;i<8;i++)
	p->c2 += p->profile[i];

    setCoverage(p->right, d);
  }
}

void freeProfiles(ProfileCollection *pc){
  int i, j;
  for(i=0;i<pc->numSet;i++){
    free(pc->cov[i]);
    free(pc->pos[i]);
    for(j=0;j<pc->setSize[i];j++){
      free(pc->intDat[i][j]);
      free(pc->strDat[i][j]);
    }
    free(pc->intDat[i]);
    free(pc->strDat[i]);
  }
  free(pc->cov);
  free(pc->pos);
  free(pc->intDat);
  free(pc->strDat);
  free(pc->setSize);
  free(pc);
}

/* addTree: add key to tree */
Node *addTree(Node *node, char *key, int *profile1, int *profile2, int d){
  int cond;

  if(node == NULL)  /* new key has arrived */
    node = newNode(key, profile1, profile2, d);
  else if((cond = strcmp(key, node->key)) == 0)
    node->n++;      /* repeated key */
  else if(cond < 0) /* descend into left subtree */
    node->left = addTree(node->left, key, profile1, profile2, d);
  else              /* descend into right subtree */
    node->right = addTree(node->right, key, profile1, profile2, d);
  return node;
}

/*newNode: generate and initialize new node */
Node *newNode(char *key, int *profile1, int *profile2, int d){
  Node *node;
  int i, j;

  node = (Node *)emalloc(sizeof(Node));
  node->key = strdup2(key);
  for(i=0;i<4;i++)
    node->profile[i] = profile1[i];
  if(d){
    j = 0;
    for(i=4;i<8;i++)
      node->profile[i] = profile2[j++];
  }else
    for(i=4;i<8;i++)
      node->profile[i] = 0;
  node->n = 1;
  node->left = NULL;
  node->right = NULL;

  return node;
}

/* printTree: traverse tree and print profile and count for every node */
void printTree(FILE *fp, Node *node){
  int i;

  if(node != NULL){
    printTree(fp, node->left);

    fprintf(fp,"%d",node->profile[0]);
    for(i=1;i<8;i++)
      fprintf(fp,"\t%d",node->profile[i]);
    fprintf(fp,"\t%d\n",node->n);

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

void printProfiles(ProfileCollection *pc){
  int i, j, k;

  for(i=0;i<pc->numSet;i++){
    printf(">Contig_%d\n",i+1);
    for(j=0;j<pc->setSize[i];j++){
      if(pc->cov[i][j]){
	printf("%d",pc->pos[i][j]);
	for(k=0;k<4;k++)
	  printf("\t%d",pc->intDat[i][j][k]);
	printf("\t%s\n",pc->strDat[i][j]);
      }
    }
  }
}

void setTestMode(){
  testMode  = 1;
}

void adjustContig(ProfileCollection *pc, int ii, int jj){
  pc->intDat[ii] = (int **)erealloc(pc->intDat[ii],pc->setSize[ii]*sizeof(int *));
  pc->strDat[ii] = (char **)erealloc(pc->strDat[ii],pc->setSize[ii]*sizeof(char *));
  pc->cov[ii] = (int *)erealloc(pc->cov[ii],pc->setSize[ii]*sizeof(int));
  pc->pos[ii] = (int *)erealloc(pc->pos[ii],pc->setSize[ii]*sizeof(int));
  pc->intDat[ii][jj] = (int *)emalloc(5*sizeof(int));
}

void iniContig(ProfileCollection *pc, int ii){
  if(pc->setSize == NULL)
    pc->setSize = (int *)emalloc(pc->numSet*sizeof(int));
  else
    pc->setSize = (int *)erealloc(pc->setSize,pc->numSet*sizeof(int));
  if(pc->intDat == NULL)
    pc->intDat = (int ***)emalloc(pc->numSet*sizeof(int **));
  else
    pc->intDat = (int ***)erealloc(pc->intDat,pc->numSet*sizeof(int **));
  if(pc->strDat == NULL)
    pc->strDat = (char ***)emalloc(pc->numSet*sizeof(char **));
  else
    pc->strDat = (char ***)erealloc(pc->strDat,pc->numSet*sizeof(char **));
  if(pc->cov == NULL)
    pc->cov = (int **)emalloc(pc->numSet*sizeof(int *));
  else
    pc->cov = (int **)erealloc(pc->cov,pc->numSet*sizeof(int *));
  if(pc->pos == NULL)
    pc->pos = (int **)emalloc(pc->numSet*sizeof(int *));
  else
    pc->pos = (int **)erealloc(pc->pos,pc->numSet*sizeof(int *));

  pc->setSize[ii] = 0;
  pc->intDat[ii] = (int **)emalloc(pc->setSize[ii]*sizeof(int *));
  pc->strDat[ii] = (char **)emalloc(pc->setSize[ii]*sizeof(char *));
  pc->pos[ii] = (int *)emalloc(pc->setSize[ii]*sizeof(int));
  pc->cov[ii] = (int *)emalloc(pc->setSize[ii]*sizeof(int));
}

ProfileCollection *newProfileCollection(){
  ProfileCollection *pc;

  pc = (ProfileCollection *)emalloc(sizeof(ProfileCollection));
  pc->setSize = NULL;
  pc->intDat = NULL;
  pc->strDat = NULL;
  pc->cov = NULL;
  pc->pos = NULL;
  
  return pc;
}

/* pairFreq: compute the frequency of homo- and
 * heterozygous pairs of profiles with distance d
 */
void pairFreq(ProfileCollection *pc, int d){
  int i, j, k;
  int  c1, c2, numPairs;

  pc->homPairFreq = 0;
  pc->hetPairFreq = 0;
  numPairs = 0;
  for(i=0;i<pc->numSet;i++){
    for(j=0;j<pc->setSize[i]-d;j++){
      numPairs++;
      c1 = 0;
      c2 = 0;
      for(k=0;k<4;k++){
	if(pc->intDat[i][j][k])
	  c1++;
	if(pc->intDat[i][j+d][k])
	  c2++;
      }
      if(c1>1 && c2>1)
	pc->hetPairFreq++;
      if(c1==1 && c2==1)
	pc->homPairFreq++;
    }
  }
  pc->hetPairFreq /= (double)numPairs;
  pc->homPairFreq /= (double)numPairs;
}
