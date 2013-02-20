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
#include "eprintf.h"
#include "tab.h"
#include "stringUtil.h"
#include "interface.h"
#include "profileTree.h"
#include "queue.h"

Node *getPairwProfiles(Node *root, int fd, Args *args, int d);
Node *getIndivProfiles(Node *root, int fd, Args *args);
Node *getSummarizedProfiles(Node *root, int fd, Args *args);

char testMode = 0;
double numPos;

Node *getProfileTree(int fd, Args *args, int d)
{
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
	root = getPairwProfiles(root, fd, args, d+i);
    }else
      root = getPairwProfiles(root, fd, args, d);
  }else
    root = getIndivProfiles(root, fd, args);

  return root;
}

Node *getSummarizedProfiles(Node *root, int fd, Args *args){
  FILE *fp;
  int count, i, status, pair;
  char **strArr, *key;

  pair = 0;
  strArr = (char **)emalloc(4*sizeof(char *));
  for(i=0;i<4;i++)
    strArr[i] = (char *)emalloc(10*sizeof(char));
  key = (char *)emalloc(100*sizeof(char));
  fp = fdopen(fd,"r");
  if(fp == NULL){
    printf("ERROR: could not open file.\n");
    exit(-1);
  }
  while((status = fscanf(fp,"%d\t%s\t%s\t%s\t%s\n",&count,strArr[0],strArr[1],strArr[2],strArr[3])) != EOF){
    key[0] = '\0';
    for(i=0;i<3;i++){
      strcat(key,strArr[i]);
      strcat(key, " ");
    }
    strcat(key,strArr[3]);
    if(coverage(key, 0) >= args->c){
      numPos += count;
      root = addTree(root,key,count,pair);
    }
  }
  for(i=0;i<4;i++)
    free(strArr[i]);
  free(strArr);
  free(key);
  fclose(fp);
  return root;
}

/* getProfileTree: get tree of profiles for specific profile distance, d, from file */
Node *getPairwProfiles(Node *root, int fd, Args *args, int d)
{
  char *buf, *key, *line;
  QueueItem *end, *start, **qiStore;
  int n, i, j, l, numStore, pair;
  short headerOpen;

  pair = 1;
  lseek(fd, 0L, 0);
  buf = (char *)emalloc(args->b*sizeof(char));
  line = (char *)emalloc(100*sizeof(char));
  key = (char *)emalloc(100*sizeof(char));
  headerOpen = 0;
  l = 0;
  qiStore = (QueueItem **)emalloc((d+1)*sizeof(QueueItem *));
  for(i=0;i<d+1;i++)
    qiStore[i] = newQueueItem();
  queueInit(d+1);
  numStore = d + 1; 
  while((n = read(fd, buf, args->b)) > 0){
    for(i=0; i<n; i++){
      if(buf[i] == '>'){
	headerOpen = 1;
	while(!queueEmpty())
	    qiStore[numStore++] = queueGet();
      }else if(headerOpen && buf[i] == '\n'){  /* reach end of header */
	headerOpen = 0;
      }else if(!headerOpen){
	if(buf[i] != '\n')
	  line[l++] = buf[i];
	else{                                  /* line filled */
	  line[l] = '\0';
	  l = 0;
	  if(numStore > 0)
	      end = qiStore[--numStore];
	  else{
	      printf("WARNING: profile positions not in ascending order.\n");
	      end = newQueueItem();
	  }
	  fillQueueItem(line,end);
	  queuePut(end);
	  start = queuePeek();
	  while((start->pos+d) <= end->pos){
	    start = queueGet();
	    qiStore[numStore++] = start;
	    if(end->pos == start->pos+d){
	      key[0] = '\0';
	      for(j=0;j<4;j++){
		strcat(key,start->profile[j]);
		strcat(key, " ");
	      }
	      for(j=0;j<3;j++){
		strcat(key,end->profile[j]);
		strcat(key, " ");
	      }
	      strcat(key,end->profile[j]);
	      if(coverage(key,d) >= args->c){
		if(args->T){
		  if(start->pos % 2 != 0){
		    numPos++;
		    root = addTree(root,key,1,pair);
		  }
		}else{
		  numPos++;
		  root = addTree(root,key,1,pair);
		}
	      }
	    }
	    start = queuePeek();
	  }
	}
      }
    }
  }
  freeQueue();
  while(numStore > 0)
      freeQueueItem(qiStore[--numStore]);
  free(qiStore);
  free(buf);
  free(line);
  free(key);

  return root;
}

Node *getIndivProfiles(Node *root, int fd, Args *args)
{
  char *buf, *key, *line;
  QueueItem *qi;
  int n, i, j, l, pair;
  short headerOpen;
  
  root = NULL;
  buf = (char *)emalloc(args->b*sizeof(char));
  line = (char *)emalloc(100*sizeof(char));
  key = (char *)emalloc(100*sizeof(char));
  qi = newQueueItem();
  headerOpen = 0;
  pair = 0;
  l = 0;
  while((n = read(fd, buf, args->b)) > 0){
    for(i=0; i<n; i++){
      if(buf[i] == '>'){
	headerOpen = 1;
      }else if(headerOpen && buf[i] == '\n'){  /* reach end of header */
	headerOpen = 0;
      }else if(!headerOpen){
	if(buf[i] != '\n')
	  line[l++] = buf[i];
	else{                                  /* line filled */
	  line[l] = '\0';
	  l = 0;
	  fillQueueItem(line,qi);
	  key[0] = '\0';
	  for(j=0;j<3;j++){
	    strcat(key,qi->profile[j]);
	    strcat(key, " ");
	  }
	  strcat(key,qi->profile[j]);
	  if(coverage(key, 0) >= args->c){
	    numPos++;
	    root = addTree(root,key,1,pair);
	  }
	}
      }
    }
  }
  freeQueueItem(qi);
  free(buf);
  free(line);
  free(key);

  return root;
}


int coverage(char *key, int d){
  int profile1[4], profile2[4];
  int c1, c2, i, minCov;

  c1 = 0;
  c2 = 0;
  if(d){
    sscanf(key,"%d %d %d %d %d %d %d %d",&profile1[0],&profile1[1],&profile1[2],&profile1[3], \
	   &profile2[0],&profile2[1],&profile2[2],&profile2[3]);
    for(i=0;i<4;i++){
      c1 += profile1[i];
      c2 += profile2[i];
    }
    minCov = (c1 < c2) ? c1 : c2;
  }else{
    sscanf(key,"%d %d %d %d",&profile1[0],&profile1[1],&profile1[2],&profile1[3]);
    for(i=0;i<4;i++)
      c1 += profile1[i];
    c2 = 0;
    minCov = c1;
  }

  return minCov;
}

/* addTree: add key to tree */
Node *addTree(Node *node, char *key, int count, int pair){
  int cond;

  if(node == NULL)  /* new key has arrived */
    node = newNode(key,count,pair);
  else if((cond = strcmp(key, node->key)) == 0)
    node->n++;      /* repeated key */
  else if(cond < 0) /* descend into left subtree */
    node->left = addTree(node->left, key, count, pair);
  else              /* descend into right subtree */
    node->right = addTree(node->right, key, count, pair);

  return node;
}

/*newNode: generate and initialize new node */
Node *newNode(char *key, int count, int pair){
  Node *node;
  int n, i;

  node = (Node *)emalloc(sizeof(Node));
  n = strlen(key) + 1;
  node->key = (char *)emalloc(n*sizeof(char));
  node->key = strncpy(node->key,key,n);
  node->n = count;
  node->left = NULL;
  node->right = NULL;
  node->profile1 = (int *)emalloc(sizeof(int)*4);
  if(pair){
    node->profile2 = (int *)emalloc(sizeof(int)*4);
    sscanf(node->key,"%d %d %d %d %d %d %d %d",&node->profile1[0],&node->profile1[1],&node->profile1[2],&node->profile1[3], \
	   &node->profile2[0],&node->profile2[1],&node->profile2[2],&node->profile2[3]);
  }else{
    node->profile1 = (int *)emalloc(sizeof(int)*4);
    sscanf(node->key,"%d %d %d %d",&node->profile1[0],&node->profile1[1],&node->profile1[2],&node->profile1[3]);
  }
  node->c1 = 0;
  node->c2 = 0;
  for(i=0;i<4;i++)
    node->c1 += node->profile1[i];
  if(pair){
    for(i=0;i<4;i++)
      node->c2 += node->profile2[i];
  }
  return node;
}

/* printTree: traverse tree and print profile & count for every node */
void printTree(FILE *fp, Node *node){

  if(node != NULL){
    printTree(fp, node->left);

    fprintf(fp,"%s %d\n",node->key,(int)strlen(node->key));

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
