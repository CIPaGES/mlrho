/***** profileTreeBam.c ***************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sat Sep  8 04:50:05 2012
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
#include "sam.h"

double numPos;
int minCov;
Node *root;

void pileup(const char *inFile, Args *args, int d);

double getNumPosBam(){
  return numPos;
}

Node *getProfileTreeBam(char *inFile, Args *args, int d){
  /* int i; */

  root = NULL;
  numPos = 0;
  if(d){
    printf("ERROR: can only deal with single site data in this version\n");
    exit(0);
    /* if(args->L){ */
    /*   for(i=0;i<args->S;i++) */
    /* 	root = getPairwProfiles(root, fd, args, d+i); */
    /* }else */
    /*   root = getPairwProfiles(root, fd, args, d); */
  }else
    pileup(inFile, args, d);

  return root;
}

static int fetchFunc(const bam1_t *b, void *data){
  bam_plbuf_t *buf = (bam_plbuf_t *)data;
  bam_plbuf_push(b, buf);
  return 0;
}

char *profileToKey(int *profile, char *key){
  int i;
  char buf[20];

  key[0] = '\0';
  for(i=0;i<4;i++){
    itoa(profile[i],buf);
    strcat(key,buf);
  }
  return key;
}

static int pileupFuncIndiv(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data){
  int i, c, cov, profile[4], count[256], pair;
  tmpstruct_t *tmp;
  char *dic = "ACGT";
  char *key = "                                                                                    ";

  pair = 0;
  tmp = (tmpstruct_t *) data;
  if((int)pos >= tmp->beg && (int)pos < tmp->end && !pl->is_del){
    for(i=0;i<4;i++)
      count[(int)dic[i]] = 0;
    for(i=0;i<n;i++){
      const bam_pileup1_t *p = pl + i;
      c = bam1_seqi(bam1_seq(p->b), p->qpos);
      count[(int)bam_nt16_rev_table[c]]++;
    }
    cov = 0;
    for(i=0;i<4;i++){
      profile[i] = count[(int)dic[i]];
      cov += profile[i];
    }
    if(cov >= minCov){
      key = profileToKey(profile, key);
      numPos++;
      root = addTree(root,key,1,pair);
    }
  }
  return 0;
}


static int pileupFuncPair(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data){
  return 0;
}

void pileup(const char *inFile, Args *args, int d){
  tmpstruct_t tmp;
  int ref;
  bam_index_t *idx;
  bam_plbuf_t *buf;

  /* tmp->sf = samopen(inFile, "rb", 0); */
  tmp.in = NULL;
  tmp.beg = 0;
  tmp.end = 0x7fffffff;
  printf("inFile: %s\n",inFile);
  tmp.in = samopen(inFile, "rb", NULL);
  if(tmp.in == NULL){
    fprintf(stderr, "ERROR: failed to open BAM file %s.\n",inFile);
    exit(-1);
  }
  if(args->g == NULL){
    if(d)
      sampileup(tmp.in, -1, pileupFuncPair, &tmp);
    else
      sampileup(tmp.in, -1, pileupFuncIndiv, &tmp);
  }else{
    idx = bam_index_load(inFile);
    if(idx == NULL){
      fprintf(stderr, "ERROR: failed to locd BAM index.\n");
      exit(-1);
    }
    bam_parse_region(tmp.in->header, args->g, &ref, &(tmp.beg), &(tmp.end));
    if(ref < 0){
      fprintf(stderr, "ERROR: invalid region %s.\n", args->g);
      exit(-1);
    }
    if(d)
      buf = bam_plbuf_init(pileupFuncPair, &tmp);
    else
      buf = bam_plbuf_init(pileupFuncIndiv, &tmp);
    bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetchFunc);
    bam_plbuf_push(0,buf);
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);
  }
  samclose(tmp.in);
}
