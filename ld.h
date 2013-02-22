/***** ld.h ***************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Oct 25 08:54:30 2012
 **************************************************/
#ifndef LDHEADER
#define LDHEADER

typedef struct position{
  int pos;
  int pro;
}Position;

typedef struct contigDescr{
  int n;                 /* number of contigs */
  int *len;              /* contig lengths */
  Position *posBuf;      /* buffer of positions */
}ContigDescr;

ContigDescr *getContigDescr();
FILE *iniLdAna(Args *args);
ContigDescr *getContigDescr();
#endif
