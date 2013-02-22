/***** pairs.h ************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Oct 26 15:46:50 2012
 **************************************************/
typedef pairs{
  int n;             /* number of profiles */
  int *singleCounts; /* singleCounts[i]: number of profiles pairs (i,j) */
  int **ids;         /* ids[i][j]: identity the j-th pair (i,.) */
  int **counts;      /* counts[i][j]: count of the j-th pair (i,.) */
  int *maxSizes;     /* maxSizes[i]: max. number of entries in pairIds[i] and pairCounts[i] */
}Pairs;
