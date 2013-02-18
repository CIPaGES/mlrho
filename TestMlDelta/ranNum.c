/***** ranNum.c **************************************************************
 * Description: Random number generation for mlDivEst
 * Author: Michael Lynch
 * Modified by: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Jan 23 16:08:53 2009.
 *****************************************************************************/
#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>


#define hetero   0.0001			/* average nucleotide heterozygosity */
#define phi      0.00046		/* read-error rate */
#define coverage 4			/* average sequencing coverage per site */
#define nsites	 100000			/* number of informative sites */
#define niters	 1000			/* number of individuals iterated */



extern long ignbin(long n,double pp);

/* Definitions for the binomial generator. */
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/* Definitions for the random number generator. */
#define TRUE -1
#define FALSE 0
#define boolean int

static double u[98], c, cd, cm;
static int i97, j97;
static boolean test = FALSE;

double ranmar(void)
/*
  C This is the random number generator proposed by George Marsaglia in 
  C Florida State University Report: FSU-SCRI-87-50
  C It was slightly modified by F. James to produce an array of pseudorandom
  C numbers.
*/
{
  double uni;
        
  if (test==FALSE) {
    puts("Call the init routine rmarin() before calling ranmar().");
    exit(2);
  }
  uni = u[i97] - u[j97];
  if (uni < 0.0) uni += 1.0;
  u[i97] = uni;
  i97--;
  if (i97==0) i97 = 97;
  j97--;
  if (j97==0) j97 = 97;
  c -= cd;
  if (c<0.0) c += cm;
  uni -= c;
  if (uni<0.0) uni += 1.0;
  return uni;
}




/* Seed for the random number generator. */

void rmarin(int ij,int kl) 
{
  /*
    C This is the initialization routine for the random number generator RANMAR()
    C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
    C                                                      0 <= KL <= 30081
    C The random number sequences created by these two seeds are of sufficient 
    C length to complete an entire calculation with. For example, if sveral 
    C different groups are working on different parts of the same calculation,
    C each group could be assigned its own IJ seed. This would leave each group
    C with 30000 choices for the second seed. That is to say, this random 
    C number generator can create 900 million different subsequences -- with 
    C each subsequence having a length of approximately 10^30.
    C 
    C Use IJ = 1802 & KL = 9373 to test the random number generator. The
    C subroutine RANMAR should be used to generate 20000 random numbers.
    C Then display the next six random numbers generated multiplied by 4096*4096
    C If the random number generator is working properly, the random numbers
    C should be:
    C           6533892.0  14220222.0  7275067.0
    C           6172232.0  8354498.0   10633180.0
  */
  int i, j, k, l, ii, jj, m;
  double s, t;
        
  if (ij<0 || ij>31328 || kl<0 || kl>30081) {
    puts("The first random number seed must have a value between 0 and 31328.");
    puts("The second seed must have a value between 0 and 30081.");
    exit(1);
  }
        
  i = (ij/177)%177 + 2;
  j = ij%177 + 2;
  k = (kl/169)%178 + 1;
  l = kl%169;
        
  for (ii=1; ii<=97; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj=1; jj<=24; jj++) {
      m = (((i*j)%179)*k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l + 1) % 169;
      if ((l*m)%64 >= 32) s += t;
      t *= 0.5;
    }
    u[ii] = s;
  }
        
  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;
        
  i97 = 97;
  j97 = 33;
        
  test = TRUE;
}




/* ****************************************************************************************** */

/* binomial random number generator */

long ignbin(long n,double pp)
/*
**********************************************************************
long ignbin(long n,double pp)
GENerate BINomial random deviate
Function
Generates a single random deviate from a binomial
distribution whose number of trials is N and whose
probability of an event in each trial is P.
Arguments
n  --> The number of trials in the binomial distribution
from which a random deviate is to be generated.
p  --> The probability of an event in each trial of the
binomial distribution from which a random deviate
is to be generated.
ignbin <-- A random deviate yielding the number of events
from N independent trials, each of which has
a probability of event P.
Method
This is algorithm BTPE from:
Kachitvichyanukul, V. and Schmeiser, B. W.
Binomial Random Variate Generation.
Communications of the ACM, 31, 2
(February, 1988) 216.
**********************************************************************
SUBROUTINE BTPEC(N,PP,ISEED,JX)
BINOMIAL RANDOM VARIATE GENERATOR
MEAN .LT. 30 -- INVERSE CDF
MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
(SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
USABLE ALGORITHM.
REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
"BINOMIAL RANDOM VARIATE GENERATION,"
COMMUNICATIONS OF THE ACM, FORTHCOMING
WRITTEN:  SEPTEMBER 1980.
LAST REVISED:  MAY 1985, JULY 1987
REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
GENERATOR
ARGUMENTS
N : NUMBER OF BERNOULLI TRIALS            (INPUT)
PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
VARIABLES
PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
M:  INTEGER VALUE OF THE CURRENT MODE
FM:  FLOATING POINT VALUE OF THE CURRENT MODE
XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
P1:  AREA OF THE TRIANGLE
C:  HEIGHT OF THE PARALLELOGRAMS
XM:  CENTER OF THE TRIANGLE
XL:  LEFT END OF THE TRIANGLE
XR:  RIGHT END OF THE TRIANGLE
AL:  TEMPORARY VARIABLE
XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
P2:  AREA OF THE PARALLELOGRAMS
P3:  AREA OF THE LEFT EXPONENTIAL TAIL
P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
FROM THE REGION
V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
(REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
REJECT THE CANDIDATE VALUE
IX:  INTEGER CANDIDATE VALUE
X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
K:  ABSOLUTE VALUE OF (IX-M)
F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
ALSO USED IN THE INVERSE TRANSFORMATION
R: THE RATIO P/Q
G: CONSTANT USED IN CALCULATION OF PROBABILITY
MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
OF F WHEN IX IS GREATER THAN M
IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
CALCULATION OF F WHEN IX IS LESS THAN M
I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
YNORM: LOGARITHM OF NORMAL BOUND
ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
USED IN THE FINAL ACCEPT/REJECT TEST
QN: PROBABILITY OF NO SUCCESS IN N TRIALS
REMARK
IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
ARE NOT INVOLVED.
ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/

{
  static float psave = -1.0;
  static long nsave = -1;
  static long ignbin2,i,ix,ix1,k,m,mp,T1;
  static double al,alv,amaxp,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u2,v,w,w2,x,x1,
    x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

  if(pp != psave) goto S10;
  if(n != nsave) goto S20;
  if(xnp < 30.0) goto S150;
  goto S30;
 S10:
  /*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
*/
  psave = pp;
  p = min(psave,1.0-psave);
  q = 1.0-p;
 S20:
  xnp = n*p;
  nsave = n;
  if(xnp < 30.0) goto S140;
  ffm = xnp+p;
  m = ffm;
  fm = m;
  xnpq = xnp*q;
  p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
  xm = fm+0.5;
  xl = xm-p1;
  xr = xm+p1;
  c = 0.134+20.5/(15.3+fm);
  al = (ffm-xl)/(ffm-xl*p);
  xll = al*(1.0+0.5*al);
  al = (xr-ffm)/(xr*q);
  xlr = al*(1.0+0.5*al);
  p2 = p1*(1.0+c+c);
  p3 = p2+c/xll;
  p4 = p3+c/xlr;
 S30:
  /*
*****GENERATE VARIATE
*/
  u2 = ranmar()*p4;
  v = ranmar();
  /*
    TRIANGULAR REGION
  */
  if(u2 > p1) goto S40;
  ix = xm-p1*v+u2;
  goto S170;
 S40:
  /*
    PARALLELOGRAM REGION
  */
  if(u2 > p2) goto S50;
  x = xl+(u2-p1)/c;
  v = v*c+1.0-ABS(xm-x)/p1;
  if(v > 1.0 || v <= 0.0) goto S30;
  ix = x;
  goto S70;
 S50:
  /*
    LEFT TAIL
  */
  if(u2 > p3) goto S60;
  ix = xl+log(v)/xll;
  if(ix < 0) goto S30;
  v *= ((u2-p2)*xll);
  goto S70;
 S60:
  /*
    RIGHT TAIL
  */
  ix = xr-log(v)/xlr;
  if(ix > n) goto S30;
  v *= ((u2-p3)*xlr);
 S70:
  /*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
  k = ABS(ix-m);
  if(k > 20 && k < xnpq/2-1) goto S130;
  /*
    EXPLICIT EVALUATION
  */
  f = 1.0;
  r = p/q;
  g = (n+1)*r;
  T1 = m-ix;
  if(T1 < 0) goto S80;
  else if(T1 == 0) goto S120;
  else  goto S100;
 S80:
  mp = m+1;
  for(i=mp; i<=ix; i++) f *= (g/i-r);
  goto S120;
 S100:
  ix1 = ix+1;
  for(i=ix1; i<=m; i++) f /= (g/i-r);
 S120:
  if(v <= f) goto S170;
  goto S30;
 S130:
  /*
    SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
  */
  amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
  ynorm = -(k*k/(2.0*xnpq));
  alv = log(v);
  if(alv < ynorm-amaxp) goto S170;
  if(alv > ynorm+amaxp) goto S30;
  /*
    STIRLING'S FORMULA TO MACHINE ACCURACY FOR
    THE FINAL ACCEPTANCE/REJECTION TEST
  */
  x1 = ix+1.0;
  f1 = fm+1.0;
  z = n+1.0-fm;
  w = n-ix+1.0;
  z2 = z*z;
  x2 = x1*x1;
  f2 = f1*f1;
  w2 = w*w;
  if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
								     (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
																     (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
																								    (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
																																    -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
  goto S30;
 S140:
  /*
    INVERSE CDF LOGIC FOR MEAN LESS THAN 30
  */
  qn = pow(q,(double)n);
  r = p/q;
  g = r*(n+1);
 S150:
  ix = 0;
  f = qn;
  u2 = ranmar();
 S160:
  if(u2 < f) goto S170;
  if(ix > 110) goto S150;
  u2 -= f;
  ix += 1;
  f *= (g/ix-r);
  goto S160;
 S170:
  if(psave > 0.5) ix = n-ix;
  ignbin2 = ix;
  return ignbin2;
}

