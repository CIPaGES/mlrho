/***** interface.h **********************************************************
 * Description: Header file for user interface.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun 20 13:07:28 2004.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/ 

#define DEFAULT_S 100000
#define DEFAULT_P 0.00046
#define DEFAULT_H 0.0001
#define DEFAULT_C 4
#define DEFAULT_I 1000


/* define argument container */
typedef struct args{
  char *d; /* data file */
  float P; /* error rate, phi */
  float H; /* heterozygosity */
  int c;   /* coverage */
  int s;   /* number of sites */
  int i;   /* number of iterations */
  char p;  /* print program information */
  char h;  /* help message? */
  char e;  /* error message? */
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);