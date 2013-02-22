/***** profile.h **********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Oct 24 12:30:23 2012
 **************************************************/
#ifndef PROFILE
#define PROFILE

typedef struct profile{  /* profile written to disk: */
  int profile[4];        /* profile */
  int n;                 /* number of occurrences */
}Profile;

void readProfiles(char *baseName);
Profile *getProfiles();
int getNumProfiles();

#endif
