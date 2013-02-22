/***** profile.c **********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Oct 24 12:24:34 2012
 **************************************************/
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "eprintf.h"
#include "profile.h"

void setNumProfiles(int numProfiles);
void setProfiles(Profile *profiles);

Profile *thisProfiles;
int thisNumProfiles;

void readProfiles(char *baseName){
  char *fileName, tag;
  int i, numProfiles, numRead, sop;
  Profile *profiles;
  FILE *fp;

  fileName = (char *)emalloc(256*sizeof(char));
  fileName = strcpy(fileName,baseName);
  fileName = strcat(fileName,".sum");

  fp = efopen(fileName,"rb");
  for(i=0;i<3;i++)
    numRead = fread(&tag,sizeof(char),1,fp);
  assert(numRead == 1);
  numRead = fread(&numProfiles,sizeof(int),1,fp);
  assert(numRead == 1);
  sop = sizeof(Profile);
  profiles = (Profile *)emalloc(numProfiles*sop);
  for(i=0;i<numProfiles;i++){
    numRead = fread(&profiles[i],sop,1,fp);
    assert(numRead == 1);
  }
  setProfiles(profiles);
  setNumProfiles(numProfiles);
  fclose(fp);
  free(fileName);
}

void setProfiles(Profile *profiles){
  thisProfiles = profiles;
}

Profile *getProfiles(){
  return thisProfiles;
}

void setNumProfiles(int numProfiles){
  thisNumProfiles = numProfiles;
}

int getNumProfiles(){
  return thisNumProfiles;
}
