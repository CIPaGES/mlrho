#!/usr/bin/awk -f
# pro2sum.awk
# Purpose: Summarize profiles of nucleotide counts
# Usage: pro2sum.awk test.pro
# Author: Bernhard Haubold, haubold@evolbio.mpg.de
# Date: 12 September 2012
BEGIN{
  version = "0.1";
  if(ARGV[1] == "-h"){
    print "Usage: pro2sum.awk [input.pro]";
    print "Summarize profiles of nucleotide counts; distributed as part of mlRho.";
    print "Example: pro2sum.awk test.pro";
    exit;
  }
}
{
  if(!/>/){
    key = $2;
    for(i=3;i<=NF;i++)
      key = key "\t" $i;
    count[key]++;
  }
}END{
  for(a in count)
    print count[a] "\t" a;
 }