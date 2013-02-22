BEGIN{print ">Contig";
}{
  if($2 !~ /[^acgtACGT]/){
    s["A"]=0;
    s["C"]=0;
    s["G"]=0;
    s["T"]=0;
    s["a"]=0;
    s["c"]=0;
    s["g"]=0;
    s["t"]=0;
    for(i=1;i<=length($2);i++)
      s[substr($2,i,1)]++;
    print $1 "\t" s["A"]+s["a"] "\t" s["C"]+s["c"] "\t" s["G"]+s["g"] "\t" s["T"]+s["t"];
  }
}
