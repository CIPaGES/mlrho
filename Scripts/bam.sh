samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 1:1-20000000 | 
samtools mpileup - |
cut -f 2,5 |
awk -f ./Scripts/bam2pro.awk |
formatPro
mlRho -M 0 -I
mlRho -m 1000 -M 1005

