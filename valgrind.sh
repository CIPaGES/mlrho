#valgrind --leak-check=full --show-reachable=yes --log-file="mlRho.val" --leak-check=full ./mlRho -g 17:7512445-7513455 -B -M 3 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00154/alignment/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam 
#valgrind --leak-check=full --show-reachable=yes --log-file="mlRho.val" --leak-check=full ./mlRho -u -M 0 test.sum
valgrind --log-file="mlRho.val" ./mlRho test.pro -m 1000 -M 1005 -l 

