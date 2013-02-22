for i in $(seq 1000)
do
    ms 2 1 -t 1000 -r 500 100000    |
    ms2dna                          | 
    sequencer -c 4 -p               
done                                |
formatPro
mlRho -M 0 -I
mlRho -m 1000 -M 1005
