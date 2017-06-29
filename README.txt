1) Untar HOGFasta.tar.gz
2) HOGFasta, orthoxml file, *.txt and pipeline(_py2).sh have to be in the same folder
3) To run (example):
./pipeline.sh Traes_1DS 50 0.1 5 100 65 0.01
Traes_1DS = this is a unique substring with which starts every wheat gene on chr 1DS //inferring split genes on wheat chr 1DS
50 = min length of candidate genes (fragments)
0.1 = max overlap in the alignment
5 = min size of the gene family (2 candidates + 3 references)
100 = number of bootstrap samples
65 = threshold for collapsing (0.65)
0.01 = significance of the likelihood ratio test
