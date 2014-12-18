###  Description
This is part of the supplementary material of the paper: 

*Prokaryotic ancestry of a dual localized peroxiredoxin 5 in malaria parasites.*
 Carine F. Djuika, Jaime Huerta-Cepas, Jude M. Przyborski, Sophia Deil, Cecilia
 P. Sanchez, Tobias Doerks, Peer Bork, Michael Lanzer, Marcel
 Deponte. **Microbial Cell (2014). In press***

The repository contains all the data and scripts used for phylogenetic analyses.

```
hyp_testing                         -> Analysis of 96 constraint topologies using CONSEL
mrbayes_raw                         -> Source and configuration files used in MrBayes analysis
n-terminus.alg.fa                   -> Multipl sequence alignment of 3 prx5 sequences, PfAOP including the N-terminus region
prx5.alg.fa                         -> Multiple sequence alignment of the 84 prx5 sequences
prx5.py                             -> Python script for tree analyses
raxml_100_seeds                     -> Raxml test for suboptimal topologies
raxml_raw                           -> Raxml full scan test, with 100 bootstraps (raw files)
raxml_tree.nw                       -> tree from Raxml full scan 
mrbayes_tree.nw                     -> tree obtained from Bayesian full scan 
raxml_trimmed_alg                   -> Raxml scan using a trimmed version of the alignment
suboptimal_topologies.summary.txt   -> Condensed summary of all suboptimal topologies found 
bootstrapped_topologies.summary.txt -> Condensed summary of all bootstrapped topologies found
```

### Notes:

- All text files use unix new line separator
- Large intermediate files obtained in some of the analyses are not included 
