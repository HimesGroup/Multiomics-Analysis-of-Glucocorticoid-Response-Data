# Multiomics-Analysis-of-Glucocorticoid-Response-Data
Supporting codes and data of the manuscript "Multiomics Analysis Identifies BIRC3 as a Novel Glucocorticoid Response-Associated Gene"

## Description
In this paper, we ranked inhaled corticosteroid (ICS) response-associated variants by assigning them multicomis integrative scores derived from five omics layers (GR-binding, GRE, RNAP II-binding, differential expression and eQTL) and determined significance based on permutaitons. We identified four SNPs near the gene *BIRC3* having significant scores that may influence ICS response in people with asthma.

## Files

* GSK.score.txt: scores of 4,468 GSK variants based on 1000 Genome Project EUR reference panel in each omics layer

* GSK.score.perm.res: permutation results of 4,468 GSK variants based on 1000 Genome Project EUR reference panel

* GSK_cos.score.txt: scores of 2,822 GSK variants based on 1000 Genome Project cosmopolitan reference panel in each omics layer

* GSK_cos.score.perm.res: permutation results of 2,822 GSK variants based on 1000 Genome Project cosmopolitan reference panel

* GSK.cand.score.txt: scores of 4,468 GSK and 4,271 non-GSK variants in each omics layer

* GSK.cand.score.perm.res: permutation results of 4,468 GSK and 4,271 non-GSK variants

* multiomics.score.permutation.R: R codes to generate multiomics integrative scores and perform permutations

## Run codes
Rscript multiomics.score.permutation.R --score_fn [input score file] --out_fn [result file] --n_core 20 --nperm 1000000
