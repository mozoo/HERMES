# HERMES
HERMES is a straightforward index which tries to summarize the mitochondrial evolution pace using a single number.
Several mitogenomic features are evaluated in a factor analysis framework; namely, in the current version, they
are %URs, AMIGA, SU skew, root-to-tip distance and pairwise ML distance.

The HERMES index is a method of quantifying molecular evolution of mitochondrial genomes
in different species and clusters; this method was originally proposed in Plazzi et al.
(2016). The index relies on maximum likelihood factor analysis to summarize different measures
that are typically found to be linked with evolutionary rates; it is intended to be computed
a posteriori, i.e. after a phylogenetic and genomic analysis. As different empirical measures
are merged together in a single score, it is a “hyper-empirical” index; moreover, it is a
relative measure, because all species are compared with an outgroup: therefore, it was called
Hyper-Empirical Relative Mitochondrial Evolutionary Speed (HERMES) index.

The present Python script performs data collection and then calls a dedicated R script to
complete the factor analysis. The mitogenomic features that are currently implemented in
HERMES.py are:

the percentage of Unassigned Regions (URs);

the absolute value of the Strand Usage (SU) skew;

the Amount of Mitochondrial Identical Gene Arrangements (AMIGA) index;

the root-to-tip distance computed over a given phylogenetic tree;

the Maximum Likelihood (ML) pairwise distance from a given outgroup.

Normalization and varimax rotation are used, factor scores are found using correlation
preserving, and correlations are found using the Pearson method; given the possible presence
of missing values, missing data are set to be imputed using the median.

All the variables are pooled together for each species into the value of a single loading: we
define this score as the HERMES score of a given species.

The best-performing variable set and the goodness-of-fit of the analysis is assessed following
the recommendations of Hu and Bentler (1999): Tucker-Lewis Index (TLI) greater than 0.95; root
mean square of the residuals (SRMR) smaller than 0.08; root mean squared error of approximation
(RMSEA) less than 0.06; moreover, the  Kaiser-Meyer-Olkin index (KMO) was taken into account on
this regard.
