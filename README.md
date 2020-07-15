---
title: "pleiotest package: Sequential test for detecting genetic pleiotropy"
author: "Fernando Aguate"
date: "7/15/2020"
---

pleiotest is a package that provides tools for detecting genetic associations with multiple traits (i.e. pleiotropy).

It performs a multi-trait genome-wide association analysis with seemingly unrelated regressions based on generalized least squares (GLS). Results from this model are then used to run a sequential Wald test with some considerations to formally test for pleiotropic effects. 
The package offers some computational advantages that allows it to handle large and unbalanced data sets. In addition, it has functions and arguments to include covariates, subset the data, save the results or plot them.


```R
library(pleiotest)
```

This code creates fake balanced data to use as an example

```R
# Set seed to obtain consistent results
set.seed(2020)
# Number of traits, individuals and SNPs
n_traits <- 3
n_individuals <- 10000
n_snps <- 5000
# Create pheno numeric matrix with traits in columns
pheno <- matrix(rnorm(n_traits * n_individuals), ncol = n_traits)
# Set fake names for the fake traits
colnames(pheno) <- c('tic', 'tac', 'toe')
# Creeate geno numeric matrix with SNPs in columns
geno <- sapply(1:n_snps, function(i) rbinom(n_individuals, 2, runif(1, 0.01, .49)))
# the row names of geno and pheno must be matching IDs
rownames(geno) <- rownames(pheno) <- paste0('ind', 1:n_individuals)
# Column names of geno can be rs SNP IDs
colnames(geno) <- paste0('rsid', 1:ncol(geno))
# The pleioR function takes a "melted" version of pheno
library(reshape2)
pheno_m <- melt(pheno)
```

Run the analysis with pleioR

```R
pleio_object <- pleioR(pheno = pheno_m, geno = geno)
# To obtain by trait estimates of the multi-trait GWAS use mt_gwas
mt_result <- mt_gwas(pleio_object)
# head of the table with results for the 2nd trait with fake name "tac"
head(mt_result$tac)
```

Sequential test for pleiotropy

```R
pleio_results <- pleio_test(pleio_object)
# pleio_results is a list of p values, indices indicating traits in association, and labels for each index
head(pleio_results$pValues)
```
```R
# If p1 in pValues is significant, ind1 indicates which is the associated trait
# If p2 in pValues is significant, ind2 indicates which are the two associated traits
head(pleio_results$Index)
```
```R
# Labels for traits
pleio_results$traits
```

It's possible to indexing SNPs and individuals in geno (useful when running in parallel)

```R
# Here we index 2 SNPs and 50% of individuals in geno
pleio_object2 <- pleioR(pheno = pheno_m, geno = geno,
                           i = sample(x = seq_len(nrow(geno)), size = nrow(geno) * .5),
                           j = 1:2)
# Use save_at to save results in a folder. 
mt_result2 <- mt_gwas(pleio_object2, save_at = '~/Documents/')
# This creates a file named mt_gwas_result_x.rdata
file.exists('~/Documents/mt_gwas_result_1.rdata')
```

```R
pleio_result2 <- pleio_test(pleio_object2, save_at = '~/Documents/')
# results will be saved in pleio_test_result_x.rdata
file.exists('~/Documents/pleio_test_result_1.rdata')
# p values of the indexed SNPs
head(pleio_result2$pValues)
```

pleioR can also deal with unbalanced data

```R
# sample 10% of data to remove
missing_values <- sample(1:nrow(pheno_m), round(nrow(pheno_m) * .1))
pheno_m2 <- pheno_m[-missing_values,]
# The randomly generated unbalance created the following sub-sets:
identify_subsets(trait = pheno_m2$Var2, id = pheno_m2$Var1)[[1]]
```

Fit the model with unbalanced data

```R
# It's possible to drop sub-sets with less than x observations (drop_subsets = x)
pleio_object3 <- pleioR(pheno = pheno_m2, geno = geno, drop_subsets = 200)
# to save computation time, stop the sequence for p-values larger than loop_breaker
system.time({pleio_result3 <- pleio_test(pleio_object3, loop_breaker = .99)})
system.time({pleio_result3 <- pleio_test(pleio_object3, loop_breaker = .05)})
system.time({pleio_result3 <- pleio_test(pleio_object3, loop_breaker = .01)})
```

It's also possible to add covariates (optional)

```R
covariates <- matrix(rnorm(n_individuals * 3, 1, 2), ncol = 3)
rownames(covariates) <-  rownames(geno)
pleio_object4 <- pleioR(pheno = pheno_m2, geno = geno, covariates = covariates)
pleio_result4 <- pleio_test(pleio_object4)
head(pleio_result4$pValues)
```

Plotting the results of pleio_test using base pair positions

```R
# Fake base pair positions and centromeres positions
bp_positions <- data.frame('chr' = rep(1:22, length.out = nrow(pleio_result4[[1]])), 
                           'pos' = seq(1, 1e8, length.out = nrow(pleio_result4[[1]])), 
                           row.names = rownames(pleio_result4[[1]]))
centromeres = aggregate(pos ~ chr, data = bp_positions, FUN = function(x) mean(x) / 1e6)

# The function pleio_plot also returns a table with the significant SNPs
pleio_plot(pleio_res = pleio_result4, n_traits = 2, alpha = .05, bp_positions = bp_positions, chr_spacing = 1000)
```
```R
# The function pleio_ideogram also returns a table with the plotted regions
pleio_ideogram(pleio_res = pleio_result4, n_traits = 2, alpha = .05, bp_positions = bp_positions, window_size = 1e6, centromeres = centromeres, set_ylim_prop = 1.3)
```
