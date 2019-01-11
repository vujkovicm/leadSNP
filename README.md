# leadSNP

Identify list of potential novel SNPs associated with disease. Requires 1) clumped GWAS output, 2) list of established loci (from literature), 3) plink LD out for all clumped SNPs (experimental + established)

## Example on how to preformat the input files from using prep.sh

# R functions 
clump.import, snp.block, snp.novel, ld.annotate, snp.annotate, snp.coordinates
