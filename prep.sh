# Take the output from the METAL file (or MR-MEGA)

# extract chr:bp and p-value, and add header
echo "SNP P" > example.toClump
awk '{if(NR>1) print $1" "$6}' mydata.1.TBL >> example.toClump

# clump SNPs from experimental GWAS output
/software/PLINK/PLINK_1.9/plink \
        --bfile ../ldref/EUR.ref \
        --clump example.toClump \
        --clump-kb 500 \
        --clump-p1 0.00000005 \
        --clump-p2 0.00001 \
        --clump-r2 0.05 \
        --clump-field P \
        --out example.distinct.locus

# clean up the clumped file
awk -F' ' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' example.distinct.locus.clumped > example.distinct.locus.txt

# combine established and novel SNPs into one file for LD calculation
awk 'FNR>1 {print $1}' example.establishedSNPs.txt > example.establishedSNPs.and.experimentalSNPs.txt
awk 'FNR>1 {print $3}' example.distinct.locus.txt >> example.establishedSNPs.and.experimentalSNPs.txt

# get LD of all established and potentially novel SNPs
/software/PLINK/PLINK_1.9/plink \
       --bfile ../ldref/EUR.ref \
       --r2 dprime \
       --ld-window-r2 0.05 \
       --extract example.establishedSNPs.and.experimentalSNPs.txt \
       --out example.established.and.experimental

