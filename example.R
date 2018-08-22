options(warn = -1)
.Machine$double.eps = .Machine$double.xmin

# leadSNP functions 
source("leadFunctions.R")

### INPUT ###

# 1. clumped gwas output
trans.in = "../data/META.TRANS.clumped"
eur.in   = "../data/META.EUR.clumped"
afr.in   = "../data/META.AFR.clumped"
sas.in   = "../data/META.SAS.clumped"

# 2. established snps (format chr1:54872)
ref = read.table("../ref/DIAMANTE_established_loci.txt", F, stringsAsFactors = F)

# 3. LD of all SNPs together (experimental plus established loci)
ld = read.table("META.REF.snps.ld", T, stringsAsFactors = F)

### PROCESS ###

# 1. import file
trans = clump.import(trans.in, "TRANS")
eur   = clump.import(eur.in, "EUR")
afr   = clump.import(afr.in, "AFR")
sas   = clump.import(sas.in, "SAS")

# 2. get coordinates
ref = snp.coordinates(ref)

# 3. register where all the SNPs are from (populations, established loci)
ld = ld.annotate(ld, df.ref = ref, df.trans = trans, df.eur = eur, df.afr = afr, df.sas = sas)

### ANALYZE ###

# 1. get block (region) of SNPs based on physical distance and p-value threshold
trans.block = snp.block(trans, half.window = 500000, p.threshold = 5e-8)
eur.block   = snp.block(eur,   half.window = 500000, p.threshold = 5e-8)
afr.block   = snp.block(afr,   half.window = 500000, p.threshold = 5e-8)
sas.block   = snp.block(sas,   half.window = 500000, p.threshold = 5e-8)

# 2. are SNPs novel as compared to established set of SNPs based on phsyical distance and LD
trans.novel = snp.novel(trans.block, df.ref = ref, df.ld = ld, half.window = 500000)
eur.novel   = snp.novel(eur.block,   df.ref = ref, df.ld = ld, half.window = 500000)
sas.novel   = snp.novel(sas.block,   df.ref = ref, df.ld = ld, half.window = 500000)
afr.novel   = snp.novel(afr.block,   df.ref = ref, df.ld = ld, half.window = 500000)

# 3. get GENE name(s) for the SNPs
trans.anno  = snp.annotate(trans.novel, range.1 = 100000, range.2 = 200000)
eur.anno    = snp.annotate(eur.novel,   range.1 = 100000, range.2 = 200000)
afr.anno    = snp.annotate(afr.novel,   range.1 = 100000, range.2 = 200000)
sas.anno    = snp.annotate(sas.novel,   range.1 = 100000, range.2 = 200000)

### EXPORT ###

write.table(trans.anno, "t2d.trans.novel.anno.txt", row.names = F, col.names = T, quote = F, sep = " ")
write.table(eur.anno,   "t2d.eur.novel.anno.txt",   row.names = F, col.names = T, quote = F, sep = " ")
write.table(afr.anno,   "t2d.afr.novel.anno.txt",   row.names = F, col.names = T, quote = F, sep = " ")
write.table(sas.anno,   "t2d.sas.novel.anno.txt",   row.names = F, col.names = T, quote = F, sep = " ")
