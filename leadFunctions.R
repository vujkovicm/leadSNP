#=====================================================================================
#
#           Filename :  leadFunctions.R
#
#       Description  :  identify list of potential novel SNPs associated with disease
#                       requires * clumped GWAS output (from meta-analysis)
#                                * list of established loci (from literature)
#                                * plink LD output for all clumped SNPs (from analysis + established)
#
#        R functions :  clump.import, snp.block, snp.novel, ld.annotate, snp.annotate, snp.coordinates
#
#            Version :  1.0
#            Created :  21-Aug-2018
#           Revision :  none
#  Last modification :  21-Aug-2018
#
#             Author :  Marijana Vujkovic
#        Modified by :
#        Institution :  University of Pennsylvania, Perelman School of Medicine
#              Email :  vujkovic@pennmedicine.upenn.edu
#            License :  GPL (>=3)
#
#=====================================================================================

# import clumped file
clump.import = function(filename = "", pop) {
   if(file.exists(filename)) {
      library("data.table")
      df <- fread(filename, select = c("CHR", "SNP", "BP", "P", "S001", "S0001"))
      df[which(df$P == 0), "P"] = 1e-315
      df$POP = pop
      return(as.data.frame(df, stringsAsFactors = F))
   }
}


# get coordinates based on format chr1:58984
snp.coordinates = function(df){
        names(df) = "CHRCBP"
        tmp = as.data.frame(do.call(rbind, strsplit(df$CHRCBP, '\\:')), stringsAsFactors = F)
        colnames(tmp) = c("CHR", "BP")
        df = cbind(df, tmp)
        df$CHR = as.numeric(substring(df$CHR, 4))
        df$BP  = as.numeric(df$BP)
        return(df)
}

snp.block = function(df, half.window = 500000, p.threshold = 5e-8) {
  df = df[order(df$CHR, df$BP), ]
  out = NULL
  for (iCHR in 1:22)  {
     df.chr = df[df[, "CHR"] == iCHR & df[, "P"] < p.threshold, c("SNP", "BP", "P", "S001", "S0001", "POP")]
     if (nrow(df.chr) > 0) {
        for (i in 1:nrow(df.chr)) {
           current.bp  = df.chr[i, "BP"]
           current.df  = df.chr[df.chr[, "BP"] >= current.bp - half.window & df.chr[, "BP"] <= current.bp + half.window, ]
           block.start = min(current.df[, "BP"])
           block.end   = max(current.df[, "BP"])
           block.p     = min(current.df[, "P"])
           block.bp    = current.df[current.df[, "P"] == block.p, "BP"]
           block.S001  = sum(current.df[, "S001"])
           block.S0001 = sum(current.df[, "S0001"])
           block.snp   = paste0("chr", iCHR, ":", block.bp)
           tmp = data.frame(CHR = iCHR, BP = current.bp, SNP =  df.chr[i, "SNP"], P =  df.chr[i, "P"], S001 = df.chr[i, "S001"], S0001 = df.chr[i, "S0001"],
                POP = df.chr[i, "POP"], BLOCK.SNP = block.snp, BLOCK.P = block.p, BLOCK.S001 = block.S001, BLOCK.S0001 = block.S0001,
                BLOCK.START = block.start, BLOCK.END = block.end, stringsAsFactors = F)
           out = rbind(out, tmp)
         }
      }
   }
   out$BLOCK.PRIMARY   = ifelse(out$SNP == out$BLOCK.SNP, 1, 0)
   out$BLOCK.SECONDARY = ifelse(out$SNP != out$BLOCK.SNP, 1, 0)
   return(out)
}


# is a SNP novel as compared to a list of established SNPs?
snp.novel = function(df, chr = "CHR", bp = "BP", snp = "SNP", df.ref, df.ld, half.window = 500000) {
   out = NULL
   for (iCHR in 1:22) {
      df.chr   = df[df[, chr] == iCHR, ]
      ref.chr  = df.ref[df.ref[, "CHR"] == iCHR, ]
      if (nrow(df.chr) > 0) {
        df.chr$SENTINEL = NA
        df.chr$ISOLATED = 0
        for (i in 1:nrow(df.chr)) {
           if(nrow(ref.chr) == 0) {
              df.chr$SENTINEL[i] = NA
              df.chr$ISOLATED[i] = 0
           }
           else {
              for(j in 1:nrow(ref.chr)) {
                 if((df.chr[i, bp] >= ref.chr[j, bp] - half.window) & (df.chr[i, bp] <= ref.chr[j, bp] + half.window)) {
                    df.chr$SENTINEL[i] = ref.chr[j, "CHRCBP"]
                 }
                 else if (df.chr[i, snp] %in% df.ld$SNP_A) {
                    if(ref.chr[j, "CHRCBP"] %in% df.ld[which(df.ld$SNP_A == df.chr[i, snp]), "SNP_B"]) {
                       df.chr$SENTINEL[i] = ref.chr[j, "CHRCBP"]
                    }
                 }
                 else if (df.chr[i, snp] %in% df.ld$SNP_B) {
                    if(ref.chr[j, "CHRCBP"] %in% df.ld[which(df.ld$SNP_B == df.chr[i, snp]), "SNP_A"]) {
                       df.chr$SENTINEL[i] = ref.chr[j, "CHRCBP"]
                    }
                 }
              }
           }
           if (is.na(df.chr$SENTINEL[i]) == T) {
              if (df.chr$S0001[i] <= df.chr$S001[i]) {
                 df.chr$ISOLATED[i] = 1
                 df.chr$SENTINEL[i] = "."
              }
           }
        }
        out = rbind(out, df.chr)
      }
   }
   out$SENTINEL = ifelse(is.na(out$SENTINEL), "-", out$SENTINEL)
   out$NOVEL    = ifelse(out$SENTINEL == "-", 1, 0)
   return(out)
}


## LD: which SNPs are stored where
ld.annotate = function(df.ld, df.ref, df.trans = F, df.eur = F, df.sas = F, df.afr = F, df.amr = F, df.asn = F){
   df.ld[, "REF_A"] = ifelse(df.ld$SNP_A %in% df.ref[, "CHRCBP"], 1, 0)
   df.ld[, "REF_A"] = ifelse(df.ld$SNP_A %in% df.ref[, "CHRCBP"], 1, 0)
   if(is.data.frame(df.trans)){
      df.ld[, "TRANS_A"] = ifelse(df.ld$SNP_A %in% df.trans[, "SNP"], 1, 0)
      df.ld[, "TRANS_B"] = ifelse(df.ld$SNP_B %in% df.trans[, "SNP"], 1, 0)
   }
   if(is.data.frame(df.eur)){
      df.ld[, "EUR_A"] = ifelse(df.ld$SNP_A %in% df.eur[, "SNP"], 1, 0)
      df.ld[, "EUR_B"] = ifelse(df.ld$SNP_B %in% df.eur[, "SNP"], 1, 0)
   }
   if(is.data.frame(df.afr)){
      df.ld[, "AFR_A"] = ifelse(df.ld$SNP_A %in% df.afr[, "SNP"], 1, 0)
      df.ld[, "AFR_B"] = ifelse(df.ld$SNP_B %in% df.afr[, "SNP"], 1, 0)
   }
   if(is.data.frame(df.sas)){
      df.ld[, "SAS_A"] = ifelse(df.ld$SNP_A %in% df.sas[, "SNP"], 1, 0)
      df.ld[, "SAS_B"] = ifelse(df.ld$SNP_B %in% df.sas[, "SNP"], 1, 0)
   }
   if(is.data.frame(df.amr)){
      df.ld[, "AMR_A"] = ifelse(df.ld$SNP_A %in% df.amr[, "SNP"], 1, 0)
      df.ld[, "AMR_B"] = ifelse(df.ld$SNP_B %in% df.amr[, "SNP"], 1, 0)
   }
   if(is.data.frame(df.asn)){
      df.ld[, "ASN_A"] = ifelse(df.ld$SNP_A %in% df.asn[, "SNP"], 1, 0)
      df.ld[, "ASN_B"] = ifelse(df.ld$SNP_B %in% df.asn[, "SNP"], 1, 0)
   }
   return(df.ld)
}


# check if the variants from race-specific meta-analysis are in the overall-transethnic meta-analyisis
# if not, they are distinct variants
snp.distinct = function(df, df.ref = trans, snp.name = "SNP", pop = "EUR", ref = "TRANS") {
   pop.A = paste0(pop, "_A")
   pop.B = paste0(pop, "_B")
   flag.distinct = paste0(pop, ".distinct")
   for(i in 1:nrow(df)) {
      if (df[i, snp.name] %in% df.ref[, snp.name]){
         df$TRANS.SENTINEL[i] = df[i, snp.name]
      }
      else {
         if ((df[i, snp.name] %in% ld$SNP_A) && (ld[which(ld$SNP_A == df[i, snp.name]), pop.A] == 1) && (ld[which(ld$SNP_A == df[i, snp.name]), "TRANS_B"] == 1)) {
            df$TRANS.SENTINEL[i] = ld[which(ld$SNP_A == df[i, snp.name]), "SNP_B"]
         }
         else {
            if ((df[i, snp.name] %in% ld$SNP_B) && (ld[which(ld$SNP_B == df[i, snp.name]), pop.B] == 1) && (ld[which(ld$SNP_B == df[i, snp.name]), "TRANS_A"] == 1)) {
               df$TRANS.SENTINEL[i] = ld[which(ld$SNP_B == df$CHRCBP[i]), "SNP_A"]
            }
            else {
               df$TRANS.SENTINEL[i] = NA
            }
         }
       }
    }
    df[, flag.distinct] = ifelse(is.na(df$TRANS.SENTINEL) == T, 1, 0)
    return(df)
}


# get nearest gene for each SNP
# requires dataframe that contains at least chromosome (1) and basepair location (57865)
snp.annotate = function(df, chr = "CHR", bp = "BP", range.1 = 100000, range.2 = 200000) {
   .packages = c("biomaRt")
   .inst <- .packages %in% installed.packages()
   if(length(.packages[!.inst]) > 0) {
      install.packages(.packages[!.inst])
   }
   lapply(.packages, require, character.only = TRUE)
   # load the respository
   mart.hs  = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
   # initiate search, is SNP in gene?
   df$query = paste(gsub("chr", '', df[, chr]), df[, bp], df[, bp], sep = ":")
   df$ENSID = NA
   df$GENE = ""
   for(i in 1:nrow(df)) {
      out = getBM(
      attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'chromosomal_region', values = df$query[i], mart = mart.hs)
      out = subset(out, !grepl("\\.", out$external_gene_name))
      out = subset(out, !grepl("LINC", out$external_gene_name))
      if(nrow(out) > 0) {
         df$ENSID[i] = out$ensembl_gene_id
         df$GENE[i] = out$external_gene_name
      }
      else { # if SNP not in gene, extend range
        df$query[i] = paste(gsub("chr", '', df[i, chr]), df[i, bp] - range.1, df[i, bp] + range.1, sep = ":")
        out = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'chromosomal_region', values = df$query[i], mart = mart.hs )
        out = subset(out, !grepl("\\.", out$external_gene_name))
        out = subset(out, !grepl("LINC", out$external_gene_name))
        if(nrow(out) > 0) {
           df$ENSID[i] = NA
           df$GENE[i]  =  paste(unique(unlist(out$external_gene_name)), collapse = ';')
        }
        else { # if not in first range, then try second range
           df$query[i] = paste(gsub("chr", '', df[i, chr]), df[i, bp] - range.2, df[i, bp] + range.2, sep = ":")
           out = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'chromosomal_region', values = df$query[i], mart = mart.hs)
           out = subset(out, !grepl("\\.", out$external_gene_name))
           out = subset(out, !grepl("LINC", out$external_gene_name))
           df$ENSID[i] = NA
           df$GENE[i]  =  paste(unique(unlist(out$external_gene_name)), collapse = ';')
        }
      }
   }
   return(df)
}
