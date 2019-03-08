#=====================================================================================
#
#           Filename :  leadFunctions.R
#
#       Description  :  identify list of potential novel SNPs associated with disease
#                       requires * clumped GWAS output (from meta-analysis)
#                                * list of established loci (from literature)
#                                * plink LD output for all clumped SNPs (from analysis + established)
#
#        R functions :  clump.import, snp.block, snp.novel, ld.annotate, snp.annotate, snp.coordinates, ancestry.tag
#                          overlapping.cojo.snps, updated.primary.lead, blocks.revisited, clean.cojo, clean.vep
#
#            Version :  1.1
#            Created :  21-Aug-2018
#           Revision :  new function ancestry.restricted, fixed snp.novel
#  Last modification :  08-Mar-2019
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
      df <- fread(filename, select = c("CHR", "SNP", "BP", "P"))
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
        tmp$CHR = ifelse(tmp$CHR == "chrX", "chr24", tmp$CHR)
        df = cbind(df, tmp)
        df$CHR = as.numeric(substring(df$CHR, 4))
        df$BP  = as.numeric(df$BP)
        return(df)
}

# define independend regions based on LD and physical distance
snp.block = function(df, df.ld, half.window = 500000, p.threshold = 5e-8) {
  df = df[order(df$CHR, df$BP), ]
  out = NULL
  for (iCHR in 1:22)  {
     #print(paste("chr", iCHR, sep = ":"))
     df.chr = df[df[, "CHR"] == iCHR & df[, "P"] < p.threshold, c("SNP", "BP", "P", "POP")]
     if (nrow(df.chr) > 0) {
        for (i in 1:nrow(df.chr)) {
           #print(i)
           current.bp  = df.chr[i, "BP"]
           current.snp = df.chr[i, "SNP"]
           current.df  = df.chr[df.chr[, "BP"] >= current.bp - half.window & df.chr[, "BP"] <= current.bp + half.window, ]
           # add SNPs that are in R2 with current.bp
           r2.A        = df.ld[which(df.ld$SNP_A == current.snp), "SNP_B"]
           r2.B        = df.ld[which(df.ld$SNP_B == current.snp), "SNP_A"]
           r2          = unique(c(r2.A, r2.B))
           if(length(r2) > 0) {
                r2.df  = df.chr[which(df.chr$SNP %in% r2),]
                # merge region based on physical distance and r2
                current.df  = unique(rbind(current.df, r2.df))
           }
           block.start = min(current.df[, "BP"])
           block.end   = max(current.df[, "BP"])
           block.p     = min(current.df[, "P"])
           block.bp    = current.df[current.df[, "P"] == block.p, "BP"]
           block.snp   = paste0("chr", iCHR, ":", block.bp)
           tmp = data.frame(CHR = iCHR, BP = current.bp, SNP =  df.chr[i, "SNP"], P =  df.chr[i, "P"], POP = df.chr[i, "POP"],
                BLOCK.SNP = block.snp, BLOCK.BP = block.bp, BLOCK.P = block.p, BLOCK.START = block.start, BLOCK.END = block.end, stringsAsFactors = F)
           out = rbind(out, tmp)
         }
      }
   }
   out$BLOCK.PRIMARY   = ifelse(out$SNP == out$BLOCK.SNP, 1, 0)
   out$BLOCK.SECONDARY = ifelse(out$SNP != out$BLOCK.SNP, 1, 0)
   # clean up the blocks
   count = 0
   for (i in 1:nrow(out)) {
        if (out$BLOCK.PRIMARY[i] == 1) {
                count = count + 1
                out$BLOCK.N[i] = count
        } else {
                out$BLOCK.N[i] = NA
        }
   }
   # now merge blocks that are knd of close and give them an id
   for (i in 2:nrow(out)) {
        if (is.na(out$BLOCK.N[i]) == T) {
                prv = out[which((out$BLOCK.PRIMARY == 1) & (out$BLOCK.N == max(out$BLOCK.N[1:i], na.rm = T))), "SNP"]
                if(i < nrow(out)) {
                        nxt = out[which((out$BLOCK.PRIMARY == 1) & (out$BLOCK.N == min(out$BLOCK.N[i:nrow(out)], na.rm = T))), "SNP"]
                }
                if((nxt == prv) & (i < nrow(out))) {
                        out$BLOCK.N[i] = out$BLOCK.N[i - 1]
                }
                else {
                if(out[i, "BLOCK.SNP"] == out[which(out$SNP == prv), "BLOCK.SNP"]) {
                        out$BLOCK.N[i] = out$BLOCK.N[i - 1]
                } else if(out[i, "BLOCK.SNP"] == out[which(out$SNP == nxt), "BLOCK.SNP"]) {
                        out$BLOCK.N[i] = out$BLOCK.N[i - 1] + 1
                } else  if (out[i, "CHR"] - out[i - 1, "CHR"] != 0) {
                        out$BLOCK.N[i] = out$BLOCK.N[i - 1] + 1
                } else if (i < nrow(out)) {
                        if (out[i, "CHR"] - out[i + 1, "CHR"] != 0) {
                                out$BLOCK.N[i] = out$BLOCK.N[i - 1]
                        } else {
                                prev.diff = abs(out[i, "BP"] - out[which(out$SNP == prv), "BP"])
                                next.diff = abs(out[i, "BP"] - out[which(out$SNP == nxt), "BP"])
                                if(prev.diff < next.diff) {
                                        out$BLOCK.N[i] = out$BLOCK.N[i - 1]
                                } else {
                                        out$BLOCK.N[i] = out$BLOCK.N[i - 1] + 1
                                }
                                }
                        }
                }
        }
   }
   # update according to the new block information the positions
   out.final = NULL
   for(i in 1:max(out$BLOCK.N)) {
        block             = out[which(out$BLOCK.N == i), ]
        block$BLOCK.SNP   = block[which(block$BLOCK.P == min(block$BLOCK.P)), "BLOCK.SNP"][[1]] # SNP with minimum P-value
        if(length(block[which(block$SNP == tail(names(sort(table(block$BLOCK.SNP))), 1)[[1]]), "BP"]) > 0) {
                block$BLOCK.BP    = block[which(block$BLOCK.PRIMARY == 1), "BP"]
                block$BLOCK.P     = block[which(block$BLOCK.PRIMARY == 1), "P"]
                block$BLOCK.START = min(block$BP)
                block$BLOCK.END   = max(block$BP)
                out.final         = rbind(out.final, block)
        } else if(i < out[which(out$SNP == block$BLOCK.SNP), "BLOCK.N"]){ # leadSNP is in next block
                out[which(out$BLOCK.N == i), "BLOCK.PRIMARY"]   = 0
                out[which(out$BLOCK.N == i), "BLOCK.SECONDARY"] = 1
                out[which(out$BLOCK.N == i), "BLOCK.P"]         = out[which(out$SNP == block$BLOCK.SNP), "P"]
                out[which(out$BLOCK.N == i), "BLOCK.BP"]        = out[which(out$SNP == block$BLOCK.SNP), "BP"]
                out[which(out$BLOCK.N == i), "BLOCK.N"]         = out[which(out$SNP == block$BLOCK.SNP), "BLOCK.N"]
                out[which(out$BLOCK.N == i), "BLOCK.SNP"]       = out[which(out$SNP == block$BLOCK.SNP), "BLOCK.SNP"]
        } else {
                out.final = subset(out.final, out.final$BLOCK.N != (i - 1)) # remove previous entry from out.final
                out[which(out$BLOCK.N == i), "BLOCK.PRIMARY"]   = 0
                out[which(out$BLOCK.N == i), "BLOCK.SECONDARY"] = 1
                out[which(out$BLOCK.N == i), "BLOCK.P"]         = out[which(out$SNP == block$BLOCK.SNP), "P"]
                out[which(out$BLOCK.N == i), "BLOCK.BP"]        = out[which(out$SNP == block$BLOCK.SNP), "BP"]
                out[which(out$BLOCK.N == i), "BLOCK.N"]         = out[which(out$SNP == block$BLOCK.SNP), "BLOCK.N"]
                out[which(out$BLOCK.N == i), "BLOCK.SNP"]       = out[which(out$SNP == block$BLOCK.SNP), "BLOCK.SNP"]
                i = i - 1
        }
    }
    # clean up the block ids
    count = 0
    for (i in 1:nrow(out.final)) {
        if (out.final$BLOCK.PRIMARY[i] == 1) {
                count = count + 1
                out.final$BLOCK.N[i] = count
        } else {
                out.final$BLOCK.N[i] = NA
        }
    }
    for (i in 1:nrow(out.final)){ # fill in the missings
        if(is.na(out.final$BLOCK.N[i]) == T) {
                out.final$BLOCK.N[i] = out.final[which(out.final$SNP == out.final$BLOCK.SNP[i]), "BLOCK.N"]
        }
    }
    return(out.final)
}

# need to put this as seperate, because need to merge with COJO
overlapping.cojo.snps = function(df) {
        tmp = unique(df[which(df$SecondarySNP %in% df[duplicated(df$SecondarySNP), "SecondarySNP"] & is.na(df$SecondarySNP) == F), "SecondarySNP"])
        if(length(tmp) == 0) {
                return(0)
        } else {
                return(tmp)
        }
}

# get the new primary
updated.primary.lead = function(df) {
        tmp = unique(df[which(df$SecondarySNP %in% df[duplicated(df$SecondarySNP), "SecondarySNP"] & is.na(df$SecondarySNP) == F), "SecondarySNP"])
        if(length(tmp) == 0) {
                return(0)
        } else {
                out = NULL
                for (i in 1:length(tmp)) {
                        df.tmp    = df[which(df$SecondarySNP == tmp[i]), c("BLOCK.SNP", "BLOCK.P")]
                        new.lead  = df.tmp[which(df.tmp$BLOCK.P == min(df.tmp$BLOCK.P)), "BLOCK.SNP"]
                        out       = c(out, new.lead)
                }
                return(unique(out))
        }
}

# return the loci that COJO should be rerun on
updated.primary.lead = function(df) {
        tmp = unique(df[which(df$SecondarySNP %in% df[duplicated(df$SecondarySNP), "SecondarySNP"] & is.na(df$SecondarySNP) == F), "SecondarySNP"])
        if(length(tmp) == 0) {
                return(0)
        } else {
                out = NULL
                for (i in 1:length(tmp)) {
                        df.tmp    = df[which(df$SecondarySNP == tmp[i]), c("BLOCK.SNP", "BLOCK.P")]
                        new.lead  = df.tmp[which(df.tmp$BLOCK.P == min(df.tmp$BLOCK.P)), "BLOCK.SNP"]
                        out       = c(out, new.lead)
                }
                return(unique(out))
        }
}

# clean output from GCTA cojo
clean.cojo = function(df, prefix = F) {
        df$LeadSNP   = gsub(pattern = ".jma.cojo*", replacement = "", df$LeadSNP)
        if(prefix != FALSE) {
        df$LeadSNP   = gsub(pattern = prefix,       replacement = "", df$LeadSNP)
        }
        df = df[order(df$Chr, df$bp), ]
        for(i in unique(df$LeadSNP)) {
                if(!(i %in% df$SecondarySNP)) {
                        # set the minimum SNP to lead SNP
                        min.p = min(df[which(df$LeadSNP == i), "pJ"])
                        df[which(df$LeadSNP == i & df$pJ == min.p), "SecondarySNP"] = i
                }
        }
        df = subset(df, df$LeadSNP != df$SecondarySNP)
        return(df)
}

# clean VEP output
clean.vep = function(df) {
        names(df) = c("SNP", "CHR", "BP", "REA", "CONSEQUENCE", "GENE", "RSID", "SIFT", "POLYPHEN", "AFR_RAF", "AMR_RAF", "EAS_RAF", "EUR_RAF", "SAS_RAF", "AA_RAF", "EA_RAF")
        df$RSID   = gsub(pattern = ",.*", replacement = "", x = df$RSID)
        df$SNP    = paste0("chr", df$CHR, ":", df$BP)
        df$CONSEQUENCE = gsub(pattern = "_VARIANT", replacement = "", x = df$CONSEQUENCE)
        return(df)
}

# merge regions which share conditionally independent SNPs
blocks.revisited = function(df, location = "primary") {
        tmp = overlapping.cojo.snps(df)
        if((tmp == 0)[1]) {
                return(df)
        } else {
                df.merged  = NULL
                block.list = NULL
                for (i in 1:length(tmp)) {
                        if(!(df[which(df$SecondarySNP == tmp[i]), "BLOCK.SNP"][1] %in% block.list)) {
                                old.blocks = df[which(df$SecondarySNP == tmp[i]), "BLOCK.SNP"]
                                df.tmp = df[which(df$BLOCK.SNP %in% old.blocks), ]
                                df.tmp$COJO.SNP       = df.tmp[which(df.tmp[, "P"] == min(df.tmp[, "P"])), "BLOCK.SNP"][1]
                                df.tmp$COJO.BP        = df.tmp[which(df.tmp[, "P"] == min(df.tmp[, "P"])), "BP"][1]
                                df.tmp$COJO.P         = df.tmp[which(df.tmp[, "P"] == min(df.tmp[, "P"])), "P"][1]
                                df.tmp$COJO.SECONDARY = ifelse(df.tmp$SNP == df.tmp$COJO.SNP[1], 0, 1)
                                df.tmp$COJO.PRIMARY   = ifelse(df.tmp$SNP == df.tmp$COJO.SNP[1], 1, 0)
                                df.tmp$COJO.START     = df.tmp[which(df.tmp[, "BLOCK.START"] == min(df.tmp[, "BLOCK.START"])), "BLOCK.START"][1]
                                df.tmp$COJO.END       = df.tmp[which(df.tmp[, "BLOCK.END"] == max(df.tmp[, "BLOCK.END"])), "BLOCK.END"][1]
                                df.tmp$COJO.N         = min(df.tmp[, "BLOCK.N"]) #hmmm put in the block of the new COJO.SNP
                                # remove one of the entries
                                dup.snp = row.names(df.tmp[which((df.tmp$SNP != df.tmp$COJO.SNP) & (df.tmp$SecondarySNP == tmp[i])), ])
                                # don't remove, just set fields to missing
                                df.tmp[which(row.names(df.tmp) == dup.snp), c("SecondarySNP", "bp", "refA", "freq", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se", "pJ", "LD_r")] = rep(NA, 13)
                                df.merged = rbind(df.merged, df.tmp)
                                block.list = c(block.list,  old.blocks)
                        }
                }
                # keep the regions that were correct as is
                df.resto = df[which(!(df$BLOCK.SNP %in% block.list)), ]
                df.resto$COJO.SNP       = df.resto$BLOCK.SNP
                df.resto$COJO.BP        = df.resto$BLOCK.BP
                df.resto$COJO.P         = df.resto$BLOCK.P
                df.resto$COJO.SECONDARY = df.resto$BLOCK.SECONDARY
                df.resto$COJO.PRIMARY   = df.resto$BLOCK.PRIMARY
                df.resto$COJO.START     = df.resto$BLOCK.START
                df.resto$COJO.END       = df.resto$BLOCK.END
                df.resto$COJO.N         = df.resto$BLOCK.N
                # merge and correct the block id's
                out = rbind(df.resto, df.merged)
                out = out[order(out$CHR, out$BP), ]
                for(i in 2:nrow(out)) {
                        if(out$COJO.N[i] != out$COJO.N[i-1]) {
                                out[which(out$COJO.N == out$COJO.N[i]), "COJO.N"] =  out$COJO.N[i - 1] + 1
                        }
                }
                return(out)
        }
}

# LD: which SNPs are stored where
ld.annotate = function(df.ld, df.ref, df.trans = F, df.eur = F, df.sas = F, df.afr = F, df.amr = F, df.eas = F){
   df.ld[, "REF_A"] = ifelse(df.ld$SNP_A %in% df.ref[, "CHRCBP"], 1, 0)
   df.ld[, "REF_B"] = ifelse(df.ld$SNP_B %in% df.ref[, "CHRCBP"], 1, 0)
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
   if(is.data.frame(df.eas)){
      df.ld[, "EAS_A"] = ifelse(df.ld$SNP_A %in% df.eas[, "SNP"], 1, 0)
      df.ld[, "EAS_B"] = ifelse(df.ld$SNP_B %in% df.eas[, "SNP"], 1, 0)
   }
   return(df.ld)
}

# is the SNP novel?
snp.novel = function(df, chr = "CHR", bp = "BP", snp = "BLOCK.SNP", df.ref, df.ld, half.window = 500000, missing.ref.snps = F) {
   out = NULL
   for (iCHR in 1:22) {
   #print(iCHR)
      df.chr   = df[df[, chr] == iCHR, ]
      ref.chr  = df.ref[df.ref[, "CHR"] == iCHR, ]
      if (is.data.frame(missing.ref.snps)){
            missing.chr = subset(missing.ref.snps, missing.ref.snps$CHR == iCHR)
      }
      if (nrow(df.chr) > 0) {
        df.chr$SENTINEL = NA
        for (i in 1:nrow(df.chr)) {
           #print(i)
           if(nrow(ref.chr) == 0) {
              df.chr$SENTINEL[i] = NA
           }
           else {
              for(j in 1:nrow(ref.chr)) {
                 #print(j)
                 if((df.chr[i, "BLOCK.BP"] >= ref.chr[j, bp] - half.window) & (df.chr[i, "BLOCK.BP"] <= ref.chr[j, bp] + half.window)) {
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
             if (is.data.frame(missing.ref.snps)){
                if(nrow(missing.chr) > 0) {
                   for(j in 1:nrow(missing.chr)) {
                      if((df.chr[i, "BLOCK.BP"] >= missing.chr[j, bp] - half.window) & (df.chr[i, "BLOCK.BP"] <= missing.chr[j, bp] + half.window & (is.na(df.chr$SENTINEL[i]) == T))) {
                          df.chr$SENTINEL[i] = missing.chr[j, "CHRCBP"]
                      }
                   }
                }
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

# flag ancestry-restricted SNPs (not in trans-ethnic meta)
ancestry.restricted = function(df.race, df.trans, df.ld, pop = "EUR", chr = "CHR", bp = "BP", snp = "LeadSNP", trans.snp = "LeadSNP", half.window = 500000) {
   ld.tmp    = subset(df.ld, (df.ld$TRANS > 0) & (df.ld[, pop] > 0)) # take only SNPs that compare a EUR SNP against a TRANS SNP
   trans.ref = snp.coordinates(as.data.frame(df.trans[, trans.snp], stringsAsFactors = F))
   out = NULL
   for (iCHR in 1:22) {
      df.chr   = df.race[df.race[, chr] == iCHR, ]       # split comparisons by chromosome
      pd.chr   = trans.ref[trans.ref[, "CHR"] == iCHR, ] # check physical distance (pd)
      ld.chr   = ld.tmp[ld.tmp$CHR_A == iCHR, ]          # check linkage disequilibrium (ld)
      if (nrow(df.chr) > 0) {
        df.chr$TRANS.SNP   = NA
        df.chr$TRANS.EXACT = NA
        df.chr$TRANS.R2    = NA
        df.chr$TRANS.DIST  = NA
        for (i in 1:nrow(df.chr)) {
           if(nrow(pd.chr) + nrow(ld.chr) == 0) {   # if there are no snps in that chromosome in trans or trans.ld then set to missing
              df.chr$TRANS.SNP[i] = NA
           }
           else {
              if(nrow(pd.chr) != 0) {
                 for(j in 1:nrow(pd.chr)) {  # check physical distance +- 0.5 MB
                    #print(j)
                    if((df.chr[i, bp] >= pd.chr[j, bp] - half.window) & (df.chr[i, bp] <= pd.chr[j, bp] + half.window)) {
                       df.chr$TRANS.SNP[i] = pd.chr[j, "CHRCBP"]
                       df.chr$TRANS.DIST[i]  = abs(pd.chr[j, "BP"] - df.chr[i, bp])
                       if(df.chr$TRANS.SNP[i] == df.chr[i, snp]) {
                          df.chr$TRANS.EXACT[i] = 1
                          df.chr$TRANS.R2[i]    = 1
                       next;
                       }
                    }
                 }
              }
              if((nrow(ld.chr) != 0) & (is.na(df.chr$TRANS.EXACT[i]))) {   # now check ld, onlly when TRANS.SNP has not been found
                 for(k in 1:nrow(ld.chr)) { # check fi in LD
                    if (df.chr[i, snp] %in% ld.chr$SNP_A) {
                       # take trans.snp that has the highest r2
                       r.tmp =  ld.chr[which(ld.chr$SNP_A == df.chr[i, snp]), c("SNP_B", "R2", "BP_B")]
                       df.chr$TRANS.SNP[i]   = r.tmp[which(r.tmp$R2 == max(r.tmp$R2)), "SNP_B"]
                       df.chr$TRANS.R2[i]    = r.tmp[which(r.tmp$R2 == max(r.tmp$R2)), "R2"]
                       df.chr$TRANS.DIST[i]  = abs(r.tmp[which(r.tmp$R2 == max(r.tmp$R2)), "BP_B"] -  df.chr[i, bp])
                    }
                    if (df.chr[i, snp] %in% ld.chr$SNP_B) {
                       r.tmp =  ld.chr[which(ld.chr$SNP_B == df.chr[i, snp]), c("SNP_A", "R2", "BP_A")]
                       if((max(r.tmp$R2) > df.chr$TRANS.R2[i] & is.na(df.chr$TRANS.R2[i]) == F) | (is.na(df.chr$TRANS.R2[i]) == T)) { # overwrite if r2 is larger
                          df.chr$TRANS.SNP[i]  = r.tmp[which(r.tmp$R2 == max(r.tmp$R2)), "SNP_A"]
                          df.chr$TRANS.R2[i]   = r.tmp[which(r.tmp$R2 == max(r.tmp$R2)), "R2"]
                          df.chr$TRANS.DIST[i] = abs(r.tmp[which(r.tmp$R2 == max(r.tmp$R2)), "BP_A"] -  df.chr[i, bp])
                       }
                    }
                 }
              }
           }
        }
        out = rbind(out, df.chr)
      }
   }
   out$TRANS.SNP   = ifelse(is.na(out$TRANS.SNP), "-", out$TRANS.SNP)
   out$TRANS.EXACT = ifelse((out$TRANS.SNP != "-") & (is.na(out$TRANS.EXACT) == T), 0, out$TRANS.EXACT)
   out$RACE.ONLY   = ifelse(out$TRANS.SNP == "-", 1, 0)
   return(out)
}

# ancestry tag for transethnic meta-analysis
tag.ancestry = function(ld.trans, df.eur = F, df.sas = F, df.amr = F, df.eas = F, r2.high = 0.8, r2.low = 0.05, exact.high.combined = T) {
        if(data.frame(df.eur)) { ld.eur  = data.frame(SNP = unique(df.eur[, c("SNP")]), EUR = 1, AFR = 0, SAS = 0, AMR = 0, EAS = 0) } else { ld.eur = data.frame(SNP=NULL, EUR=NULL) }
        if(data.frame(df.amr)) { ld.amr  = data.frame(SNP = unique(df.amr[, c("SNP")]), EUR = 0, AFR = 0, SAS = 0, AMR = 1, EAS = 0) }
        if(data.frame(df.sas)) { ld.sas  = data.frame(SNP = unique(df.sas[, c("SNP")]), EUR = 0, AFR = 0, SAS = 1, AMR = 0, EAS = 0) }
        if(data.frame(df.eas)) { ld.eas  = data.frame(SNP = unique(df.eas[, c("SNP")]), EUR = 0, AFR = 0, SAS = 0, AMR = 0, EAS = 1) }
        if(data.frame(df.afr)) { ld.afr  = data.frame(SNP = unique(df.afr[, c("SNP")]), EUR = 0, AFR = 1, SAS = 0, AMR = 0, EAS = 0) }
        # gotta make sure that LD trans always has EAS AMR, even if it means they are all 0
        ld.tmp1 = ld.trans[,c("SNP_A", "EUR", "AFR", "SAS", "AMR", "EAS")]
        ld.tmp2 = ld.trans[,c("SNP_B", "EUR", "AFR", "SAS", "AMR", "EAS")]
        names(ld.tmp1)[1] = names(ld.tmp2)[1] = "SNP"
        ld.tmp = unique(rbind(ld.tmp1, ld.tmp2, ld.eur, ld.amr, ld.sas, ld.eas, ld.afr))
        ld.full = ddply(ld.tmp, "SNP", numcolwise(sum))
        ld.tmp = ld.full[,2:6]
        ld.tmp[ld.tmp > 1] = 1
        ld.full = cbind(ld.full$SNP, ld.tmp)
        names(ld.full)[1] = "SNP"
        ld.full$AncestryTag = paste(ld.full$EUR, ld.full$AFR, ld.full$SAS, ld.full$AMR, ld.full$EAS, sep=".")
        return(ld.full)
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
      out = subset(out, !grepl("\\.",  out$external_gene_name))
      out = subset(out, !grepl("LINC", out$external_gene_name))
      out = subset(out, !grepl("_",    out$external_gene_name))
      if(nrow(out) > 0) {
         df$ENSID[i] = out$ensembl_gene_id
         df$GENE[i] = out$external_gene_name
      }
      else { # if SNP not in gene, extend range
        df$query[i] = paste(gsub("chr", '', df[i, chr]), df[i, bp] - range.1, df[i, bp] + range.1, sep = ":")
        out = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'chromosomal_region', values = df$query[i], mart = mart.hs )
        out = subset(out, !grepl("\\.",  out$external_gene_name))
        out = subset(out, !grepl("LINC", out$external_gene_name))
        out = subset(out, !grepl("_",    out$external_gene_name))
        if(nrow(out) > 0) {
           df$ENSID[i] = "-"
           df$GENE[i]  = paste(unique(unlist(out$external_gene_name)), collapse = ';')
        }
  }
        else { # if not in first range, then try second range
           df$query[i] = paste(gsub("chr", '', df[i, chr]), df[i, bp] - range.2, df[i, bp] + range.2, sep = ":")
           out = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'chromosomal_region', values = df$query[i], mart = mart.hs)
           out = subset(out, !grepl("\\.",  out$external_gene_name))
           out = subset(out, !grepl("LINC", out$external_gene_name))
           out = subset(out, !grepl("_",    out$external_gene_name))
           df$ENSID[i] = "-"
           df$GENE[i]  = paste(unique(unlist(out$external_gene_name)), collapse = ';')
        }
      }
   }
   return(df)
}
