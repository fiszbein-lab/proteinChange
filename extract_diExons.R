#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(methods)
  library(ggpubr)
  library(xlsx)
  library(ggsignif)
  library(tidyverse)
  library(ComplexHeatmap)
  library(tidyverse)
  library(purrr)
  library(foreach)
  library(doParallel)
  library(hypeR)
})
args = commandArgs(trailingOnly=TRUE)
outname <- args[1]#"mpc_cmc"
file_path <- args[2]#"/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/mpc_cpc"
out_dir <- args[3]#"/projectnb2/evolution/zwakefield/proteinChange/pipeline/"
di_baseline <- as.numeric(args[4])/100#.5
exon_type <- args[5]#"AFE"
numOutliers_acc <- as.numeric(args[6])#2
nC <- as.numeric(args[7])#3
nE <- as.numeric(args[8])#3



# outname <- "mpc_cmc"
# file_path <- "/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/mpc_cpc"
# out_dir <- "/projectnb2/evolution/zwakefield/proteinChange/pipeline/"
# nC <- 3
# nE <- 3
# di_baseline <- .5
# exon_type <- "AFE"
# numOutliers_acc <- 2

GO.bp <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
GO.cc <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
GO.mf <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
genesetsC2 <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
genesetsH <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
geneset_BIOCARTA <- msigdb_gsets("Homo sapiens", "C2", "CP:BIOCARTA", clean=TRUE)


lfc <- function(de_df, numCont, numExp, exon_type) {
  de_df <- de_df[de_df$type == exon_type,]
  samps <- colnames(de_df)
  # unlist(lapply(strsplit(colnames(de_df)[-c(1:10, (10+numCont+numExp+1):dim(de_df)[2])], split = "[.]"), "[[", 5))
  
  lfc <- list()
  col <- list()
  

  for (i in 1:length(de_df$gene)) {
    outlier <- which(unlist(lapply(samps, function(x) grepl(x, de_df$outlier[i]))))+10
    
    cont <- c(11:(11-1+numCont))
    cont <- cont[!(cont %in% outlier)]
    exp <- c((numCont+10+1):(11-1+numCont+numExp))
    exp <- exp[!(exp %in% outlier)]
    
    
    
    if (length(exp) == 0) {
      lfc_t <- -20.0
    } 
    if (length(cont) == 0) {
      lfc_t <- 20.0
    }
    if (length(exp) == 0 & length(cont) == 0) {
      lfc_t <- 0.0
    }
    if (length(exp) != 0 & length(cont) != 0) {
      lfc_t <- log((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.00001)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.00001))
    }
    
    
    lfc[[i]] <- lfc_t
  }
  
  de_df$lfc <- unlist(lfc)
  de_df <- de_df[de_df$p_value >= 0,]
  col <- list()
  for (i in 1:length(de_df$gene)) {
    col[[i]] <- "#A7A9AC"#natparks.pals("Acadia", 15)[7]
      if (de_df$lfc[i] <= -1.0 & de_df$p_value[i] < .001) {
        col[[i]] <- "#FE9234" #natparks.pals("Acadia", 15)[15]
      }
    if (de_df$lfc[i] >= 1.0 & de_df$p_value[i] < .001) {
      col[[i]] <- "#00A79D"# natparks.pals("Acadia", 15)[1]
    }
    
  }
  de_df$col <- unlist(col)
  
  p.scDE <- de_df
  p.scDE$numOutliers <- str_count(p.scDE$outlier, "[.]")
  return (de = p.scDE)
}


make_lfcPlot <- function(lfc_df, num_thresh = 10) {
  hg38.conv <- readRDS("/projectnb2/evolution/zwakefield/proteinChange/pipeline/hg38_geneRef_conv.RDS")
  lfc_df$hgnc <- unlist(lapply(lfc_df$gene, function(x) {ifelse(unlist(lapply(strsplit(x, split = "[.]"), "[[", 1)) %in% hg38.conv$ens, 
                                                                unique(hg38.conv$hgnc[hg38.conv$ens == unlist(lapply(strsplit(x, split = "[.]"), "[[", 1))]), x)
    
  }))
  hg38.conv$hgnc[hg38.conv$ens %in% unlist(lapply(strsplit(lfc_df$gene, split = "[.]"), "[[", 1))]
  
  
  lab_thresh <- lfc_df %>% arrange(desc(abs(lfc)), p_value)
  lfc_df$hgnc[-log(lfc_df$p_value) <= -log(lab_thresh$p_value[num_thresh]) | abs(lfc_df$lfc) < abs(lab_thresh$lfc[num_thresh])]  <- ""
  
  (deExons <- ggplot(lfc_df, aes(x = lfc, y = -log(p_value), color = col, label = hgnc)) + geom_point(aes(shape = type), size = 1.5, color = lfc_df$col) +
      theme_classic() + ylab("-Log Adj P Value") + xlab("Log2FC")
    # +coord_cartesian(xlim = c(-100, 20))
    + geom_text(hjust=.2, vjust=0, size = 3)
  )
  print("lfc done!")
  return(deExons)
}


hypeRenr <- function(lfc_df, tot_df, fdr_in = .05, pval_in = .05, lfc_in = 1, sd = FALSE) {
  
  hg38.conv <- readRDS("/projectnb2/evolution/zwakefield/proteinChange/pipeline/hg38_geneRef_conv.RDS")
  
  # lfc_df <- lfc_out
  # tot_df <- read.delim(file_path, sep = " ")
  # fdr_in = .05
  # pval_in = .05
  # lfc_in = 1.5
  # sd = FALSE
  
  upPre <- hg38.conv$hgnc[hg38.conv$ens %in% unlist(lapply(strsplit(lfc_df$gene[lfc_df$p_value <= pval_in & lfc_df$lfc > lfc_in], split = "[.]"), "[[", 1))]
  downPre <- hg38.conv$hgnc[hg38.conv$ens %in% unlist(lapply(strsplit(lfc_df$gene[lfc_df$p_value <= pval_in & lfc_df$lfc < -lfc_in], split = "[.]"), "[[", 1))]
  print("up and down made!")
  
  if (sd == TRUE) {
    enrList <- list(up = unique(setdiff(upPre, downPre)),
                    down = unique(setdiff(downPre, upPre)),
                    shared = unique(intersect(upPre, downPre))
    )
  } else {
    
    enrList <- list(up = unique(upPre),
                    down = unique(downPre),
                    shared = unique(intersect(upPre, downPre))
    )
  }
  fdrUse <- fdr_in
  backgroundGenes <- hg38.conv$hgnc[hg38.conv$ens %in% unlist(lapply(strsplit(unique(c(tot_df$gene)), "[.]"), "[[", 1))]
  print("prep done!")
  hypGOcc <- hypeR(enrList, GO.cc, background = backgroundGenes, test="hypergeometric")
  hypGOmf <- hypeR(enrList, GO.mf, background = backgroundGenes, test="hypergeometric")
  hypGObp <- hypeR(enrList, GO.bp, background = backgroundGenes, test="hypergeometric")
  hypKEGG <- hypeR(enrList, genesetsC2, background = backgroundGenes, test="hypergeometric")
  hypHall <- hypeR(enrList, genesetsH, background = backgroundGenes, test="hypergeometric")
  hypBIO <- hypeR(enrList, geneset_BIOCARTA, background = backgroundGenes, test="hypergeometric")
  print("hyp objects made...")
  
  GOcc.dots <- hyp_dots(hypGOcc, fdr = fdrUse, title = "GO Cell Comp Enrichment", merge = TRUE)
  GOmf.dots <- hyp_dots(hypGOmf, fdr = fdrUse, title = "GO Mol Fxn Enrichment", merge = TRUE)
  GObp.dots <- hyp_dots(hypGObp, fdr = fdrUse, title = "GO Biol Proc Enrichment", merge = TRUE)
  KEGG.dots <- hyp_dots(hypKEGG, fdr = fdrUse, title = "KEGG Enrichment", merge = TRUE)
  Hall.dots <- hyp_dots(hypHall, fdr = fdrUse, title = "Hallmark Enrichment", merge = TRUE)
  BIO.dots <- hyp_dots(hypBIO, fdr = fdrUse, title = "BIOCARTA Enrichment", merge = TRUE)
  print("dots made...")
  gocc <- GOcc.dots + theme_bw()
  gomf <- GOmf.dots + theme_bw()
  gobp <- GObp.dots + theme_bw()
  kegg <- KEGG.dots + theme_bw()
  hall <- Hall.dots + theme_bw()
  bio <- BIO.dots + theme_bw()
  
  return(list(gocc = gocc,
              gomf = gomf,
              gobp = gobp,
              kegg = kegg,
              hall = hall,
              bio = bio
  ))
}

e.2p <- function(df, d_psi_cutoff, exon_type, numOut = numOutliers_acc) {
  df <- df[df$X.001 == "True" & df$numOutliers < numOut & abs(df$delta_PSI) >= d_psi_cutoff & df$type == exon_type,]
  df <- df[df$gene %in% names((table(df$gene) > 1)[table(df$gene) > 1]),]
  keep.g <- list()
  for (i in 1:length(unique(df$gene))) {
    psi_split <- df$delta_PSI[df$gene == unique(df$gene)[i]]
    if (TRUE %in% unique(psi_split > 0) & FALSE %in% unique(psi_split > 0)) {
      keep.g[[length(keep.g)+1]] <- unique(df$gene)[i]
    }
  }
  keep.g <- unlist(keep.g)
  df <- df[df$gene %in% keep.g,]
  return(df)
}


makeExcel <- function(df, filename, out, ex_type = exon_type) {
  o.df <- df %>% dplyr::select(gene, exon, p_value, colnames(df)[11:16], outlier, delta_PSI, lfc)
  colnames(o.df)[4:9] <- unlist(lapply(strsplit(colnames(o.df)[4:9], split = "out."), "[[", 2))
  if (ex_type == "AFE") {
    o.df$FE_ann <- meta$FE[paste(meta$gene, meta$loc, sep=";") %in% paste(o.df$gene, o.df$exon, sep=";")] > 0
  } else {
    o.df$LE_ann <- meta$LE[paste(meta$gene, meta$loc, sep=";") %in% paste(o.df$gene, o.df$exon, sep=";")] > 0
    
  }
  
  # write.table(o.df, file = paste("/projectnb2/evolution/zwakefield/proteinChange/pipeline/", filename, ".csv", sep=""), quote = FALSE, row.names = FALSE, sep = ",")
  write.xlsx(o.df, paste(out_dir, filename, ".xlsx", sep=""), row.names = FALSE)
             
  }



df <- read.delim(file_path, sep = " ")
df.l <- lfc(df, numCont = nC, numExp = nE, exon_type = exon_type)
write.table(df.l, paste(out_dir, outname, "_lfc_", exon_type, ".csv", sep = ""), col.names = T, quote = F, sep = ",")
lfcPlot <- make_lfcPlot(df.l)

pdf(file = paste(out_dir, outname, "_volcanoPlot_", exon_type, ".pdf", sep = ""))
print(lfcPlot)
dev.off()

enrPlots <- hypeRenr(df.l, df)

pdf(file = paste(out_dir, outname, "_enrichmentPlots_", exon_type, ".pdf", sep = ""))
print(enrPlots)
dev.off()

p.as.df <- e.2p(df.l, di_baseline, exon_type)


meta <- read.delim("/projectnb2/evolution/zwakefield/proteinChange/pipeline/hg38_gen_meta.txt", header = FALSE)
meta$loc <- unlist(lapply(strsplit(meta$V4, split = "[;]"), "[[", 1))
meta$gene <- unlist(lapply(strsplit(meta$V4, split = "[;]"), "[[", 2))

meta$FE <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(meta$V4, split = "[;]"), "[[", 4)), split = "[:]"), "[[", 2)))
meta$IE <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(meta$V4, split = "[;]"), "[[", 5)), split = "[:]"), "[[", 2)))
meta$LE <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(meta$V4, split = "[;]"), "[[", 6)), split = "[:]"), "[[", 2)))
meta$SE <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(meta$V4, split = "[;]"), "[[", 7)), split = "[:]"), "[[", 2)))

makeExcel(p.as.df, filename=paste(outname, "_", exon_type, "_", 100*di_baseline, sep = ""))



# one.exon.commdf <- data.frame(gene = unlist(lapply(strsplit(names(sort(table(c(unique(paste(p.as.d0_d1$gene, p.as.d0_d1$exon, sep=";")),
# unique(paste(p.as.psc_mpc$gene, p.as.psc_mpc$exon, sep=";")),
# unique(paste(p.as.psc_cpc$gene, p.as.psc_cpc$exon, sep=";")),
# unique(paste(p.as.mpc_cpc$gene,p.as.mpc_cpc$exon, sep=";")))), decreasing = TRUE)[sort(table(c(unique(paste(p.as.d0_d1$gene, p.as.d0_d1$exon, sep=";")),
# unique(paste(p.as.psc_mpc$gene, p.as.psc_mpc$exon, sep=";")),
# unique(paste(p.as.psc_cpc$gene, p.as.psc_cpc$exon, sep=";")),
# unique(paste(p.as.mpc_cpc$gene,p.as.mpc_cpc$exon, sep=";")))), decreasing = TRUE)>0]), split = "[;]"), "[[", 1)),
# exon = unlist(lapply(strsplit(names(sort(table(c(unique(paste(p.as.d0_d1$gene, p.as.d0_d1$exon, sep=";")),
# unique(paste(p.as.psc_mpc$gene, p.as.psc_mpc$exon, sep=";")),
# unique(paste(p.as.psc_cpc$gene, p.as.psc_cpc$exon, sep=";")),
# unique(paste(p.as.mpc_cpc$gene,p.as.mpc_cpc$exon, sep=";")))), decreasing = TRUE)[sort(table(c(unique(paste(p.as.d0_d1$gene, p.as.d0_d1$exon, sep=";")),
# unique(paste(p.as.psc_mpc$gene, p.as.psc_mpc$exon, sep=";")),
# unique(paste(p.as.psc_cpc$gene, p.as.psc_cpc$exon, sep=";")),
# unique(paste(p.as.mpc_cpc$gene,p.as.mpc_cpc$exon, sep=";")))), decreasing = TRUE)>0]), split = "[;]"), "[[", 2)), data.frame(t(data.frame(lapply(names(sort(table(c(unique(paste(p.as.d0_d1$gene, p.as.d0_d1$exon, sep=";")),
# unique(paste(p.as.psc_mpc$gene, p.as.psc_mpc$exon, sep=";")),
# unique(paste(p.as.psc_cpc$gene, p.as.psc_cpc$exon, sep=";")),
# unique(paste(p.as.mpc_cpc$gene,p.as.mpc_cpc$exon, sep=";")))), decreasing = TRUE)[sort(table(c(unique(paste(p.as.d0_d1$gene, p.as.d0_d1$exon, sep=";")),
# unique(paste(p.as.psc_mpc$gene, p.as.psc_mpc$exon, sep=";")),
# unique(paste(p.as.psc_cpc$gene, p.as.psc_cpc$exon, sep=";")),
# unique(paste(p.as.mpc_cpc$gene,p.as.mpc_cpc$exon, sep=";")))), decreasing = TRUE)>0]),
# 
# 
# function(x) unlist(lapply(list(p.as.d0_d1,
# p.as.psc_mpc,
# p.as.psc_cpc,
# p.as.mpc_cpc), function(y) if (x %in% unique(paste(y$gene, y$exon, sep = ";"))) {y$delta_PSI[unique(paste(y$gene, y$exon, sep=";")) %in% x]} else {FALSE}))
# ))), row.names = NULL))
# colnames(one.exon.commdf)[3:6] <- c("d0_d1", "psc_mpc", "psc_cmc", "mpc_cmc")

# one.exon.commdf
# one.exon.commdf$LE_ann <- meta$LE[paste(meta$gene, meta$loc, sep=";") %in% paste(one.exon.commdf$gene, one.exon.commdf$exon, sep=";")] > 0
# write.xlsx(one.exon.commdf, "/projectnb2/evolution/zwakefield/proteinChange/table_ALE.xlsx", row.names = F)




