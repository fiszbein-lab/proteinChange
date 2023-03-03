#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methods)
  library(tidyverse)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(msa)
  library(pracma)
  library(stringr)
  library(ORFik)               
  library(GenomicFeatures)
  library(data.table)
  library(Biostrings) 
})

args = commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
res_name <- args[2]



scores <- read.delim(paste(out_dir, "/proteinOut.txt", sep = ""), header = F)
scores$id <- paste(paste(paste(scores$V1, scores$V2, sep = "#"), scores$V3, sep = ";"), paste(scores$V4, scores$V5, sep = "-"), sep = ";")

fa <- read_lines(paste(out_dir, "/outFast.fa", sep = ""))
nFa <- fa[grep(">", fa)]
gN <- gsub(">", "", nFa)
geneL <- unlist(lapply(strsplit(unlist(lapply(strsplit(fa[grep(">", fa)], split = ";"), "[[", 1)), split = "#"), "[[", 2))

ips <- read_tsv(paste(out_dir, "/outFast.fa.tsv", sep = ""), col_names = F) %>% separate(X13, into = c("X13", "X14"), sep = '\t', fill = 'right')
ips$X14[is.na(ips$X14)] <- rep("-", length(which(is.na(ips$X14))))
ips <- ips %>% dplyr::filter(X4 %in% c("Pfam", "ProSitePatterns", "CDD", "PANTHER"))

protInf <- list()
protDom <- list()
protInt <- list()
upS <- seq(from=1,to=(length(gN)-1), by=2)
for (i in upS) {
  protInf[[i]] <- ips$X6[ips$X1 == gN[i]]
  protInf[[i+1]] <- ips$X6[ips$X1 == gN[i+1]]
  protInt[[i]] <- intersect(protInf[[i]], protInf[[i+1]])
  protInt[[i+1]] <- intersect(protInf[[i]], protInf[[i+1]])
  protDom[[i]] <- setdiff(protInf[[i]], protInf[[i+1]])
  protDom[[i+1]] <- setdiff(protInf[[i+1]], protInf[[i]])
  
  
  protInf[[i]] <- paste(protInf[[i]], collapse = ';')
  protInf[[i+1]] <- paste(protInf[[i+1]], collapse = ';')
  protInt[[i]] <- paste(protInt[[i]], collapse = ';')
  protInt[[i+1]] <- paste(protInt[[i+1]], collapse = ';')
  protDom[[i]] <- paste(protDom[[i]], collapse = ';')
  protDom[[i+1]] <- paste(protDom[[i+1]], collapse = ';')
}

protInf.o <- unlist(protInf)
protIntersect.o <- unlist(protInt)
protSetDiff.o <- unlist(protDom)


protInf.o[protInf.o == ""] <- "none"
protIntersect.o[protIntersect.o == ""] <- "none"
protSetDiff.o[protSetDiff.o == ""] <- "none"


finOut <- data.frame(exon = gN,
                     protein = unlist(lapply(gN, function(x) scores$V6[scores$id == x])),
                     alignScore = unlist(lapply(gN, function(x) scores$V8[scores$id == x])),
                     alignType = unlist(lapply(gN, function(x) scores$V9[scores$id == x])),
                     protInfor = protInf.o,
                     protIntersect = protIntersect.o,
                     protSetDiff = protSetDiff.o)
write.table(finOut, paste(out_dir, "/", res_name, "_completeOut.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = '\t')  
print("Pipe Complete!")


