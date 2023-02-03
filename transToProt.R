#!/usr/bin/env Rscript
library(methods)
args = commandArgs(trailingOnly=TRUE)
bed_path <- args[1]
out_dir <- args[2]

library(tidyverse)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(msa)
library(pracma)
cd(out_dir)
# bed <- read.delim("/projectnb2/evolution/zwakefield/proteinChange/pipeline/cardiacDiff/cardiacDiff.bed", header = F)
bed <- read.delim(bed_path, header = F)

possT <- unlist(lapply(strsplit(bed$V4, split = ";"), "[[", 1))

edb <- EnsDb.Hsapiens.v86
prts <- proteins(edb,
                 return.type = "AAStringSet")


transToProt <- data.frame(prts@elementMetadata@listData)
protCode <- unlist(lapply(possT, function(y) {if (length(as.data.frame(prts[prts@ranges@NAMES == transToProt[transToProt$tx_id == y,]$protein_id])$x) == 0) {
  "None"
} else {
  as.data.frame(prts[prts@ranges@NAMES == transToProt[transToProt$tx_id == y,]$protein_id])$x
}
}))
# protCode

protAlign <- list()
protC <- c()
pMatch <- c()
alignType <- c()
# out_dir <- "/projectnb2/evolution/zwakefield/proteinChange/pipeline/cardiacDiff/"
cd(paste(out_dir, "prettyAlignments/", sep=""))
for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
  if (protCode[i] == "None" | protCode[i+1] == "None") {
    protC <- c(protC, "None")
    protAlign[[i]] <- "None"
    pMatch <- c(pMatch, 0)
    alignType <- c(alignType, "None")
  } else if (protCode[i] == protCode[i+1]) {
    protC <- c(protC, "Same")
    protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
    pMatch <- c(pMatch, 1)
    alignType <- c(alignType, "Match")
  } else {
    protC <- c(protC, "Different")
    protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)
    msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
                   , file = paste(out_dir, "prettyAlignments/", possT[i], "_", possT[i+1], "_prettyAlignment.pdf", sep = ""), output = "pdf")
    pMatch <- c(pMatch, table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/(sum(table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?"))))))
    minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
    if (bed$V6[i] == "+") {
      if (nchar(paste(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > 0], collapse = "")) > .5*minPc)  {
        alignType <- c(alignType, "PartialMatch")
      } else {
        alignType <- c(alignType, "FrameShift")
      }
    } else {
      if (nchar(paste(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > 0], collapse = "")) > .5*minPc) {
        alignType <- c(alignType, "PartialMatch")
      } else {
        alignType <- c(alignType, "FrameShift")
      }
    }
  }
}

print(table(protC))
bed$protein <- protCode
bed$match <- rep(protC, each = 2)
bed$prop <- rep(pMatch, each = 2)
# Filled Density Plot
gdf <- ggplot(data.frame(dens = pMatch), aes(x = dens)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "lightblue", bins = 7) +
  geom_density(color="darkblue") + theme_classic()
pdf(file = paste(out_dir, "alignScores.pdf", sep = ""))
print(gdf)
dev.off()

bed$matchType <- rep(alignType, each = 2)
print(table(alignType))
# write.table(bed, "/projectnb2/evolution/zwakefield/proteinChange/pipeline/proteinOut.txt", quote = F, row.names = F, col.names = F, sep = '\t')                     
write.table(bed, paste(out_dir, "proteinOut.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = '\t')                     
