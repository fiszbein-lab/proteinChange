#!/usr/bin/env Rscript
library(methods)
args = commandArgs(trailingOnly=TRUE)
bed_path <- args[1]
out_dir <- args[2]

library(ensembldb, verbose = FALSE)
library(EnsDb.Hsapiens.v86, verbose = FALSE)
library(msa, verbose = FALSE)
library(pracma, verbose = FALSE)
cd(out_dir)
# bed <- read.delim("/projectnb2/evolution/zwakefield/proteinChange/pipeline/cardiacDiff2/cardiacDiff2.bed", header = F)
bed <- read.delim(bed_path, header = F)

possT <- bed$V4

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
cd(paste(out_dir, "prettyAlignments/", sep=""))
for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
  if (protCode[i] == "None" | protCode[i+1] == "None") {
    protC <- c(protC, "None")
    protAlign[[i]] <- "None"
    pMatch <- c(pMatch, 0)
  } else if (protCode[i] == protCode[i+1]) {
    protC <- c(protC, "Same")
    protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
    pMatch <- c(pMatch, 1)
  } else {
    protC <- c(protC, "Different")
    protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)
    msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
                   , file = paste(out_dir, "prettyAlignments/", bed$V4[i], "_", bed$V4[i+1], "_prettyAlignment.pdf", sep = ""), output = "pdf")
    pMatch <- c(pMatch, table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/(sum(table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?"))))))
  }
}

print(table(protC))
bed$protein <- protCode
bed$match <- rep(protC, each = 2)
bed$prop <- pMatch
# write.table(bed, "/projectnb2/evolution/zwakefield/proteinChange/pipeline/proteinOut.txt", quote = F, row.names = F, col.names = F, sep = '\t')                     
write.table(bed, paste(out_dir, "proteinOut.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = '\t')                     
