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
bed_path <- args[1]
fa_path <- args[2]
out_dir <- args[3]
pipe_dir <- args[4]



# cd(out_dir)
# bed <- read.delim("/projectnb2/evolution/zwakefield/proteinChange/pipeline/cardiacDiff/cardiacDiff.bed", header = F)
bed <- read.delim(bed_path, header = F)
gcpc_trans <- read_lines(paste(pipe_dir, '/gencode.v43.pc_translations.fa', sep = ""))
index <- grep(">", gcpc_trans)
endInd <- c(index[-1], (length(gcpc_trans)+1))
new_trans <- list()
for (i in 1:length(index)) {
  new_trans[[i]] <- c(gcpc_trans[index[i]])
  new_trans[[i]] <- c(new_trans[[i]], paste(gcpc_trans[(index[i]+1):(endInd[i]-1)], collapse = ""))
}
c_trans <- unlist(new_trans)
# "/projectnb2/evolution/zwakefield/proteinChange/pipeline/cardiacDiff/cardiacDiff.bed"
bed <- read.delim(bed_path, header = F)
trans <- unlist(lapply(strsplit(unique(bed$V4), "#"), "[[", 1))
# bed <- read.delim(bed_path, header = F)

possT <- unlist(lapply(strsplit(bed$V4, "#"), "[[", 1))


protCode <- unlist(lapply(trans, function(x) if (sum(grepl(x, c_trans)) > 0) {
  c_trans[(grep(x, c_trans)+1)[1]]}
  else{
    "none"
  }))
proBed <- data.frame(id = unique(bed$V4), prot = protCode) %>% separate(id, c("transcript", "id"), "#") %>% separate("id", c("gene", "chr", "start", "stop"), ";")
protAlign <- list()
protC <- c()
pMatch <- c()
alignType <- c()
cate <- c()
# out_dir <- "/projectnb2/evolution/zwakefield/proteinChange/pipeline/cardiacDiff/"
cd(paste(out_dir, "prettyAlignments/", sep=""))
for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
  if (protCode[i] == "none" | protCode[i+1] == "none") {
    if (protCode[i] == "none" & protCode[i+1] != "none") {
      protC <- c(protC, "nonPC", "PC")
      protAlign[[i]] <- "none"
      protAlign[[i]] <- "onePC"
      alignType <- c(alignType, "onePC")
    } else if(protCode[i] != "none" & protCode[i+1] == "none") {
      protC <- c(protC, "PC", "nonPC")
      protAlign[[i]] <- "onePC"
      alignType <- c(alignType, "onePC")
    } else {
      protC <- c(protC, "nonPC", "nonPC")
      protAlign[[i]] <- "none"
      alignType <- c(alignType, "noPC")
    }
    pMatch <- c(pMatch, 0)
  } else if (protCode[i] == protCode[i+1]) {
    protC <- c(protC, "Same", "Same")
    protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
    pMatch <- c(pMatch, 1)
    alignType <- c(alignType, "Match")
  } else {
    protC <- c(protC, "Different", "Different")
    protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)
    
    minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
    pMatch <- c(pMatch, table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/min(nchar(protCode[i]), nchar(protCode[i+1])))
    # table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/(sum(table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?"))))))
    if (nchar(paste(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > (.1*minPc)], collapse = "")) > .2*minPc)  {
      alignType <- c(alignType, "PartialMatch")
      msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
                     , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_prettyAlignment.pdf", sep = ""), output = "pdf")
    } else {
      alignType <- c(alignType, "FrameShift")
      msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
                     , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_prettyAlignment.pdf", sep = ""), output = "pdf")
    }
  }
}

print(table(protC))
# bed$protein <- protCode
proBed$match <- protC
proBed$prop <- rep(pMatch, each = 2)
# Filled Density Plot

# (gdf <- ggplot(data.frame(dens = as.numeric(pMatch), type = alignType), aes(x = dens, fill = type)) + 
#     geom_histogram(aes(y=..count../sum(..count..)), colour = 1,
#                    bins = 20) + geom_density(aes(y=.05*..count..), color = 'black', fill = "coral2", bw = .1, alpha = .3) + 
#     scale_fill_manual(values=c("azure4", "#E69F00", "#56B4E9", "darkblue", "deeppink4")) +
#     theme_classic() + xlab("Alignment Score") + ylab("Count"))

(gdf <- ggplot(data.frame(dens = as.numeric(pMatch), type = alignType), aes(x = dens, fill = type)) + 
    geom_histogram(aes(y=..count../sum(..count..)), colour = 1,
                   bins = 20) + geom_density(aes(y=.0005*..count..), color = 'black', fill = "coral2", bw = .1, alpha = .3) + 
    scale_fill_manual(values=c("azure4", "#E69F00", "#56B4E9", "darkblue", "deeppink4")) +
    theme_classic() + xlab("Alignment Score") + ylab("Fraction"))
pdf(file = paste(out_dir, "alignScores.pdf", sep = ""))
print(gdf)
dev.off()

proBed$matchType <- rep(alignType, each = 2)
print(table(alignType))
# write.table(bed, "/projectnb2/evolution/zwakefield/proteinChange/pipeline/proteinOut.txt", quote = F, row.names = F, col.names = F, sep = '\t')                     
write.table(proBed, paste(out_dir, "proteinOut.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = '\t')  

proFast <- c()
for (i in seq(from=1,to=(length(proBed[,1])-1), by=2)) {
  if (proBed$prot[i] != "none" & proBed$prot[i+1] != "none") {
    proFast <- c(proFast, paste(">", proBed$transcript[i], "#", proBed$gene[i], ";", proBed$chr[i], ";", proBed$start[i], "-", proBed$stop[i], sep = ""), 
                 proBed$prot[i], paste(">", proBed$transcript[i+1], "#", proBed$gene[i+1], ";", proBed$chr[i+1], ";", proBed$start[i+1], "-", proBed$stop[i+1], sep = ""), 
                 proBed$prot[i+1])
  }
}
write_lines(proFast, paste(out_dir, "outFast.fa", sep = ""))
