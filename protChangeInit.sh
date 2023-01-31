## Protein Change Bash Script pipeline

## SCC Settings
cd /projectnb2/evolution/zwakefield/proteinChange/pipeline
## Load Modules
module load python3 bedtools tensorflow R/4.1.2

outname="cardiacDiff"
mkdir ./$outname/
cd ./$outname/
mkdir ./prettyAlignments/
pwd
## Run Python Script
  ## Inputs: GTF, Exon files from downstream analysis diExons
  ## Outputs: transcript names + fasta files
python3 ../protChange_pipe.py '/projectnb/evolution/zwakefield/Annotations/hg38_gencode/gencode.v41.annotation.gtf' '/projectnb2/evolution/zwakefield/proteinChange/full_sheet.xlsx' $outname
#/projectnb/evolution/zwakefield/Annotations/new_hg38_ensembl/hg38.ensGene.gtf
echo "py done"
## getfasta
bedtools getfasta -fi /projectnb2/evolution/zwakefield/Annotations/hg38_gencode/GRCh38.p13.genome.fa -bed $outname.bed -fo ./$outname.fa -s
echo "bedtools done"

## Add ORF obtaining py script, then run transToProt.R again

## Run R script
  ## Proteins
  ## Also which aren't in bioMart

Rscript ../transToProt.R "/projectnb2/evolution/zwakefield/proteinChange/pipeline/$outname/$outname.bed" "/projectnb2/evolution/zwakefield/proteinChange/pipeline/$outname/"
echo "R done"
cd ./prettyAlignments/
mkdir ./tex/
mv *.tex ./tex/



## Python script
  ## any unidentifiable transcripts, perform ORF finder -> protein translation
  ## any non annotated exon sequences, perform ORF finder -> protein translation
