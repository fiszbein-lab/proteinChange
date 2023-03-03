#!/bin/bash -l

#########################################################################
# Name: Zach Wakefield
# Date: Sep 22, 2022
# Description: protein change due to diExons in Cardiac Differentiation
#########################################################################

# Set SCC project to charge
#$ -P hybrids

# Specify hard time limit for the job.
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.
#   You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=5:00:00


# Name job
#$ -N exonProtChange

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Specify email
#$ -M zachpwakefield@gmail.com

# Join error and output streams in one file
#$ -j y

## Protein Change Bash Script pipeline

## SCC Settings
cd /projectnb2/evolution/zwakefield/proteinChange/pipeline
## Load Modules


module load python3 bedtools tensorflow R/4.1.2
outname="cardiacDifferentiation"
mkdir ./$outname/
cd ./$outname/
pheno_names=("d0_d1" "psc_mpc" "psc_cpc" "mpc_cpc")
phenotypes=("/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/stemcellsDE" "/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/psc_mpc" "/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/psc_cpc" "/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/mpc_cpc")
pipe_dir="/projectnb2/evolution/zwakefield/proteinChange/pipeline"
nC=3
nE=3
nO=2
di="50"
exonTypes=("AFE" "ALE")


# for i in ${!phenotypes[@]}
# do
#   module load python3 bedtools tensorflow R/4.1.2
#   mkdir ${pheno_names[i]}_AFE
#   cd ${pheno_names[i]}_AFE
#   Rscript $pipe_dir/extract_diExons.R ${pheno_names[i]} ${phenotypes[i]} $pipe_dir/$outname/${pheno_names[i]}_AFE/ $di "AFE" $nO $nC $nE
#   mkdir ./prettyAlignments/
#   python3 $pipe_dir/protChange_pipe.py '/projectnb/evolution/zwakefield/Annotations/hg38_gencode/hg38_for_df.txt' ${pheno_names[i]}_AFE_$di.xlsx ${pheno_names[i]} "AFE"
#   echo "py done"
#   bedtools getfasta -fi /projectnb2/evolution/zwakefield/Annotations/hg38_gencode/GRCh38.p13.genome.fa -bed ${pheno_names[i]}.bed -fo ./${pheno_names[i]}.fa -s -name
#   echo "bedtools done"
#   Rscript $pipe_dir/transToProt.R "$pipe_dir/$outname/${pheno_names[i]}_AFE/${pheno_names[i]}.bed" "$pipe_dir/$outname/${pheno_names[i]}_AFE/${pheno_names[i]}.fa" "$pipe_dir/$outname/${pheno_names[i]}_AFE/" "$pipe_dir"
#   echo "R done"
#   cd ./prettyAlignments/
#   mkdir ./tex/
#   mv *.tex ./tex/

#   module unload python3
#   module load java/16.0.2 gcc/8.3.0 python3/3.8.10

#   cd ../
#   mkdir ./plots
#   mv *.pdf* ./plots/
 
#   /projectnb2/evolution/zwakefield/tools/my_interproscan/interproscan-5.60-92.0/interproscan.sh -i ./outFast.fa -f tsv -goterms -T $TMPDIR

#   Rscript $pipe_dir/addDomainInfo.R "$pipe_dir/$outname/${pheno_names[i]}_AFE" "${pheno_names[i]}"
#   module unload python3/3.8.10 java/16.0.2 gcc/8.3.0
#   echo "Domain done"
#   cd ../
# done


for i in ${!phenotypes[@]}
do
  for j in ${!exonTypes[@]}
  do
      module load python3 bedtools tensorflow R/4.1.2
      mkdir ${pheno_names[i]}_${exonTypes[j]}
      cd ${pheno_names[i]}_${exonTypes[j]}
      Rscript $pipe_dir/extract_diExons.R ${pheno_names[i]} ${phenotypes[i]} $pipe_dir/$outname/${pheno_names[i]}_${exonTypes[j]}/ $di ${exonTypes[j]} $nO $nC $nE
      mkdir ./prettyAlignments/
      python3 $pipe_dir/protChange_pipe.py '/projectnb/evolution/zwakefield/Annotations/hg38_gencode/hg38_for_df.txt' ${pheno_names[i]}_${exonTypes[j]}_$di.xlsx ${pheno_names[i]} ${exonTypes[j]}
      echo "py done"
      bedtools getfasta -fi /projectnb2/evolution/zwakefield/Annotations/hg38_gencode/GRCh38.p13.genome.fa -bed ${pheno_names[i]}.bed -fo ./${pheno_names[i]}.fa -s -name
      echo "bedtools done"
      Rscript $pipe_dir/transToProt.R "$pipe_dir/$outname/${pheno_names[i]}_${exonTypes[j]}/${pheno_names[i]}.bed" "$pipe_dir/$outname/${pheno_names[i]}_${exonTypes[j]}/${pheno_names[i]}.fa" "$pipe_dir/$outname/${pheno_names[i]}_${exonTypes[j]}/" "$pipe_dir"
      echo "R done"
      cd ./prettyAlignments/
      mkdir ./tex/
      mv *.tex ./tex/

      module unload python3
      module load java/16.0.2 gcc/8.3.0 python3/3.8.10

      cd ../
      mkdir ./plots
      mv *.pdf* ./plots/

      /projectnb2/evolution/zwakefield/tools/my_interproscan/interproscan-5.60-92.0/interproscan.sh -i ./outFast.fa -f tsv -goterms -T $TMPDIR

      Rscript $pipe_dir/addDomainInfo.R "$pipe_dir/$outname/${pheno_names[i]}_${exonTypes[j]}" "${pheno_names[i]}_${exonTypes[j]}"
      module unload python3/3.8.10 java/16.0.2 gcc/8.3.0
      echo "Domain done"
      cd ../
  done
done

