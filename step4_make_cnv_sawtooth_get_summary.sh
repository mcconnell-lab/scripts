#!/bin/bash

# Written by Will Chronister and Inusah Diallo
# 2019-09-06

#"Sawtooth" plot is simply the cumulative CNV plot showing locations of all detected CNVs, which resembles a saw
#Plots cumulative CNVs both combined and separated by del/dup type, for each of three thresholds (lenient, stringent, new)
#Also calculates summary statistics about individuals

module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3
#Rivanna stuff (can be run on cluster if desired)
THREADS="4"
PARTITION="standard"
MEM="80GB"
TIME="300"
Fin=`cat outputs_dir.txt`

#Set up CNVs for bedtools to do genomecov later
Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/organize_CNVs_for_bedtools_genomecov.R Dir_Names.txt $Fin

if [[ $Fin = *"human"* ]]
then
  echo Species is human
  GENOME="/nv/vol97/mcconnell_lab/cnvpipe/software/bedtools-2.17.0/genomes/human.hg19.ONLY_chr1_to_chr22.genome"
elif [[ $Fin = *"mouse"* ]]
then
  echo Species is mouse
  GENOME="/nv/vol97/mcconnell_lab/cnvpipe/software/bedtools-2.17.0/genomes/mouse.mm9.ONLY_chr1_to_chr19.genome"
fi

find $Fin/cnv_results/*/*bed | sort -u | xargs -n 1 > $Fin/cnv_results/cnv_beds.txt
if [[ `cat $Fin/cnv_results/cnv_beds.txt | grep '/' | wc -l` -gt 0 ]]
then
  while read bed; do
    outname=$(echo "$bed" | sed -e 's/_sorted.bed/_genomecov.bed/')
    echo Generating "$outname"
    bedtools genomecov -bga -i $bed -g $GENOME > $outname
  done < $Fin/cnv_results/cnv_beds.txt
  rm $Fin/cnv_results/cnv_beds.txt

  Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/plot_genomecov_beds.R Dir_Names.txt $Fin
  #Beds no longer needed at this point
  rm $Fin/cnv_results/new_thresholds/*.bed
  rm $Fin/cnv_results/lenient_thresholds/*.bed
  rm $Fin/cnv_results/stringent_thresholds/*.bed
  rm Rplots.pdf
  Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/compute_final_stats.R Dir_Names.txt $Fin

  echo "Step4 of the pipeline is complete"
else
  echo No CNVs were found
fi
#rm Dir_Names.txt
#rm outputs_dir.txt

## run by entering "ls *bed | xargs -n 1 sh coverage.sh"
#chrname=$(echo "$1" | sed -e 's/cnvs.bed/cnvs_chr.bed/')
#echo "Creating $chrname"
#sed 's/^/chr/' "$1" | tr -d '\r' > $chrname
##cat "$chrname" | tr -d '\r' > "$chrname"
##dos2unix "$chrname"
#finalname=$(echo "$chrname" | sed -e 's/chr.bed/genomecov.bed/')
#echo "Creating $finalname"
#bedtools genomecov -bga -i $chrname -g /sfs/nfs/blue/wdc8rf/software/bedtools-2.17.0/genomes/human.hg19.genome > $finalname
