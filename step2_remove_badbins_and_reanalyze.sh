#!/bin/bash

# Written by Will Chronister and Inusah Diallo
# 2019-09-06
#Rivanna stuff -- if you want to run on cluster rather than frontend
module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3
THREADS="4"
PARTITION="standard"
MEM="80GB"
TIME="300" #minutes

#Set up output subdirectories
Fin=`cat outputs_dir.txt`
mkdir $Fin/pngs_badbinsOUT
mkdir $Fin/full_longs_badbinsOUT
mkdir $Fin/sos_CN_mult_badbinsOUT
#echo '#!/bin/bash' > submit.sh

if [ `ls slurm-*out | wc -l` -gt 0 ]
then
  mv `ls slurm-*out | tail -1` ./last_job_slurm_output
  rm slurm-*out
fi

while read directory; do

  Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/Remove_bad_bins_Segment_Plot.R $directory $Fin

done < Dir_Names.txt

# cat submit.sh | sbatch -p $PARTITION -n $THREADS --mem=$MEM --time=$TIME -A cnv_seq_pipeline
#rm submit.sh
#rm Rplots.pdf
