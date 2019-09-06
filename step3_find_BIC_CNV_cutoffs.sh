#!/bin/bash

# Written by Will Chronister and Inusah Diallo
# 2019-09-06
# Step3 involves finding BIC cutoffs (ensures high quality cells)
# and new CNV thresholds (ensures bona fide CNVs) for the dataset using mixed Gaussian models
# In practice, this step may need to be done manually to ensure the Gaussian models do, in fact,
# find the best cutoffs for the given dataset
# Step 3 also utilizes previously defined CNV cutoffs known as lenient (1.34, 2.60) and stringent (1.14, 2.80)
# (see Chronister et al. 2019 for more info on lenient/stringent CNV cutoffs)

module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3
#Set number of Gaussian distributions to be fit to BIC data
#1 works best if uniform normal distribution; 2 works best if long tail of bad cells
GAUSS="2"
#Rivanna stuff (can be run on frontend too)
THREADS="4"
PARTITION="standard"
MEM="80GB"
TIME="300"
Fin=`cat outputs_dir.txt`
#echo '#!/bin/bash' > submit.sh
mkdir $Fin/cnv_results
mkdir $Fin/cnv_results/new_thresholds
mkdir $Fin/cnv_results/lenient_thresholds
mkdir $Fin/cnv_results/stringent_thresholds

Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/BIC_n_CNV_cutoff_finder.R $Fin $GAUSS Dir_Names.txt

# cat submit.sh | sbatch -p $PARTITION -n $THREADS --mem=$MEM --time=$TIME -A cnv_seq_pipeline
#rm Rplots.pdf
