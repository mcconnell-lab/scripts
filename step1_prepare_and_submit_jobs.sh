#!/bin/bash

# Written by Will Chronister and Inusah Diallo
# 2019-09-06
#This script will take FASTQs or BAMs as input and run them through pipeline
#FASTQs have been tested most recently; BAMs might not work right anymore
#Bad bins are not removed in this script but are identified
#Use next script "step2" to remove bad bins and re-segment, re-plot each cell

module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3

# this portion of the code takes the user input parameters and finds directories
# that are then used to create an output folder name
# Parameters can be run name, species, WGA method, and/or individual
# Example: ./step1_prepare_and_submit_jobs.sh human 26yroNeuNpos,26yroNeuNneg
# Multiple criteria for each parameter must be seperated by a comma, e.g. species1,species2 individual1,individual2...ect


param1=$1
param2=$2
param3=$3
param4=$4

# a parameter that is not used must be indicated using 'none' in the command line
theDate=$(date +"%Y")$(date +"%m")$(date +"%d") #prints out the date in the correct format
finalOutput=$theDate"_""output""_"
all=$(echo "$1","$2","$3","$4")
all2=$(echo $all | sed 's/,/_/g')

#Name output directory ("Fin")
Fin="/nv/vol97/mcconnell_lab/cnvpipe/pipeline_results""/"$finalOutput$all2
echo $Fin > outputs_dir.txt
if [[ $1 == "none" ]]
then
    param1_final="_"
else
param1_final=$( echo $param1 | sed 's[,[ [g' )
fi

if [[ $2 == "none" ]]
then
    param2_final="_"
else
param2_final=$( echo $param2 | sed 's[,[ [g' )
fi

if [[ $3 == "none" ]]
then
    param3_final="_"
else
param3_final=$( echo $param3 | sed 's[,[ [g' )
fi

if [[ $4 == "none" ]]
then
    param4_final="_"
else
param4_final=$( echo $param4 | sed 's[,[ [g' )
fi

#param1_spaced=$( echo $param1 | sed 's[,[ [g' )
#param2_spaced=$( echo $param2 | sed 's[,[ [g' )
#param3_spaced=$( echo $param3 | sed 's[,[ [g' )
#param4_spaced=$( echo $param4 | sed 's[,[ [g' )

#Find directories that match criteria used to select directories for analysis
echo "Found Directories: "> Names.txt
for i in $param1_final
do
  if [ $i != "_" ]; then
    p1=$p1$i"_"
  fi
   for j in $param2_final
    do
        if [ $j != "_" ]; then
        p2=$p2$j"_"
        fi
        for k in $param3_final
        do
                if [ $k != "_" ]; then
                p3=$p3$k"_"
                fi
            for g in $param4_final
            do
                find /nv/vol97/mcconnell_lab/cnvpipe/wgs_data/ -maxdepth 1 -mindepth 1 -type d -name "*$i*" > o1.txt
                #echo 'find /nv/vol97/mcconnell_lab/cnvpipe/wgs_data/ -maxdepth 1 -type d -name '"*$i*"' > o1.txt'
                cat o1.txt | grep "$j" > o2.txt
                cat o2.txt | grep "$k" > o3.txt
                cat o3.txt | grep "$g" >> o4.txt
                if [ $g != "_" ]; then
                p4=$p4$g"_"
                fi
            done
        done
    done
done


sort o4.txt >> Names.txt
cat Names.txt
cat Names.txt | grep "/" > Dir_Names.txt
rm o1.txt
rm o2.txt
rm o3.txt
rm o4.txt
rm Names.txt
Dir_num=`cat Dir_Names.txt | wc -l`

# This portion of the code is carries out the actual analysis of the data and storing of the results

#Variables for running on cluster (can be tuned as needed)
THREADS="4"
PARTITION="standard"
MEM="80GB"
TIME="300" #minutes

#Software/tool locations
BWADIR="/nv/vol97/mcconnell_lab/cnvpipe/software/bin"
BEDDIR="/nv/vol97/mcconnell_lab/cnvpipe/software/bedtools-2.17.0/bin"
SAMDIR="/nv/vol97/mcconnell_lab/cnvpipe/software/samtools-1.1"
PICDIR="/nv/vol97/mcconnell_lab/cnvpipe/software/picard-tools-1.105"
STATDIR="/nv/vol97/mcconnell_lab/cnvpipe/software/bin"

#Size of genomic bins/windows
WINDOW_KB="500"

#Python pipeline
PIPELINE="/nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/CNV_pipeline_FASTQorBAM_thru_BIN_COUNTS.py"

#Make output directory and subdirectories
mkdir $Fin
mkdir $Fin/bams
mkdir $Fin/beds
mkdir $Fin/cov_files_badbinsIN
mkdir $Fin/pngs_badbinsIN
mkdir $Fin/sample_stats
mkdir $Fin/sos_CN_mult_badbinsIN
mkdir $Fin/full_longs_badbinsIN
mkdir $Fin/bin_stuff
mkdir $Fin/bic_outputs
mkdir $Fin/segment_data

#Counters initialized
Overall=0
Dir_count=0

while read inputs
do


  ((Dir_count++))

  if [[ $inputs = *"human"* ]]
  then
    echo Species is human
    SPECIES="human"
    GENOME="/nv/vol97/mcconnell_lab/cnvpipe/genomes/hg19/genome.fa"
    NONUNIQSITES="/nv/vol97/mcconnell_lab/cnvpipe/b37.40mer.nonUniqSites.merged.bed"
  elif [[ $inputs = *"mouse"* ]]
  then
    echo Species is mouse
    SPECIES="mouse"
    GENOME="/nv/vol97/mcconnell_lab/cnvpipe/genomes/mm9_new/genome.fa"
    NONUNIQSITES="/nv/vol97/mcconnell_lab/cnvpipe/MM9.b37.40mer.nonUniqSites.merged.bed"
  fi

  #Check if there are no FASTQs or BAMs
  if [ `find $inputs/ -name *_R[1-2]*fastq* | wc -l` -eq 0 ] && [ `find $inputs/ -name *bam | wc -l` -eq 0 ]
  then
    echo Could not find any FASTQ files named in "_R1" or "_R2" style or BAM files
    #If FASTQs do not contain R1/R2 info, use rename command to modify FASTQ filenames into R1/R2 format
  #Otherwise, check if there are >0 FASTQ AND >0 BAM files
  elif [ `find $inputs/ -name *_R[1-2]*fastq* | wc -l` -gt 0 ] && [ `find $inputs/ -name *bam | wc -l` -gt 0 ]
  then
    echo There are both FASTQ and BAM files in the inputs directory. Please use only one input type.
  #Otherwise, check if there are >0 FASTQ files
  elif [ `find $inputs/ -name "*fastq*" | wc -l` -gt 0 ]
  then
    INPUTTYPE="fastq"
    LIFTOVER="no"
    #Determine if FASTQs are paired-end
    if [ `find $inputs/*fastq* | wc -l` -gt `find $inputs/*fastq* | awk -F"_R[1-2]" '{$0=$1}1' | sort -u | wc -l` ]
    then
      echo Inputs are paired-end FASTQs
      #Generate commands for submission to SLURM
      find $inputs/trimmed*fastq* | sort -u | xargs -n 2 | sed -e 's/ /,/g' > fastq_pairs.txt
      cell_count=0
      while read pairedfastq; do
      echo  "python -u $PIPELINE -p -f $pairedfastq -r $GENOME -w $WINDOW_KB -u $NONUNIQSITES -n $THREADS -b $BWADIR -e $BEDDIR -s $SAMDIR -t $STATDIR -c $PICDIR -sp $SPECIES -in $INPUTTYPE -lO $LIFTOVER -en $Fin" >> oneliner_jobs.txt
      ((cell_count++))
      done < fastq_pairs.txt
      rm fastq_pairs.txt
      Cell_num=0
      while read oneliner; do
        ((Cell_num++))
        ((Overall++))
        
        echo '#!/bin/bash' > single_job_submit.sh
      # echo 'module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3' >> single_job_submit.sh	
        echo $oneliner >> single_job_submit.sh
        #Check if last job (last directory, last cell)
        if [ $Dir_count == $Dir_num ] && [ $cell_count ==  $Cell_num ]
        then
          echo 'Last job; adding lines to .sh script'
          #Adding line to test if number of bin_CN_values.txt files is less than number of jobs
          echo 'while [ `ls '"$Fin"'/bams/*sorted_rmdup.bam.bai | wc -l` -lt '"$Overall"' ]; do echo Wait for it...; sleep 20; done' >> single_job_submit.sh
          #When all jobs have completed, job will proceed to next lines: moving/removing appended files

          #echo 'mv $inputs2/All_Conf_Scores.txt $Fin/sample_stats/' >> single_job_submit.sh
          #echo 'mv $inputs2/All_FlagStats.txt $Fin/sample_stats/' >> single_job_submit.sh
          #echo 'mv $inputs2/bin_CN_values.txt $Fin/bin_stuff/' >> single_job_submit.sh
          #echo 'mv $inputs2/bin_raw_counts.txt $Fin/bin_stuff/' >> single_job_submit.sh
          #echo 'mv $inputs2/BICs_and_other_info_badbinsIN.txt $Fin/bic_outputs/' >> single_job_submit.sh
          #echo 'mv $inputs2/SoS_CN_multiplier_info.txt $Fin/sos_CN_mult_badbinsIN/' >> single_job_submit.sh
          #echo 'rm $inputs2/tempres.txt' >> single_job_submit.sh
          #Add line to run bad bin detection
    
          echo 'Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/Bad_bin_finder.R '"$Fin"' '"$inputs"'' >> single_job_submit.sh
          #echo 'while [ `ls $Fin/bin_stuff/ | wc -l` -lt 3 ]; do echo Waiting for Bad_bin_finder...; sleep 20; done' >> single_job_submit.sh
          #echo 'rm Rplots.pdf' >> single_job_submit.sh
          echo 'rm '"$Fin"'/*_bin_CN_values.txt' >> single_job_submit.sh
          echo 'rm '"$Fin"'/*_bin_raw_counts.txt' >> single_job_submit.sh

          echo 'rm '"$Fin"'/*_segment_data.txt' >> single_job_submit.sh

        fi
        chmod 777 single_job_submit.sh
        cat single_job_submit.sh | sbatch -p $PARTITION -n $THREADS --mem=$MEM --time=$TIME -A cnv_seq_pipeline
        #Keep track of each job submitted
        date >> all_jobs_log.txt
        cat single_job_submit.sh >> all_jobs_log.txt
        rm single_job_submit.sh
      done < oneliner_jobs.txt
      cat oneliner_jobs.txt >> oneline_jobs.txt
      rm oneliner_jobs.txt
    #Determine if FASTQs are single-end
    elif [ `find $inputs/*fastq* | wc -l` -eq `find $inputs/*fastq* | awk -F"_R[1-2]" '{$0=$1}1' | wc -l` ]
    then
    echo Inputs are single-end FASTQs
      #Generate commands for submission to SLURM
      cell_count=0
      find $inputs/trimmed*fastq* | sort -u | xargs -n 1 > fastq_singles.txt
      while read singlefastq; do
        ((cell_count++))
        echo "python -u $PIPELINE -f $singlefastq -r $GENOME -w $WINDOW_KB -u $NONUNIQSITES -n $THREADS -b $BWADIR -e $BEDDIR -s $SAMDIR -t $STATDIR -c $PICDIR -sp $SPECIES -in $INPUTTYPE -lO $LIFTOVER -en $Fin " >> oneliner_jobs.txt
      done < fastq_singles.txt
      rm fastq_singles.txt
      Cell_num=0    
      while read oneliner; do
        ((Cell_num++))
        ((Overall++))
        echo '#!/bin/bash' > single_job_submit.sh
      # echo 'module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3' >> single_job_submit.sh
        echo $oneliner >> single_job_submit.sh
        #Check if last job (last directory, last cell)
        if [ $Dir_count == $Dir_num ] && [ $cell_count ==  $Cell_num ]
        then
          echo 'Last job; adding lines to .sh script'
          #Adding line to test if number of bin_CN_values.txt files is less than number of jobs
          echo 'while [ `ls '"$Fin"'/bams/*sorted_rmdup.bam.bai | wc -l` -lt '"$Overall"' ]; do echo Wait for it...; sleep 20; done' >> single_job_submit.sh
          #When all jobs have completed, job will proceed to next lines: moving/removing appended files

          #echo 'mv $inputs2/All_Conf_Scores.txt $Fin/sample_stats/' >> single_job_submit.sh
          #echo 'mv $inputs2/All_FlagStats.txt $Fin/sample_stats/' >> single_job_submit.sh
  #        echo 'mv $inputs2/bin_CN_values.txt $Fin/bin_stuff/' >> single_job_submit.sh
  #        echo 'mv $inputs2/bin_raw_counts.txt $Fin/bin_stuff/' >> single_job_submit.sh
          #echo 'mv $inputs2/BICs_and_other_info_badbinsIN.txt $Fin/bic_outputs/' >> single_job_submit.sh
          #echo 'mv $inputs2/SoS_CN_multiplier_info.txt $Fin/sos_CN_mult_badbinsIN/' >> single_job_submit.sh
          #echo 'rm $inputs2/tempres.txt' >> single_job_submit.sh
          #Add line to run bad bin detection

          echo 'Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/Bad_bin_finder.R '"$Fin"' '"$inputs"'' >> single_job_submit.sh
  #        echo 'while [ `ls $Fin/bin_stuff/ | wc -l` -lt 3 ]; do echo Waiting for Bad_bin_finder...; sleep 20; done' >> single_job_submit.sh
  #        echo 'rm Rplots.pdf' >> single_job_submit.sh
          echo 'rm '"$Fin"'/*_bin_CN_values.txt' >> single_job_submit.sh
          echo 'rm '"$Fin"'/*_bin_raw_counts.txt' >> single_job_submit.sh

          echo 'rm '"$Fin"'/*_segment_data.txt' >> single_job_submit.sh
        fi
        chmod 777 single_job_submit.sh
        cat single_job_submit.sh | sbatch -p $PARTITION -n $THREADS --mem=$MEM --time=$TIME -A cnv_seq_pipeline
        #Keep track of each job submitted
        date >> all_jobs_log.txt
        cat single_job_submit.sh >> all_jobs_log.txt
        rm single_job_submit.sh
      done < oneliner_jobs.txt
      cat oneliner_jobs.txt >> oneline_jobs.txt
      rm oneliner_jobs.txt

    fi
  #Otherwise, check if there are >0 BAM files
  elif [ `find $inputs/ -name *bam | wc -l` -gt 0 ]
  then
    echo Inputs are BAMs
    INPUTTYPE="bam"
    #Check if liftOver (changes coordinates from mm10-mm9) required
    if [[ $inputs = *"10X_mouse"* ]] && [[ `find $inputs/ -name *bam | wc -l` -gt 0 ]]
    then
      #We can infer that if the data is from 10X, is mouse, and is BAM,
      #it was aligned to mm10 and must be liftOver-ed
      echo Input BAMs will be converted to BED and liftOver-ed from mm10 to mm9 before proceeding with analysis.
      LIFTOVER="yes"
    fi
    #Generate commands for submission to SLURM
    cell_count=0
    find $inputs/*bam | sort -u | xargs -n 1 > bam_singles.txt
    while read singlebam; do
      ((cell_count++))
      echo "python -u $PIPELINE -f $singlebam -r $GENOME -w $WINDOW_KB -u $NONUNIQSITES -n $THREADS -b $BWADIR -e $BEDDIR -s $SAMDIR -t $STATDIR -c $PICDIR -sp $SPECIES -in $INPUTTYPE -lO $LIFTOVER -en $Fin" >> oneliner_jobs.txt
    done < bam_singles.txt
    rm bam_singles.txt
    Cell_num=0
    while read oneliner; do
        ((Cell_num++))
        ((Overall++))
        echo '#!/bin/bash' > single_job_submit.sh
      # echo 'module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3' >> single_job_submit.sh
        echo $oneliner >> single_job_submit.sh
        if [ $Dir_count == $Dir_num ] && [ cell_count ==  $Cell_num ]
        then
          echo 'Last job; adding lines to .sh script'
          #Adding line to test if number of bin_CN_values.txt files is less than number of jobs
          echo 'while [ `ls '"$Fin"'/bams/*sorted_rmdup.bam.bai | wc -l` -lt '"$Overall"' ]; do echo Wait for it...; sleep 20; done' >> single_job_submit.sh
          #When all jobs have completed, job will proceed to next lines: moving/removing appended files

          #echo 'mv $inputs2/All_Conf_Scores.txt $Fin/sample_stats/' >> single_job_submit.sh
          #echo 'mv $inputs2/All_FlagStats.txt $Fin/sample_stats/' >> single_job_submit.sh
  #        echo 'mv $inputs2/bin_CN_values.txt $Fin/bin_stuff/' >> single_job_submit.sh
  #        echo 'mv $inputs2/bin_raw_counts.txt $Fin/bin_stuff/' >> single_job_submit.sh
          #echo 'mv $inputs2/BICs_and_other_info_badbinsIN.txt $Fin/bic_outputs/' >> single_job_submit.sh
          #echo 'mv $inputs2/SoS_CN_multiplier_info.txt $Fin/sos_CN_mult_badbinsIN/' >> single_job_submit.sh
          #echo 'rm $inputs2/tempres.txt' >> single_job_submit.sh
          #Add line to run bad bin detection

          echo 'Rscript /nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/Bad_bin_finder.R '"$Fin"' '"$inputs"'' >> single_job_submit.sh
  #        echo 'while [ `ls $Fin/bin_stuff/ | wc -l` -lt 3 ]; do echo Waiting for Bad_bin_finder...; sleep 20; done' >> single_job_submit.sh
  #        echo 'rm Rplots.pdf' >> single_job_submit.sh
          echo 'rm '"$Fin"'/*_bin_CN_values.txt' >> single_job_submit.sh
          echo 'rm '"$Fin"'/*_bin_raw_counts.txt' >> single_job_submit.sh

          echo 'rm '"$Fin"'/*_segment_data.txt' >> single_job_submit.sh

      fi
          chmod 777 single_job_submit.sh
      cat single_job_submit.sh | sbatch -p $PARTITION -n $THREADS --mem=$MEM --time=$TIME -A cnv_seq_pipeline
        date >> all_jobs_log.txt
        cat single_job_submit.sh >> all_jobs_log.txt
        rm single_job_submit.sh
      done < oneliner_jobs.txt
  cat oneliner_jobs.txt >> oneline_jobs.txt
    rm oneliner_jobs.txt

  fi



done < Dir_Names.txt
#rm oneliner_jobs.txt
#chmod 775 submit_jobs.sh

#find ./rawfastqs/*fastq | sort -u

#    echo '"echo '\"asdfa $TIME,$TIME\"' | asdf"'
#          "echo "asdfa 1234,1234" | asdf"

#  echo "echo \"asdfa $TIME,$TIME\""
#        echo "asdfa 1234,1234"

