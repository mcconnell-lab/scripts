#Contributors: Will Chronister and Inusah Diallo
#2019-09-06

#Calculates summary statistics for each run+WGA+species+individual analyzed,
#across each CNV threshold (lenient, stringent, new)


#Our (slightly modified) version of a function by Sacha Epskamp for checking for the presence of a package and installing if absent.
#http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
pkgTestBC <- function(x){
        if(!require(x,character.only=TRUE)){
                source("http://bioconductor.org/biocLite.R")
                biocLite(x)
                if(!require(x,character.only=TRUE)) stop("Package not found")
        }
}

pkgTest <- function(x){
        if(!require(x,character.only=TRUE)){
                install.packages(x)
                if(!require(x,character.only=TRUE)) stop("Package not found")
        }
}
pkgTest("dplyr")
pkgTest("ggplot2")
pkgTest("tidyr")
pkgTest("MASS")

library(dplyr)
library(ggplot2)

options(stringsAsFactors = F)

#setwd("~/uva-PoolJuly2018_iGenomxRiptide_mouse_p53KOmouse3/")
args <- commandArgs(trailing=T)

directorys<- args[1]
fin <- args[2]
#all_BICs_and_other_info_badbinsOUT_passfail
#cnvs_iGenomxRiptide_bicfilter.txt
dirtable <- read.table(directorys)

dirtable_vector <- dirtable[,1]

dir_name <- dirtable_vector[1]
directory <- paste(dir_name)
rundir <- basename(directory)
run <- strsplit(rundir,split="_")[[1]][1]
wga <- strsplit(rundir,split="_")[[1]][2]
species <- strsplit(rundir,split="_")[[1]][3]
indivname <- strsplit(rundir,split="_")[[1]][4]

for(thresh in c("new_thresholds","lenient_thresholds","stringent_thresholds")){
  dirthreshcnv <- paste(fin,"/cnv_results/",thresh,"/",sep="")

  #BIC
  bics <- read.table(paste(fin,"/bic_outputs/all_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''),sep="\t",header=T) 
  summary_bic <- bics %>%
    group_by(Run,WGA,Species,Individual,BIC_cutoff) %>%
    summarize(Num_cells=n(),
              Num_pass_BIC=sum(BIC_passfail=="Pass_BIC"),
              Pct_pass_BIC=100*(sum(BIC_passfail=="Pass_BIC")/n()),
              Avg_BIC=mean(BIC),
              Avg_segments=mean(Segments),
              Avg_phi=mean(Phi),
              Avg_reads=mean(Total_reads))

  cnvs <- read.table(paste(dirthreshcnv,"cnvs_bicfilter.txt",sep=""), sep="\t", header=T)
  summary_cnv <- cnvs %>%
    group_by(Run,WGA,Species,Individual,BIC_cutoff) %>%
    summarize(Num_CNV_cells=n_distinct(Sample),
              Num_CNVs=n(),
              CNVs_per_CNV_cell=n()/n_distinct(Sample),
              Avg_CNV_size_Mb=mean(loc.end-loc.start)/1000000,
              Avg_CNV_bins=mean(num.mark))

  summary_merged <- merge(summary_bic,summary_cnv,by=c("Run","WGA","Species","Individual","BIC_cutoff"),all.x=T)
  summary_merged <- summary_merged %>%
    mutate(Pct_CNV_cells=100*Num_CNV_cells/Num_pass_BIC)

  summary_merged <- summary_merged[,c(1:5,9,6:8,10:13,18,14:17)]


  write.table(summary_merged, file = paste(fin,"/all_individuals_",thresh, "_bic_cnv_summary.txt",sep=""),
              sep="\t",col.names=T,row.names=F,quote=F)
}
