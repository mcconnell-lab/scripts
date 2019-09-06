#Contributors: Will Chronister, Inusah Diallo
#2019-09-06

#Creates bed files for use in bedtools genomecov later in step 4
#bed file = chr[tab]start[tab]stop
#Carried out for all three sets of CNV thresholds (lenient, stringent, new)

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
#pkgTest("MASS")
library(dplyr)
options(stringsAsFactors = F)
args <- commandArgs(trailing=T)

#directory <- "/nv/vol97/mcconnell_lab/cnvpipe/vandenBos2016_strandseq_human_AD01"
directory <- args[1]
fin <- args[2]
#rundir <- basename(directory)
#run <- strsplit(rundir,split="_")[[1]][1]
#wga <- strsplit(rundir,split="_")[[1]][2]
#species <- strsplit(rundir,split="_")[[1]][3]
#indivname <- strsplit(rundir,split="_")[[1]][4]
dirtable <- read.table(directory)

dirtable_vector <- dirtable[,1]
indiv_names <- c()
for(j in dirtable_vector){
  indiv_names <- c(indiv_names, strsplit(basename(j), split="_")[[1]][4])
}
indiv_names <- unique(indiv_names)

for(thresh in c("new_thresholds","lenient_thresholds","stringent_thresholds")){
  #cnvs <- read.table("vandenBos2016_strandseq_human_AD01/cnvs_strandseq_bicfilter.txt",header=T)
  cnvs <- read.table(paste(fin,"/cnv_results/",thresh,"/cnvs_bicfilter.txt",sep=""),header=T)
  
  cnvs_sorted <- cnvs %>% arrange(chrom, loc.start, loc.end)
  cnvs_sorted$chrom <- gsub("^","chr",cnvs_sorted$chrom)
  dels_sorted <- cnvs_sorted %>% filter(cnv_type=="del")
  dups_sorted <- cnvs_sorted %>% filter(cnv_type=="dup")
  
  cnvs_sorted_bed <- cnvs_sorted %>% select(chrom, loc.start, loc.end)
  dels_sorted_bed <- dels_sorted %>% select(chrom, loc.start, loc.end)
  dups_sorted_bed <- dups_sorted %>% select(chrom, loc.start, loc.end)
  cnv_indivs <- unique(cnvs_sorted$Individual)
  
  if(length(cnv_indivs) > 0){
    if(nrow(cnvs_sorted_bed)>0){
    write.table(cnvs_sorted_bed,file=paste(fin,"/cnv_results/",thresh,"/all_cnvs_sorted.bed",sep=""),quote=F, sep="\t", row.names=F, col.names=F)

    }
    if(nrow(dels_sorted_bed)>0){
     write.table(dels_sorted_bed,file=paste(fin,"/cnv_results/",thresh,"/all_dels_sorted.bed",sep=""),quote=F, sep="\t", row.names=F, col.names=F)
    }
    if(nrow(dups_sorted_bed)>0){
    write.table(dups_sorted_bed,file=paste(fin,"/cnv_results/",thresh,"/all_dups_sorted.bed",sep=""),quote=F, sep="\t", row.names=F, col.names=F)
    }

    for(i in 1:length(cnv_indivs)){
    indiv_cnvs_sorted <- cnvs_sorted %>% filter(Individual==cnv_indivs[i])
    indiv_cnvs_sorted_bed <- indiv_cnvs_sorted %>% select(chrom, loc.start, loc.end)
    write.table(indiv_cnvs_sorted_bed,file=paste(fin,"/cnv_results/",thresh,"/",cnv_indivs[i],"_cnvs_sorted.bed",sep=""),quote=F, sep="\t", row.names=F, col.names=F)
    indiv_dels_sorted_bed <- indiv_cnvs_sorted %>% 
      filter(cnv_type=="del") %>% select(chrom, loc.start, loc.end)
    if(nrow(indiv_dels_sorted_bed) > 0){
      write.table(indiv_dels_sorted_bed,file=paste(fin,"/cnv_results/",thresh,"/",cnv_indivs[i],"_dels_sorted.bed",sep=""),quote=F, sep="\t", row.names=F, col.names=F)
    }
    indiv_dups_sorted_bed <- indiv_cnvs_sorted %>% 
      filter(cnv_type=="dup") %>% select(chrom, loc.start, loc.end)
    if(nrow(indiv_dups_sorted_bed) > 0){
      write.table(indiv_dups_sorted_bed,file=paste(fin,"/cnv_results/",thresh,"/",cnv_indivs[i],"_dups_sorted.bed",sep=""),quote=F, sep="\t", row.names=F, col.names=F)
    }
    }
  }else{
    print("no cnvs detected")
  }
  
  
  
  }







