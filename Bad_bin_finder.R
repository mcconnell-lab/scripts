#Written by Will Chronister and Inusah Diallo
#Last updated 2019-09-06

#This script is run after all cells have completed main CNV analysis during step1.
#It gathers all bin value files (*bin_CN_values.txt), combines them, converts to log2 values,
#and identifies bad/outlier bins via Tukey's Outlier Method.
#Autosomal bad bins are identified from both males and female data;
#XY bad bins identified from separated male and female data.
#Bad bins are removed from each cell and reanalyzed for CNVs in step2

#Originally, this script would methodically incorporate results from other analyses to inform its outlier bin detection,
#but this feature was determined to be beneficial only in rare situations.
#Code for doing so still exists below, but is skipped.
#Now the script only uses the cells being actively analyzed for its outlier detection.

#setwd("~/r_script_tests/RunU999-Pool101_PicoPLEX_human_5401/")
#directory <- "C:/Users/Will/lab/r_script_tests/RunU999-Pool101_PicoPLEX_human_5401/"
args <- commandArgs(trailing=T)
directory <- args[1]
fin <- args[1]
indir <- args[2]
#print(directory)
rundir <- basename(indir)
run <- strsplit(rundir,split="_")[[1]][1]
wga <- strsplit(rundir,split="_")[[1]][2]
species <- strsplit(rundir,split="_")[[1]][3]
indivname <- strsplit(rundir,split="_")[[1]][4]
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

library(dplyr)
library(ggplot2)
library(tidyr)

#turn off scientific notation
options(scipen = 999)
#turn off stringsAsFactors
options(stringsAsFactors = FALSE)

#Compile all shorttable files in current directory
compiled_segment_data <- data.frame()
for(i in list.files(path=directory,pattern="_segment_data.txt", full.names = T)){
  temp_segment_data <- read.table(i, header=T, sep="\t")
  compiled_segment_data <- rbind(compiled_segment_data, temp_segment_data)
}
write.table(compiled_segment_data, file =paste(fin,"/segment_data/all_segment_data_badbinsIN.txt",sep='') , quote=F, col.names = T,sep = "\t")

#Compile all bad bin files in current directory
compiled_bin_CN_values <- data.frame()
for(i in list.files(path=directory, pattern="_bin_CN_values.txt", full.names = T)){
  temp_bin_cn_values <- read.table(i, header=T, sep = "\t")
  compiled_bin_CN_values <- rbind(compiled_bin_CN_values, temp_bin_cn_values)
}
write.table(compiled_bin_CN_values, file = paste(fin,"/bin_stuff/bin_CN_values.txt",sep=''), quote = F, sep="\t",row.names = F, col.names = T)

compiled_bin_raw_counts <- data.frame()
for(i in list.files(path=directory, pattern="_bin_raw_counts.txt", full.names = T)){
  temp_bin_raw_counts <- read.table(i, header=T, sep = "\t")
  compiled_bin_raw_counts <- rbind(compiled_bin_raw_counts, temp_bin_raw_counts)
}
write.table(compiled_bin_raw_counts, file = paste(fin,"/bin_stuff/bin_raw_counts.txt",sep=''), quote = F, sep="\t",row.names = F, col.names = T)

#setwd("~/vandenBos2016_strandseq_human_AD03/")
#Check if file was written correctly and then proceed
if(file.exists(paste(fin,"/bin_stuff/bin_CN_values.txt",sep=''))){
  binvals_current_dir <- read.table(paste(fin,"/bin_stuff/bin_CN_values.txt",sep=''),header = T,sep = "\t")
} else {
  print("Could not find bin_CN_values.txt in bin_stuff directory")
}

#Double check that every cell produced a bin_CN_values file (equal to num. of bed files)
#And check that compilation of bin_CN_values has same num. rows as num. bed files
if(length(list.files(path=directory, pattern="_bin_CN_values.txt")) == length(list.files(path=paste(fin,"/beds"), pattern="_fullread.bed")) &
   nrow(binvals_current_dir) == length(list.files(path=paste(fin,"/beds",sep=''), pattern="_fullread.bed"))){
  print("All cells are accounted for")
} else {
  print("One or more cells are missing bin_CN_values")
}

###Infer sex from median CN of X chromosome
#Get bin location info from any long file
longfile <- read.table(list.files(paste(fin,"/full_longs_badbinsIN/",sep=''),full.names = T)[1], header = T, sep="\t")
if(species == "human"){
  chrX <- 23
  chrX_binvals <- binvals_current_dir[,which(longfile$Chrom == chrX)]
  chrX_median_binval <- median(as.vector(as.matrix(chrX_binvals)))
} else if (species == "mouse") {
  chrX <- 20
  chrX_binvals <- binvals_current_dir[,which(longfile$Chrom == chrX)]
  chrX_median_binval <- median(as.vector(as.matrix(chrX_binvals)))
}
#Check if median chrX bin is under 1.5 (male); otherwise, assume female
#if(chrX_median_binval < 1.5){
 # sex <- "male"
#} else {
#  sex <- "female"
#}

#Refine this***
dirs <- list.dirs(path=paste(fin,"/..",sep=''),recursive = F)
#print(dirs)
subdirs <- c()
for(u in 1:length(dirs)){
  #print(paste(dirs[u],"/bin_stuff",sep=""))
  subdir <- list.dirs(path=paste(dirs[u],"/bin_stuff",sep=""))
  if(length(subdir) > 0){
    subdirs <- c(subdirs, subdir)
  }
}

subdirs_final <- c()
for(v in 1:length(subdirs)){
  subdir_split <- strsplit(subdirs[v],split = "/\\.\\./")[[1]]
  if(dirname(subdir_split[2]) != basename(subdir_split[1])){
    subdirs_final <- c(subdirs_final, subdirs[v])
  }
}

#print(subdirs_final)

#dirs <- grep(".*_.*_.*_.*\\/.*_.*_.*_.*\\/",
#             list.dirs(path=paste(directory,"..",sep='')), value = T)
#dirs <- dirs[grep("bin_stuff",dirs)]

#Keep only directories from same WGA method and species
#subdirs_final <- grep(glob2rx(paste("*../*",wga,"_",species,"*",sep=""),trim.head = F, trim.tail = F), subdirs_final,value = T)
#subdirs_final <- subdirs_final[grep(wga,subdirs_final)]
#Keep only directories from same species
#subdirs_final <- subdirs_final[grep(species,subdirs_final)]
subdirs_final <- c()

if(length(subdirs_final) > 455455){
  print("Eligible directories for additional bin values:")
  print(subdirs_final)
} else {
  print("No other eligible directories for additional bin values found")
}
# this is to prevent the script from looking for additional bin values 
subdirs_final <- c()
binvals <- binvals_current_dir
binvals_opposite_sex <- data.frame()
if(length(subdirs_final) > 436555645){
  for(i in 1:length(subdirs_final)){
    binfile <- list.files(subdirs_final[i],"bin_CN_values.txt",recursive = T,full.names = T)
    if(length(binfile) > 0){
      binvals_from_file <- read.table(binfile, header=T, sep = "\t")
      chrX_binvals_from_file <- binvals_from_file[,which(longfile$Chrom == chrX)]
      chrX_median_binval_from_file <- median(as.vector(as.matrix(chrX_binvals_from_file)))
      if(chrX_median_binval_from_file < 1.5){
        sex_file <- "male"
      } else {
        sex_file <- "female"
      }
      print("Incorporating bin values from:")
      print(sex_file)
      print(binfile)
      #If sex is same, add to binvals
      ##Adjust this later so everything goes to same dataframe, sorted out later
      if(sex == sex_file){
        binvals <- rbind(binvals, binvals_from_file)
      } else {
        binvals_opposite_sex <- rbind(binvals_opposite_sex, binvals_from_file)
      }
    }
  }
}

number_of_cells <- nrow(binvals)


#Reformat binvals
binvals <- as.data.frame(t(binvals))
rownames(binvals) <- seq(1:nrow(binvals))
#print(binvals[nrow(binvals),])
colnames(binvals) <- binvals[nrow(binvals),]
#Remove row with cell names
binvals <- binvals[-nrow(binvals),]
binvals$Bin <- seq(1:nrow(binvals))

#binvals <- gather(binvals,key="Cell",value="Bin_value",colnames(binvals)[1:number_of_cells])
binvals <- gather(binvals,key="Cell",value="Bin_value",1:number_of_cells)

autosome_binvals <- binvals %>% filter(Bin %in% which(longfile$Chrom < chrX))


#Reformat binvals_opposite_sex
print("Bin values available from opposite sex:")
print(nrow(binvals_opposite_sex) > 0)
if(nrow(binvals_opposite_sex) > 0){
  number_of_oppsex_cells <- nrow(binvals_opposite_sex)

  binvals_opposite_sex <- as.data.frame(t(binvals_opposite_sex))
  rownames(binvals_opposite_sex) <- seq(1:nrow(binvals_opposite_sex))
  colnames(binvals_opposite_sex) <- binvals_opposite_sex[nrow(binvals_opposite_sex),]
  #Remove row with cell names
  binvals_opposite_sex <- binvals_opposite_sex[-nrow(binvals_opposite_sex),]
  binvals_opposite_sex$Bin <- seq(1:nrow(binvals_opposite_sex))

#  binvals_opposite_sex <- gather(binvals_opposite_sex,key="Cell",value="Bin_value",colnames(binvals_opposite_sex)[1:number_of_oppsex_cells])
  binvals_opposite_sex <- gather(binvals_opposite_sex,key="Cell",value="Bin_value",1:number_of_oppsex_cells)

  #combine autosomal binvals with main
  autosome_binvals_opposite_sex <- binvals_opposite_sex %>%
    filter(Bin %in% which(longfile$Chrom < chrX))
  autosome_binvals <- rbind(autosome_binvals, autosome_binvals_opposite_sex)
}



autosome_medbinvals <- autosome_binvals %>%
  group_by(Bin) %>%
  summarize(Median_bin_value=median(as.numeric(Bin_value)))
autosome_medbinvals <- autosome_medbinvals %>%
  mutate(log2_median_bin_value=log2(Median_bin_value))
autosome_medbinvals <- cbind("Chrom"=longfile$Chrom[which(longfile$Chrom < chrX)],
                             "Start"=longfile$Start[which(longfile$Chrom < chrX)],
                             "End"=longfile$End[which(longfile$Chrom < chrX)],
                             autosome_medbinvals)

# medbinvals <- cbind("Chrom"=longfile$Chrom, "Start"=longfile$Start,
#                     "End"=longfile$End, medbinvals)
# medbinvals <- medbinvals %>%
#   mutate(log2_median_bin_value=log2(Median_bin_value))



#Use Tukey Outlier Test to find outlying bins
#Autosomes first
autosome_log2_25pct <- as.numeric(quantile(autosome_medbinvals$log2_median_bin_value)[2])
autosome_log2_75pct <- as.numeric(quantile(autosome_medbinvals$log2_median_bin_value)[4])
IQR_autosome_log2 <- IQR(autosome_medbinvals$log2_median_bin_value)

outlier_bins_autosome <- autosome_medbinvals$Bin[autosome_medbinvals$log2_median_bin_value < (autosome_log2_25pct - 1.5*IQR_autosome_log2)]
outlier_bins_autosome <- sort(c(outlier_bins_autosome,
                          autosome_medbinvals$Bin[autosome_medbinvals$log2_median_bin_value > (autosome_log2_75pct + 1.5*IQR_autosome_log2)]))


autosome_chrom_breaks <- cumsum(rle(autosome_medbinvals$Chrom)$lengths)

#Plot bad bins and thresholds
print("Plotting autosomal outlier bins")
asdf <- ggplot(autosome_medbinvals,
       mapping=aes(x=Bin, y=Median_bin_value)) +
  geom_vline(size=0.5,col="gray88",xintercept=autosome_chrom_breaks[-length(autosome_chrom_breaks)]) +
  geom_point(alpha=0.55, size=2) +
  geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(autosome_log2_25pct - 1.5*IQR_autosome_log2)) +
  geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(autosome_log2_75pct + 1.5*IQR_autosome_log2)) +
  theme_bw(base_size=24) +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.subtitle = element_text(size=12, hjust = 0.5)) +
  scale_x_continuous(expand = c(.01,.01)) +
  ylab("Median bin value") +
  ggtitle(paste("Autosomes,", length(outlier_bins_autosome),"outlier bins", sep=" ")) +
  labs(subtitle=basename(directory)) +
  ylim(0,5)
ggsave(filename = paste(fin,"/bin_stuff/median_bin_values_autosomes.png",sep=''), plot=asdf, height=4, width=12)
bictable <- read.table(paste(fin,"BICs_and_other_info_badbinsIN.txt",sep="/"),header=T)
#If sex is male, get cutoffs for X and Y at same time
male_sample_names <- c("N/A")
female_sample_names <- c("N/A")
outlier_bins_male <- c("N/A")
outlier_bins_female<- c("N/A")


if(sum(bictable$Sex=="male")>0){
  XY_binvals <- binvals %>% filter(Bin %in% which(longfile$Chrom >= chrX),Cell %in% bictable$Sample[which(bictable$Sex=="male")])
  male_sample_names <- c()
  male_sample_names <- (bictable$Sample[which(bictable$Sex=="male")])
  XY_medbinvals <- XY_binvals %>%
    group_by(Bin) %>%
    summarize(Median_bin_value=median(as.numeric(Bin_value)))
  XY_medbinvals <- XY_medbinvals %>%
    mutate(log2_median_bin_value=log2(Median_bin_value))
  XY_medbinvals <- cbind("Chrom"=longfile$Chrom[which(longfile$Chrom >= chrX)],
                               "Start"=longfile$Start[which(longfile$Chrom >= chrX)],
                               "End"=longfile$End[which(longfile$Chrom >= chrX)],
                               XY_medbinvals)
  XY_chrom_breaks <- cumsum(rle(XY_medbinvals$Chrom)$lengths) + last(autosome_chrom_breaks)


  XY_log2_25pct <- as.numeric(quantile(XY_medbinvals$log2_median_bin_value)[2])
  XY_log2_75pct <- as.numeric(quantile(XY_medbinvals$log2_median_bin_value)[4])
  IQR_XY_log2 <- IQR(XY_medbinvals$log2_median_bin_value)

  outlier_bins_XY <- XY_medbinvals$Bin[XY_medbinvals$log2_median_bin_value < (XY_log2_25pct - 1.5*IQR_XY_log2)]
  outlier_bins_XY <- sort(c(outlier_bins_XY,
                            XY_medbinvals$Bin[XY_medbinvals$log2_median_bin_value > (XY_log2_75pct + 1.5*IQR_XY_log2)]))
  print("Plotting XY outlier bins")
  sdfg <- ggplot(XY_medbinvals,
         mapping=aes(x=Bin, y=Median_bin_value)) +
    geom_vline(size=0.5,col="gray88",xintercept=XY_chrom_breaks[-length(XY_chrom_breaks)]) +
    geom_point(alpha=0.55, size=2) +
    geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(XY_log2_25pct - 1.5*IQR_XY_log2)) +
    geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(XY_log2_75pct + 1.5*IQR_XY_log2)) +
    theme_bw(base_size=18) +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          plot.subtitle = element_text(size=10, hjust = 0.5)) +
    scale_x_continuous(expand = c(.01,.01)) +
    ylab("Median bin value") +
    ggtitle(paste("XY,", length(outlier_bins_XY),"outlier bins", sep=" ")) +
    labs(subtitle=basename(directory)) +
    ylim(0,5)
  ggsave(filename = paste(fin,"/bin_stuff/male_median_bin_values_XY.png",sep=''), plot=sdfg, height=4, width=6)



  #Integrate all bad bins
  outlier_bins_all <- c(outlier_bins_autosome, outlier_bins_XY)
  outlier_bins_male <- outlier_bins_all
  #Combine all medbinvals
  medbinvals <- rbind(autosome_medbinvals, XY_medbinvals)

  print("Plotting all male outlier bins")
  dfgh <- ggplot(medbinvals,
         mapping=aes(x=Bin, y=Median_bin_value)) +
    geom_vline(size=0.5,col="gray88",xintercept=autosome_chrom_breaks[-length(autosome_chrom_breaks)]) +
    geom_point(alpha=0.55, size=2) +
    #Autosomal outlier cutoffs
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=1, xend=last(autosome_chrom_breaks),
                 y=2^(autosome_log2_25pct - 1.5*IQR_autosome_log2), yend=2^(autosome_log2_25pct - 1.5*IQR_autosome_log2)) +
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=1, xend=last(autosome_chrom_breaks),
                 y=2^(autosome_log2_75pct + 1.5*IQR_autosome_log2), yend=2^(autosome_log2_75pct + 1.5*IQR_autosome_log2)) +
    #XY outlier cutoffs
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=last(autosome_chrom_breaks)+1, xend=last(XY_chrom_breaks),
                 y=2^(XY_log2_25pct - 1.5*IQR_XY_log2), yend=2^(XY_log2_25pct - 1.5*IQR_XY_log2)) +
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=last(autosome_chrom_breaks)+1, xend=last(XY_chrom_breaks),
                 y=2^(XY_log2_75pct + 1.5*IQR_XY_log2), yend=2^(XY_log2_75pct + 1.5*IQR_XY_log2)) +
    theme_bw(base_size=24) +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          plot.subtitle = element_text(size=12, hjust = 0.5)) +
    scale_x_continuous(expand = c(.01,.01)) +
    ylab("Median bin value") +
    ggtitle(paste("All chromosomes,", length(outlier_bins_male),"outlier bins", sep=" ")) +
    labs(subtitle=basename(directory)) +
    ylim(0,5)
  ggsave(filename = paste(fin,"/bin_stuff/male_median_bin_values_all_chromosomes.png",sep=''), plot=dfgh, height=4, width=14)


} 
if(sum(bictable$Sex=="female")>0) {
  X_binvals <-  binvals %>% filter(Bin %in% which(longfile$Chrom >= chrX),Cell %in% bictable$Sample[which(bictable$Sex=="female")])
  female_sample_names <- c()
  female_sample_names <- (bictable$Sample[which(bictable$Sex=="female")])
  X_medbinvals <- X_binvals %>%
    group_by(Bin) %>%
    summarize(Median_bin_value=median(as.numeric(Bin_value)))
  X_medbinvals <- X_medbinvals %>%
    mutate(log2_median_bin_value=log2(Median_bin_value))
  X_medbinvals <- cbind("Chrom"=longfile$Chrom[which(longfile$Chrom == chrX)],
                         "Start"=longfile$Start[which(longfile$Chrom == chrX)],
                         "End"=longfile$End[which(longfile$Chrom == chrX)],
                         X_medbinvals)
  X_chrom_breaks <- cumsum(rle(X_medbinvals$Chrom)$lengths) + last(autosome_chrom_breaks)


  X_log2_25pct <- as.numeric(quantile(X_medbinvals$log2_median_bin_value)[2])
  X_log2_75pct <- as.numeric(quantile(X_medbinvals$log2_median_bin_value)[4])
  IQR_X_log2 <- IQR(X_medbinvals$log2_median_bin_value)

  outlier_bins_X <- X_medbinvals$Bin[X_medbinvals$log2_median_bin_value < (X_log2_25pct - 1.5*IQR_X_log2)]
  outlier_bins_X <- sort(c(outlier_bins_X,
                            X_medbinvals$Bin[X_medbinvals$log2_median_bin_value > (X_log2_75pct + 1.5*IQR_X_log2)]))

  print("Plotting female chrX outlier bins")
  fghj <- ggplot(X_medbinvals,
         mapping=aes(x=Bin, y=Median_bin_value)) +
    geom_point(alpha=0.55, size=2) +
    geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(X_log2_25pct - 1.5*IQR_X_log2)) +
    geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(X_log2_75pct + 1.5*IQR_X_log2)) +
    theme_bw(base_size=14) +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.subtitle = element_text(size=8, hjust = 0.5)) +
    scale_x_continuous(expand = c(.01,.01)) +
    ylab("Median bin value") +
    ggtitle(paste("X chromosome,", length(outlier_bins_X),"outlier bins", sep=" ")) +
    labs(subtitle=basename(directory)) +
    ylim(0,5)
     ggsave(filename = paste(fin,"/bin_stuff/female_median_bin_values_X.png",sep=''), plot=fghj, height=4, width=4)


  Y_binvals <- binvals %>% filter(Bin %in% which(longfile$Chrom == chrX+1),Cell %in% bictable$Sample[which(bictable$Sex=="female")])
  Y_medbinvals <- Y_binvals %>%
    group_by(Bin) %>%
    summarize(Median_bin_value=median(as.numeric(Bin_value)))
  Y_medbinvals <- Y_medbinvals %>%
    mutate(log2_median_bin_value=log2(Median_bin_value))
  Y_medbinvals <- cbind("Chrom"=longfile$Chrom[which(longfile$Chrom == chrX+1)],
                        "Start"=longfile$Start[which(longfile$Chrom == chrX+1)],
                        "End"=longfile$End[which(longfile$Chrom == chrX+1)],
                        Y_medbinvals)
  Y_chrom_breaks <- cumsum(rle(Y_medbinvals$Chrom)$lengths) + X_chrom_breaks

  Y_log2_25pct <- as.numeric(quantile(Y_medbinvals$log2_median_bin_value)[2])
  Y_log2_75pct <- as.numeric(quantile(Y_medbinvals$log2_median_bin_value)[4])
  IQR_Y_log2 <- IQR(Y_medbinvals$log2_median_bin_value)

  outlier_bins_Y <- Y_medbinvals$Bin[Y_medbinvals$log2_median_bin_value < (Y_log2_25pct - 1.5*IQR_Y_log2)]
  outlier_bins_Y <- sort(c(outlier_bins_Y,
                           Y_medbinvals$Bin[Y_medbinvals$log2_median_bin_value > (Y_log2_75pct + 1.5*IQR_Y_log2)]))
  #In female populations, sometimes entire Y chrom is outliers -- check if that's happening
  if(length(outlier_bins_Y) == sum(longfile$Chrom == chrX+1)){
    print("Warning: All chrY bins are identified as outliers")
  }

  print("Plotting female chrY outlier bins")
  ghjk <- ggplot(Y_medbinvals,
         mapping=aes(x=Bin, y=Median_bin_value)) +
    geom_point(alpha=0.55, size=2) +
    geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(Y_log2_25pct - 1.5*IQR_Y_log2)) +
    geom_abline(size=1.5, col="orange2", linetype="dashed", slope=0, intercept=2^(Y_log2_75pct + 1.5*IQR_Y_log2)) +
    theme_bw(base_size=14) +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.subtitle = element_text(size=8, hjust = 0.5)) +
    scale_x_continuous(expand = c(.01,.01)) +
    ylab("Median bin value") +
    ggtitle(paste("Y chromosome,", length(outlier_bins_Y),"outlier bins", sep=" ")) +
    labs(subtitle=basename(directory)) +
    ylim(0,5)
    ggsave(filename = paste(fin,"/bin_stuff/female_median_bin_values_Y.png",sep=''), plot=ghjk, height=4, width=4)


  #Integrate all bad bins
  outlier_bins_all <- c(outlier_bins_autosome, outlier_bins_X, outlier_bins_Y)
  outlier_bins_female <-outlier_bins_all
  #Combine all medbinvals
  medbinvals <- rbind(autosome_medbinvals, X_medbinvals, Y_medbinvals)

  print("Plotting all female outlier bins")
  hjkl <- ggplot(medbinvals,
         mapping=aes(x=Bin, y=Median_bin_value)) +
    geom_vline(size=0.5,col="gray88",xintercept=autosome_chrom_breaks[-length(autosome_chrom_breaks)]) +
    geom_point(alpha=0.55, size=2) +
    #Autosomal outlier cutoffs
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=1, xend=last(autosome_chrom_breaks),
                 y=2^(autosome_log2_25pct - 1.5*IQR_autosome_log2), yend=2^(autosome_log2_25pct - 1.5*IQR_autosome_log2)) +
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=1, xend=last(autosome_chrom_breaks),
                 y=2^(autosome_log2_75pct + 1.5*IQR_autosome_log2), yend=2^(autosome_log2_75pct + 1.5*IQR_autosome_log2)) +
    #X outlier cutoffs
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=last(autosome_chrom_breaks)+1, xend=last(X_chrom_breaks),
                 y=2^(X_log2_25pct - 1.5*IQR_X_log2), yend=2^(X_log2_25pct - 1.5*IQR_X_log2)) +
    geom_segment(size=1.5, linetype="longdash", col="orange2",
                 x=last(autosome_chrom_breaks)+1, xend=last(X_chrom_breaks),
                 y=2^(X_log2_75pct + 1.5*IQR_X_log2), yend=2^(X_log2_75pct + 1.5*IQR_X_log2)) +
    #Y outlier cutoffs
    geom_segment(size=1.5, linetype="solid", col="orange2",
                 x=last(X_chrom_breaks)+1, xend=last(Y_chrom_breaks),
                 y=2^(Y_log2_25pct - 1.5*IQR_Y_log2), yend=2^(Y_log2_25pct - 1.5*IQR_Y_log2)) +
    geom_segment(size=1.5, linetype="solid", col="orange2",
                 x=last(X_chrom_breaks)+1, xend=last(Y_chrom_breaks),
                 y=2^(Y_log2_75pct + 1.5*IQR_Y_log2), yend=2^(Y_log2_75pct + 1.5*IQR_Y_log2)) +
    theme_bw(base_size=24) +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          plot.subtitle = element_text(size=12, hjust = 0.5)) +
    scale_x_continuous(expand = c(.01,.01)) +
    ylab("Median bin value") +
    ggtitle(paste("All chromosomes,", length(outlier_bins_female),"outlier bins", sep=" ")) +
    labs(subtitle=basename(directory)) +
    ylim(0,5)
  ggsave(filename = paste(fin,"/bin_stuff/female_median_bin_values_all_chromosomes.png",sep=''), plot=hjkl, height=4, width=14)

}

#Get percent of bins identified as outliers
output_male_outliers <- as.matrix(outlier_bins_male)
output_female_outliers <- as.matrix(outlier_bins_female)

pct_outliers_all <- 100*length(outlier_bins_all)/nrow(longfile)
pct_outliers_all <- sprintf("%.2f", round(pct_outliers_all, 2))
#output_outliers <- cbind(basename(directory), output_outliers)

male_info <- as.data.frame((male_sample_names))
colnames(male_info)<- c("Name")

female_info <- as.data.frame((female_sample_names))
colnames(female_info)<- c("Name")



write.table(output_male_outliers, file = paste(fin,"/bin_stuff/male_outlier_bins_all_chr.txt",sep=''),
            quote=F, row.names = F, col.names = F, sep="\t",append=T)

write.table(output_female_outliers, file = paste(fin,"/bin_stuff/female_outlier_bins_all_chr.txt",sep=''),
            quote=F, row.names = F, col.names = F, sep="\t",append=T)
write.table(male_info, file = paste(fin,"/bin_stuff/male_sample_names.txt",sep=''),quote=F,row.names=F,col.names=T, sep="\t",append=T)

write.table(female_info, file = paste(fin,"/bin_stuff/female_sample_names.txt",sep=''),quote=F,row.names=F,col.names=T, sep="\t",append=T)

#print(paste("Total outlier bins detected: ", length(outlier_bins_all), "/",nrow(longfile), " (",pct_outliers_all,"%)", sep=""))
output_outliers <- as.matrix(outlier_bins_all)
write.table(x = output_outliers, file = paste(fin,"/bin_stuff/outlier_bins_all_chr.txt",sep=''),
            quote=F, row.names = F, col.names = F, sep="\t",append=T)
print("Step 1 of the pipeline has been completed")