#Contributors: Alex Koeppel, Will Chronister, Inusah Diallo
#2019-09-06
#Some code adapted from Ginkgo (Aboukhalil and Garvin, CSHL); see bottom section
#This script removes outlier bins identified at end of step1 and reanalyzes for CNVs using DNAcopy.

#setwd("~/final_dataset_exclude_badbins/")

#turn off scientific notation
options(scipen = 999)
#turn off stringsAsFactors
options(stringsAsFactors = FALSE)
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
args <- commandArgs(trailing=T)
directory <- args[1]
fin <- args[2]
rundir <- basename(directory)
run <- strsplit(rundir,split="_")[[1]][1]
wga <- strsplit(rundir,split="_")[[1]][2]
species <- strsplit(rundir,split="_")[[1]][3]
indivname <- strsplit(rundir,split="_")[[1]][4]
#Load packages
library(DNAcopy)
library(MASS) #MASS pkg usually built-in
#library(sqldf) #Package to use SQL-style databases (for loading big tables).
library(ggplot2)
library(dplyr)

#Open data frames
BIC_table <- data.frame()
all_shorttable_data <- data.frame()

#Set DNAcopy parameters
alphaval = 0.001
undoSDval = 0
minwidthval = 5

#Read in bad bin table
bad_bin_table <- read.table(paste(fin,"/bin_stuff/outlier_bins_all_chr.txt",sep=''), header = F) 
bad_bins <- c()
female_bad_bin_table <- read.table(paste(fin,"/bin_stuff/female_outlier_bins_all_chr.txt",sep=''), header = F) 
male_bad_bin_table <- read.table(paste(fin,"/bin_stuff/male_outlier_bins_all_chr.txt",sep=''), header = F) 

female_bad_bins <- as.vector(female_bad_bin_table[,1])
male_bad_bins <- as.vector(male_bad_bin_table[,1])

male_names<- read.table(paste(fin,"/bin_stuff/male_sample_names.txt",sep=''),header=T)
female_names<- read.table(paste(fin,"/bin_stuff/female_sample_names.txt",sep=''),header=T)



infiles <- list.files(path=paste(fin,"/cov_files_badbinsIN/",sep=''), pattern = "cov", full.names = T) 

for(i in 1:length(infiles)){

  infile <- infiles[i]
  splitfname <- strsplit(basename(infile),"\\.")[[1]]
  if(length(splitfname)==2){
    header <- splitfname[1]
  } else if(length(splitfname)>2){
    header <- splitfname[1]
    for(i in 2:(length(splitfname)-1)){
      namepiece <-splitfname[i]
      header <-paste(header, namepiece, sep=".")
    }
  }

  rev <- strsplit(header,"_")
  if(rev[[1]][1]=="trimmed"){
    patts <- paste(rev[[1]][1],rev[[1]][2],sep='_')
  }
  else{
    patts <- paste(rev[[1]][1],"_",sep='')
  }


  if( length( list.files( path = directory, pattern = patts ) ) > 0){

    #print(header)
  #if male, use male outlier bins; if female, use female outlier bins
  if( sum( grep (patts, male_names) ) >0){
    bad_bins <- male_bad_bins
  }else if( sum( grep (patts, female_names) ) >0){
    bad_bins <-female_bad_bins
  }

  sampname <- strsplit(header,"_SORT")[[1]][1]
  print(paste("Now removing",length(bad_bins),"bad bins",sep=" "))
    #Check if file produced a plot first time; if not, skip
    if(file.exists(paste(fin,"/pngs_badbinsIN/",header,"_PLOT.png",sep=''))){
      #  outfile <-paste(directory,header,"_norm.num", sep='')
      #  outLong <- paste(directory,header,"_norm.num.long",sep='')
      #  outShort <- paste(directory,header,"_norm.num.short",sep='')
      outPlot <- paste(fin,"/pngs_badbinsOUT/",header,"_PLOT_badbinsOUT.png",sep='')
      outAllBinInfo <- paste(fin,"/full_longs_badbinsOUT/",header,"_norm.num.AllBinInfo_badbinsOUT.long",sep='')
      outBICs <- paste(fin,"/bic_outputs/BICs_and_other_info_badbinsOUT.txt",sep='')
      #  outBinValues <- paste(directory,header,"_bin_CN_values.txt",sep='')
      #  outBinRawCounts <- paste(directory,header,"_bin_raw_counts.txt",sep='')
      outSoSPlot <- paste(fin,"/sos_CN_mult_badbinsOUT/",header,"_SoS_badbinsOUT.png",sep='')
      outSoSinfo <- paste(fin,"/sos_CN_mult_badbinsOUT/SoS_CN_multiplier_info_badbinsOUT.txt",sep='')


      #source("http://bioconductor.org/biocLite.R")
      #biocLite("DNAcopy")


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

      pkgTestBC("DNAcopy")
      #pkgTest("sqldf")
      pkgTest("ggplot2")

      #Load packages
      library(DNAcopy)
      library(MASS) #MASS pkg usually built-in
      #library(sqldf) #Package to use SQL-style databases (for loading big tables).
      library(ggplot2)

      #Set working directory and read in file.
      #setwd(directory)
      #setwd("/Volumes/uvabx/projects/mcconnell/CNV_seqPipeLine/Fibroblast_data")

      cov.data <- read.table(infile, header=FALSE)

      #Sometimes cov.data has 1:lastbin in its first column; remove column if so
      if(ncol(cov.data)==7){
        cov.data <- cov.data[,2:7]
      }

      names(cov.data)<-c("Chrom","Start","End","GCfrac","Nfrac","Count")

      chrXnum <- max(cov.data$Chrom) - 1
      chrYnum <- max(cov.data$Chrom)

      #Note:  These are Mike's [GC] bins designed to produce a fairly even distribution of counts
      #across the bins.  These may need to be modified for other datasets.
      binSize <- c(0,34:45,47,49,53,100)
      binSize <- binSize/100

      #Assign each row of the coverage file to a bin and add the bin id to the data frame.
      binid<-vector()
      for(i in 1:nrow(cov.data)){
        thebin <-0
        gc <- cov.data$GCfrac[i]
        for(j in 1:(length(binSize)-1)){
          low <-binSize[j]
          high <- binSize[j+1]
          if(low <= gc){
            if(gc < high){
              thebin <- j
              binid<-c(binid,thebin)
            }
          }
        }
      }

      #Adds the bin id column to the data frame
      cov.data$bin <- binid

      #Compute total number of reads
      totreads <- sum(cov.data$Count)
      millreads <- totreads/1000000

      #Loop through list of the bin ids.
      #Remove outliers from each bin.
      cov.data.final<-{}
      for(k in 1:16){
        temptable <- subset(cov.data, bin==k)
        countmad <-mad(temptable$Count)
        countmean <-mean(temptable$Count)
        countmax <- countmean + (countmad * 4)
        countmin <- countmean - (countmad * 4)
        indexlist <-vector()
        for(m in 1:length(temptable$Count)){
          thecount <- temptable$Count[m]
          if(thecount < countmin){
            indexlist <- c(indexlist,m)
          }
          else if(thecount > countmax){
            indexlist <- c(indexlist,m)
          }
        }
        if(length(indexlist)>0){
          temptable.noout <- temptable[-indexlist,]
        } else {
          temptable.noout <-temptable
        }
        temptable.noout<-subset(temptable.noout, Chrom!=chrXnum)
        temptable.noout<-subset(temptable.noout, Chrom!=chrYnum)
        if(nrow(temptable.noout) == 0){
          print(paste("Normalization issue for cell: ",sampname,sep=""))
          print(paste("All bin values are outliers in GC bin ",k,". See below:",sep=""))
          print(temptable)
          #Solution: determine whether low counts or high counts are more realistic
          overall_median_readcount <- median(cov.data$Count)
          #If overall median reads/bin is closer to lower outlier threshold, then only remove upper outliers
          if(abs(countmin - overall_median_readcount) <= abs(countmax - overall_median_readcount)){
            print(paste("Resolving issue by excluding upper outliers from GC bin ",k,sep=""))
            indexlist <-vector()
            for(m in 1:length(temptable$Count)){
              thecount <- temptable$Count[m]
              if(thecount > countmax){
                indexlist <- c(indexlist,m)
              }
            }
            temptable.noout <- temptable[-indexlist,]
            temptable.noout<-subset(temptable.noout, Chrom!=chrXnum)
            temptable.noout<-subset(temptable.noout, Chrom!=chrYnum)
            #Otherwise, remove only lower outliers
          } else {
            print(paste("Resolving issue by excluding lower outliers from GC bin ",k,sep=""))
            indexlist <-vector()
            for(m in 1:length(temptable$Count)){
              thecount <- temptable$Count[m]
              if(thecount < countmin){
                indexlist <- c(indexlist,m)
              }
            }
            temptable.noout <- temptable[-indexlist,]
            temptable.noout<-subset(temptable.noout, Chrom!=chrXnum)
            temptable.noout<-subset(temptable.noout, Chrom!=chrYnum)
          }
        }
        #   binmedian <-median(temptable$Count)
        binmedian <-median(temptable.noout$Count)
        thestats <-fitdistr(temptable.noout$Count, "normal")
        binmean <- thestats$estimate[1]
        binsd <- thestats$estimate[2]
        zvals <-vector()
        pvals <-vector()
        l2fvals <-vector()
        rpkmvals <-vector()
        estvals <-vector()
        for(n in 1:length(temptable$Count)){
          thecount <- temptable$Count[n]
          #     countz <- (binmean - thecount) / (binsd/sqrt(length(temptable.noout$Count)))
          countz <- (thecount - binmean) / (binsd)
          #     countp <- dnorm(countz)
          countp <- 2*pnorm(-abs(countz))
          countl2f <- log2(thecount/binmean)
          #     countrpkm <- thecount / 500 / millreads
          countrpkm <- thecount / millreads
          countest <- 2*(thecount / binmedian)
          zvals <-c(zvals, countz)
          pvals <-c(pvals, countp)
          l2fvals <-c(l2fvals, countl2f)
          rpkmvals <-c(rpkmvals, countrpkm)
          estvals <-c(estvals, countest)
        }
        temptable$zscore <- zvals
        temptable$pvalue <- pvals
        temptable$Log2_foldchange <- l2fvals
        temptable$rpkm <- rpkmvals
        temptable$est_copy <- estvals
        cov.data.final <-rbind(cov.data.final, temptable)
      }

      #Sort the data frame by the chromosome, then the start position
      cov.data.final <-cov.data.final[with(cov.data.final, order(Chrom, Start)), ]
      #Write the table to a file.
      #write.table(cov.data.final, file=outfile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
      #write.table(cov.data.final, file="SRR1006431_NORM.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

      #Change bad bins to NA
      cov.data.final$est_copy[bad_bins] <- NA

      #Code for Segmenting with CBS (modified version of Mike's code)
      #Separate out the autosome estimates
      autosomeCNval <- cov.data.final$est_copy[which(!(cov.data.final$Chrom %in% c("chrX","X","chrY","Y",chrXnum,chrYnum)))]
      #Separate out the X chromosome estimates
      XCNval<- cov.data.final$est_copy[which((cov.data.final$Chrom %in% c("chrX","X",chrXnum)))]
      #Compute median and mad for autosomal values
      C1_med<-median(autosomeCNval, na.rm = T)
      C1_mad<-mad(autosomeCNval, na.rm = T)
      tC1_mad<-2*C1_mad
      #Compute median and mad for X chrom. values
      C1_xmed<-median(XCNval, na.rm = T)
      C1_xmad<-mad(XCNval, na.rm = T)
      #Read in data as CNA (copy number array) objectm and smooth.
      C1_cna=CNA(cov.data.final$est_copy, cov.data.final$Chrom,cov.data.final$Start,data.type="logratio",sampleid=sampname)
      C1_cna_smoothed<-smooth.CNA(C1_cna)
      #Segment with CBS algorithm, and summarize
      C1_seg = segment(C1_cna_smoothed,alpha=alphaval, min.width = 5, undo.splits="sdundo", undo.SD=0, verbose = 1)
      C1_seg_summary = segments.summary(C1_seg)
      outtablenames <- colnames(cov.data.final)
      for(i in 1:dim(C1_seg[[3]])[1]) {
        cov.data.final[C1_seg[[3]][i,1]:C1_seg[[3]][i,2],13] = C1_seg_summary[i,6]
        cov.data.final[C1_seg[[3]][i,1]:C1_seg[[3]][i,2],14] = C1_seg_summary[i,7]
        cov.data.final[C1_seg[[3]][i,1]:C1_seg[[3]][i,2],15] = C1_seg_summary[i,8]
        cov.data.final[C1_seg[[3]][i,1]:C1_seg[[3]][i,2],16] = C1_seg_summary[i,9]
      }
      outtablenames<-c(outtablenames,"seg.mean", "seg.sd", "seg.median", "seg.mad")
      colnames(cov.data.final)<-outtablenames
      #WHY DOES THIS NOT APPEND?
      write.table(cov.data.final, file = outAllBinInfo, col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE) 

      #Cut out extraneous columns to match Mike's output.  Comment this if we ever want the extra columns.
      longtable<-cov.data.final[,c(1:6,8:9,11:16)]
      shorttable<-cbind(C1_seg[[3]],C1_seg_summary[,2:9],C1_med,C1_mad,C1_xmed,C1_xmad)
      #write.table(longtable, file = outLong,col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE);
      #write.table(shorttable, file = outShort,col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE);

      #Calculate BIC (formula I -- uniform variance across segments)
      longtable_noNA <- filter(longtable, est_copy != "NA")
      total_var <- sum((longtable_noNA$est_copy - longtable_noNA$seg.mean)^2)/nrow(longtable_noNA)
      ln_total_var <- log(total_var)

      num_segments <- nrow(shorttable)
      num_changepoints <- num_segments - 1
      kp <- 1 + 2*(num_changepoints)

      seg_penalty <- (kp/(nrow(longtable_noNA)))*(log(nrow(longtable_noNA)))

      BIC <- ln_total_var + seg_penalty
      BIC <- sprintf("%.3f", round(BIC, 3))

      ln_total_var <- sprintf("%.4f", round(ln_total_var, 4))
      seg_penalty <- sprintf("%.4f", round(seg_penalty, 4))

      # Calculate phi for null model
      # Filter out autosomal rows where est_copy == "NA"
      library(dplyr)
      autosomal_cov.data.final_noNA <- cov.data.final %>%
        filter(as.numeric(Chrom) < max(as.numeric(cov.data.final$Chrom))-1) %>%
        filter(est_copy != "NA")
      cov.data.final_arstuff <- ar(autosomal_cov.data.final_noNA$est_copy, order.max = 1, aic=F)
      phi <- cov.data.final_arstuff$ar
      var <- cov.data.final_arstuff$var.pred

      phi <- sprintf("%.3f", round(phi, 3))
      var <- sprintf("%.3f", round(var, 3))

      # Calculate phi for alt. model (subtract seg. means from bin values)
      # Filter out autosomal rows where est_copy == "NA"
      cov.data.final_arstuff_meansubd <- ar((autosomal_cov.data.final_noNA$est_copy-autosomal_cov.data.final_noNA$seg.mean), order.max = 1, aic=F)
      phi_meansubd <- cov.data.final_arstuff_meansubd$ar
      var_meansubd <- cov.data.final_arstuff_meansubd$var.pred

      phi_meansubd <- sprintf("%.3f", round(phi_meansubd, 3))
      var_meansubd <- sprintf("%.3f", round(var_meansubd, 3))

      #Create BIC table
      #Compute total reads and mean per window and append to file.
      meanperbin <- sprintf("%.1f", round(mean(cov.data.final$Count),1))
      totalsampreads <-sum(cov.data$Count)
      #meantotoutline <-paste(sampname,meanperbin,totalsampreads, sep="\t")
      #rawmeanfilename <- "All_counts_means.txt"
      #meanfilename <-paste(directory,rawmeanfilename,sep='')
      #write(meantotoutline, file=meanfilename, append=TRUE)
      if(exists("bad_bins")){
        bad_bin_status <- "bad_bins_excluded"
      } else {
        bad_bin_status <- "bad_bins_included"
      }

      BIC_info <- as.data.frame(t(as.matrix(c(rundir, run, wga, species, indivname, bad_bin_status, sampname, BIC, ln_total_var, seg_penalty, num_segments, phi, var, phi_meansubd, var_meansubd, totalsampreads, meanperbin))))
      colnames(BIC_info) <- c("Run_dir", "Run","WGA","Species","Individual","Bad_bin_status","Sample","BIC","ln_total_var","Segment_penalty","Segments","Phi","Variance","Mean_subtracted_phi","Variance_mean_subtracted","Total_reads","Reads_per_bin")
      write.table(BIC_info, file = outBICs, quote = F, sep="\t",row.names = F, col.names = !file.exists(outBICs),append=T) 

      shorttable2 <- cbind(rundir, run, wga, species, indivname, bad_bin_status, sampname, BIC, shorttable)
      colnames(shorttable2) <- c("Run_dir","Run","WGA","Species","Individual","Bad_bin_status","Sample","BIC",colnames(shorttable))
      all_shorttable_data <- rbind(all_shorttable_data, shorttable2)

      #  bin_values <- as.data.frame(t(as.matrix(c(cov.data.final$est_copy,sampname))))
      #  colnames(bin_values) <- c(paste("Bin",sprintf('%0.4d', 1:nrow(cov.data.final)),sep="_"),"Sample")
      #  write.table(bin_values, file = outBinValues, quote = F, sep="\t",row.names = F, col.names = T)

      #  raw_counts <- as.data.frame(t(as.matrix(c(cov.data.final$Count,sampname))))
      #  colnames(raw_counts) <- c(paste("Bin",sprintf('%0.4d', 1:nrow(cov.data.final)),sep="_"),"Sample")
      #  write.table(raw_counts, file = outBinRawCounts, quote = F, sep="\t",row.names = F, col.names = T)

      #Generate column with absolute positions
      newpos<-vector()
      lastval <-0
      for(y in 1:chrYnum){
        tempchromtable <- subset(longtable_noNA, Chrom==y)
        fixpos <- tempchromtable$Start + lastval
        lastval <- lastval + tempchromtable$Start[length(tempchromtable$Start)]
        newpos<-c(newpos,fixpos)
      }

      longtable_noNA$Position <- newpos

      #Set all copy values above max to max (5)
      tflist <-longtable_noNA$est_copy > 5
      #numover <-nrow(longtable[tflist,])
      longtable_noNA$est_copy[tflist] <- 5.0

      #Generate positions for chromosome labels
      chr1pos <- median(longtable_noNA[longtable_noNA$Chrom =="1",]$Position, na.rm = T)
      chr5pos <- median(longtable_noNA[longtable_noNA$Chrom =="5",]$Position, na.rm = T)
      chr10pos <- median(longtable_noNA[longtable_noNA$Chrom =="10",]$Position, na.rm = T)
      chr15pos <- median(longtable_noNA[longtable_noNA$Chrom =="15",]$Position, na.rm = T)
      chrXpos <- median(longtable_noNA[longtable_noNA$Chrom ==as.character(chrXnum),]$Position, na.rm = T)
      chrYpos <- median(longtable_noNA[longtable_noNA$Chrom ==as.character(chrYnum),]$Position, na.rm = T)
      label.df<-data.frame(label=c("chr1","chr5","chr10","chr15","chrX","Y"),
                          positions=c(chr1pos,chr5pos,chr10pos,chr15pos,chrXpos,chrYpos))

      #NEW PLOTTING CODE
      h <- ggplot(longtable_noNA, aes(x=Position, y=est_copy, color=factor(Chrom))) +
        geom_point(size=1) +
        theme_bw(base_size = 15) +
        theme(legend.position = "none",plot.title = element_text(size=12)) +
        scale_color_manual(values=rep(c("green","lightblue"),12)) +
        geom_line(y=longtable_noNA$seg.median, col="red",size=1) +
        geom_line(y=C1_med + tC1_mad, col="black", size=1, linetype="dashed") +
        geom_line(y=C1_med - tC1_mad, col="black", size=1, linetype="dashed") +
        geom_line(y=C1_med + C1_mad, col="black", size=1, linetype="dashed") +
        geom_line(y=C1_med - C1_mad, col="black", size=1, linetype="dashed") +
        ylim(0,5) +
        ylab("Copy Number") +
        ggtitle(paste(sampname, ", phi ", phi, ", BIC ", BIC, sep="")) +
        scale_x_continuous(breaks=c(label.df$positions), labels=label.df$label, expand = c(0,0))
      ggsave(outPlot, plot = h, height = 4, width = 12) 



      # #Borrowed Ginkgo code for CN multiplier to find minimum Sum of Squares Error
      # Map breakpoints to kth sample
      cov.data.final_noNA <- cov.data.final %>% filter(est_copy != "NA")
      loc <- cov.data.final_noNA[,1:2]
      frag <- C1_seg$output[,2:3]
      len = dim(frag)[1]
      bps = array(0, len)
      for (j in 1:len){
        bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
      }
      bps = sort(bps)
      #bps[(len=len+1)] = l
      bps[(len=len+1)] = nrow(cov.data.final_noNA)

      # Track global breakpoint locations
      breaks = matrix(0,nrow(cov.data.final_noNA),1)
      breaks[bps] = 1

      # Modify bins to contain median read count/bin within each segment
      fixed = matrix(0,nrow(cov.data.final_noNA),1)
      #fixed[,1][1:bps[2]] = median(normal[,k][1:bps[2]])
      fixed[,1][1:bps[2]] = median(longtable_noNA$est_copy[1:bps[2]])

      for(i in 2:(len-1)){
        fixed[,1][bps[i]:(bps[i+1]-1)] = median(longtable_noNA$est_copy[bps[i]:(bps[i+1]-1)])
      }
      fixed[,1] = fixed[,1]/mean(fixed[,1])

      # ----------------------------------------------------------------------------
      # -- Determine Copy Number (SoS Method)
      # ----------------------------------------------------------------------------

      # Determine Copy Number
      maxPloidy = 8
      CNgrid           = seq(0.5, maxPloidy, by=0.05)
      outerRaw         = fixed[,1] %o% CNgrid
      outerRound       = round(outerRaw)
      outerDiff        = (outerRaw - outerRound) ^ 2
      outerColsums = colSums(outerDiff, na.rm = FALSE, dims = 1)
      #CN multiplier
      CNmult       = CNgrid[order(outerColsums)]
      CNerror      = round(sort(outerColsums), digits=2)

      minSoS <- outerColsums==min(outerColsums)
      CN_mult_data <- data.frame("CN_multiplier"=CNgrid,"SoS"=outerColsums,
                                minimum_SoS=minSoS)

      hh <- ggplot(CN_mult_data,aes(CN_multiplier, SoS)) +
        scale_color_manual(values=c("black","red")) +
        geom_point() + geom_line() +
        geom_point(x=CN_mult_data$CN_multiplier[CN_mult_data$minimum_SoS==T],
                  y=CN_mult_data$SoS[CN_mult_data$minimum_SoS==T],
                  color="red",size=3, pch=17) +
        xlab("Copy Number Multiplier") +
        ylab("Sum of Squares Error") +
        theme_bw(base_size = 16) +
        theme(legend.position = "none",plot.title = element_text(size=12)) +
        ggtitle(sampname)
      ggsave(outSoSPlot, hh, height=6, width=10)


      abbrev_CN_mult_data <- CN_mult_data %>% arrange(SoS) %>% head(50) %>% dplyr::select(CN_multiplier, SoS)
      abbrev_CN_mult_data <- cbind.data.frame("Cell"=sampname, "Rank"=seq(1,50), abbrev_CN_mult_data)
      colnames(abbrev_CN_mult_data)[4] <- "SoS_error"

      write.table(abbrev_CN_mult_data, file = outSoSinfo, quote = F, sep="\t",row.names = F, col.names = !file.exists(outSoSinfo), append=T) 

    }
  }


}
write.table(all_shorttable_data,paste(fin,"/segment_data/all_segment_data_badbinsOUT.txt",sep=''),quote=F,row.names=F,col.names=!file.exists(paste(fin,"/segment_data/all_segment_data_badbinsOUT.txt",sep='')),sep="\t",append=T)

print("Step 2 of the pipeline has been completed")