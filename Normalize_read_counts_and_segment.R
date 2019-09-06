#Contributed to by Alex Koeppel, Will Chronister, Inusah Diallo
#2019-09-06
#Some code adapted from Ginkgo (Aboukhalil and Garvin, CSHL); see bottom section
#Last major update 9/14/18: updated output files for new pipeline
#Full cov.data.final now output in addition to .long/.short/.num
#Improved CNV profile plotting code in various ways
#BIC calculated and output in BIC_info table (appended)

#Set arguments
args <- commandArgs(trailing=T)
#options(stringsAsFactors = F)
#11/3/16 modification here, only takes one argument
#directory <- "~/r_script_tests/"
directory <-args[1] 
path <- args[2]
#infile <- "3003__Lieber-Pool78_trimmed_Pool78_G3_S86_L002_R1_trimmed_Pool78_G3_S86_L002_R2_500kb_SORT.cov"
infile <- args[3]
#infile <-"SRR1006431_500kb.cov"




#directory<-"/.automount/nas8-s/export/vol97/mcconnell_lab/cnvpipe/Fibroblast_data_50kb/"
#infile <-"SRR1006431_50kb_SORT.cov"

splitfname <- strsplit(infile,"\\.")[[1]]
if(length(splitfname)==2){
  header <- splitfname[1]
} else if(length(splitfname)>2){
  header <- splitfname[1]
  for(i in 2:(length(splitfname)-1)){
    namepiece <-splitfname[i]
    header <-paste(header, namepiece, sep=".")
  }
}

sampname <- strsplit(header,"_SORT")[[1]][1]
#splitsampname <- strsplit(header,"_")[[1]]
#if(length(splitsampname)==2){
#  sampname <- splitsampname[1]
#} else if(length(splitsampname)>2){
#  sampname <- splitsampname[1]
#  for(i in 2:(length(splitsampname)-2)){
#    sampnamepiece <-splitsampname[i]
#    sampname <-paste(sampname, sampnamepiece, sep="_")
#  }
#}

#Set output filenames
outfile <-paste(directory,"/",header,"_norm.num", sep='')
outLong <- paste(directory,"/",header,"_norm.num.long",sep='')
outShort <- paste(directory,"/",header,"_norm.num.short",sep='')
outPlot <- paste(directory,"/",header,"_PLOT.png",sep='')
outAllBinInfo <- paste(directory,"/",header,"_norm.num.AllBinInfo.long",sep='')
outBICs <- paste(directory,"/","BICs_and_other_info_badbinsIN.txt",sep='')
outBinValues <- paste(directory,"/",header,"_bin_CN_values.txt",sep='')
outBinRawCounts <- paste(directory,"/",header,"_bin_raw_counts.txt",sep='')
outSoSPlot <- paste(directory,"/",header,"_SoS.png",sep='')
outSoSinfo <- paste(directory,"/","SoS_CN_multiplier_info.txt",sep='')
outSegmentData <- paste(directory,"/",header,"_segment_data.txt",sep='')

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
pkgTest("MASS")


#Load packages
library(DNAcopy)
library(MASS) #MASS pkg usually built-in
#library(sqldf) #Package to use SQL-style databases (for loading big tables).
library(ggplot2)

#Set working directory and read in file.
#setwd(directory)
#setwd("/Volumes/uvabx/projects/mcconnell/CNV_seqPipeLine/Fibroblast_data")

fullfilepath <-paste(directory,infile,sep='/')
#f <- file(fullfilepath)
cov.data <- read.table(fullfilepath, header=FALSE)
#cov.data <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = F, row.names = F, sep='\t'))

#Sometimes cov.data has 1:4505 in its first column; remove column if so
if(ncol(cov.data)==7){
  cov.data <- cov.data[,2:7]
}

names(cov.data)<-c("Chrom","Start","End","GCfrac","Nfrac","Count")

chrXnum <- max(cov.data$Chrom) - 1
chrYnum <- max(cov.data$Chrom)

#Eliminate sex chromosomes
#cov.data<-subset(cov.data, Chrom!='X')
#cov.data<-subset(cov.data, Chrom!='Y')

#Note:  These are Mike's bins designed to produce a fairly even distribution of counts
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
   }
   else{
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
       print(paste("Resolving issue by retaining lower outliers from GC bin ",k,sep=""))
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
       print(paste("Resolving issue by retaining upper outliers from GC bin ",k,sep=""))
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
write.table(cov.data.final, file=outfile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)  
#write.table(cov.data.final, file="SRR1006431_NORM.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

#Code for Segmenting with CBS (modified version of Mike's code)
#Separate out the autosome estimates
autosomeCNval <- cov.data.final$est_copy[which(!(cov.data.final$Chrom %in% c("chrX","X","chrY","Y",chrXnum,chrYnum)))]
#Separate out the X chromosome estimates
XCNval<- cov.data.final$est_copy[which((cov.data.final$Chrom %in% c("chrX","X",chrXnum)))]
#Compute median and mad for autosomal values
C1_med<-median(autosomeCNval)
C1_mad<-mad(autosomeCNval)
tC1_mad<-2*C1_mad
#Compute median and mad for X chrom. values
C1_xmed<-median(XCNval)
C1_xmad<-mad(XCNval)
#Read in data as CNA (copy number array) objectm and smooth.
C1_cna=CNA(cov.data.final$est_copy, cov.data.final$Chrom,cov.data.final$Start,data.type="logratio",sampleid=sampname)
C1_cna_smoothed<-smooth.CNA(C1_cna)
#Segment with CBS algorithm, and summarize
C1_seg = segment(C1_cna_smoothed,alpha=0.001, min.width = 5, undo.splits="sdundo", undo.SD=0, verbose = 2)
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

write.table(cov.data.final, file = outAllBinInfo, col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE) 

#Cut out extraneous columns to match Mike's output.  Comment this if we ever want the extra columns.
longtable<-cov.data.final[,c(1:6,8:9,11:16)]
shorttable<-cbind(C1_seg[[3]],C1_seg_summary[,2:9],C1_med,C1_mad,C1_xmed,C1_xmad)
write.table(longtable, file = outLong,col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE); 
write.table(shorttable, file = outShort,col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE); 

###NEW ADDITIONS TO SCRIPT
#Calculate BIC (formula I -- uniform variance across segments)
total_var <- sum((cov.data.final$est_copy - cov.data.final$seg.mean)^2)/nrow(cov.data.final)
ln_total_var <- log(total_var)

num_segments <- nrow(shorttable)
num_changepoints <- num_segments - 1
kp <- 1 + 2*(num_changepoints)

seg_penalty <- (kp/(nrow(cov.data.final)))*(log(nrow(cov.data.final)))

BIC <- ln_total_var + seg_penalty
BIC <- sprintf("%.3f", round(BIC, 3))
ln_total_var <- sprintf("%.4f", round(ln_total_var, 4))
seg_penalty <- sprintf("%.4f", round(seg_penalty, 4))

#Calculate phi for null model
# longtable_arstuff <- ar(longtable$est_copy[1:4272], order.max = 1, aic=F)
# Filter out autosomal rows where est_copy == "NA"
library(dplyr)
autosomal_cov.data.final <- cov.data.final %>% filter(as.numeric(Chrom) < max(as.numeric(cov.data.final$Chrom))-1)
cov.data.final_arstuff <- ar(autosomal_cov.data.final$est_copy, order.max = 1, aic=F)
phi <- cov.data.final_arstuff$ar
var <- cov.data.final_arstuff$var.pred

phi <- sprintf("%.3f", round(phi, 3))
var <- sprintf("%.3f", round(var, 3))

#Calculate phi for alt. model (subtract seg. means from bin values)
cov.data.final_arstuff_meansubd <- ar((autosomal_cov.data.final$est_copy-autosomal_cov.data.final$seg.mean), order.max = 1, aic=F)
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


dirs <- path
rundir <- basename(dirs)
run <- strsplit(rundir,split="_")[[1]][1]
wga <- strsplit(rundir,split="_")[[1]][2]
species <- strsplit(rundir,split="_")[[1]][3]
indivname <- strsplit(rundir,split="_")[[1]][4]
# Calculate median Xchr value to identify sex 
if(species == "human"){
  chrX <- 23
  chrX_binvals <- longtable[which(longtable$Chrom == chrX),10]
  print(chrX_binvals)
  chrX_median_binval <- median(as.vector(as.matrix(chrX_binvals)))
} else if (species == "mouse") {
  chrX <- 20
  chrX_binvals <- longtable[which(longtable$Chrom == chrX),10]
  chrX_median_binval <- median(as.vector(as.matrix(chrX_binvals)))
}
print(paste(sampname,chrX_median_binval,sep=":"))
#Determine if male or female based on median chrX bin value;
#This will influence outlier bins identified/removed
if(chrX_median_binval < 1.5){
  sex <- "male"
} else {
  sex <- "female"
}
BIC_info <- as.data.frame(t(as.matrix(c(rundir, run, wga, species, indivname,sex,bad_bin_status, sampname, BIC, ln_total_var, seg_penalty, num_segments, phi, var, phi_meansubd, var_meansubd, totalsampreads, meanperbin))))
colnames(BIC_info) <- c("Run_dir","Run","WGA","Species","Individual","Sex","Bad_bin_status","Sample","BIC","ln_total_var","Segment_penalty","Segments","Phi","Variance","Mean_subtracted_phi","Variance_mean_subtracted","Total_reads","Reads_per_bin")
write.table(BIC_info, file = outBICs, quote = F, sep="\t",row.names = F, col.names = !file.exists(outBICs),append=T) 

bin_values <- as.data.frame(t(as.matrix(c(cov.data.final$est_copy,sampname))))
colnames(bin_values) <- c(paste("Bin",sprintf('%0.4d', 1:nrow(cov.data.final)),sep="_"),"Sample")
write.table(bin_values, file = outBinValues, quote = F, sep="\t",row.names = F, col.names = T) 

raw_counts <- as.data.frame(t(as.matrix(c(cov.data.final$Count,sampname))))
colnames(raw_counts) <- c(paste("Bin",sprintf('%0.4d', 1:nrow(cov.data.final)),sep="_"),"Sample")
write.table(raw_counts, file = outBinRawCounts, quote = F, sep="\t",row.names = F, col.names = T) 

shorttable2 <- cbind(rundir, run, wga, species, indivname, bad_bin_status, sampname, BIC, shorttable)
colnames(shorttable2) <- c("Run_dir","Run","WGA","Species","Individual","Bad_bin_status","Sample","BIC",colnames(shorttable))
write.table(shorttable2, file = outSegmentData, quote = F, sep="\t",row.names = F, col.names = T) 
# all_shorttable_data <- rbind(all_shorttable_data, shorttable2)

#Generate column with absolute positions
newpos<-vector()
lastval <-0
for(y in 1:chrYnum){
  tempchromtable <- subset(longtable, Chrom==y)
  fixpos <- tempchromtable$Start + lastval
  lastval <- lastval + tempchromtable$Start[length(tempchromtable$Start)]
  newpos<-c(newpos,fixpos)
}

longtable$Position <- newpos

#Set all copy values above max to max (5)
tflist <-longtable$est_copy > 5
#numover <-nrow(longtable[tflist,])
longtable$est_copy[tflist] <- 5.0

#Generate positions for chromosome labels
chr1pos <- median(longtable[longtable$Chrom =="1",]$Position, na.rm = T)
chr5pos <- median(longtable[longtable$Chrom =="5",]$Position, na.rm = T)
chr10pos <- median(longtable[longtable$Chrom =="10",]$Position, na.rm = T)
chr15pos <- median(longtable[longtable$Chrom =="15",]$Position, na.rm = T)
chrXpos <- median(longtable[longtable$Chrom ==as.character(chrXnum),]$Position, na.rm = T)
chrYpos <- median(longtable[longtable$Chrom ==as.character(chrYnum),]$Position, na.rm = T)
label.df<-data.frame(label=c("chr1","chr5","chr10","chr15","chrX","Y"),
                     positions=c(chr1pos,chr5pos,chr10pos,chr15pos,chrXpos,chrYpos))

#NEW PLOTTING CODE
h <- ggplot(longtable, aes(x=Position, y=est_copy, color=factor(Chrom))) +
  geom_point(size=1) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none",plot.title = element_text(size=12)) +
  scale_color_manual(values=rep(c("green","lightblue"),12)) +
  geom_line(y=longtable$seg.median, col="red",size=1) +
  geom_line(y=C1_med + tC1_mad, col="black", size=1, linetype="dashed") +
  geom_line(y=C1_med - tC1_mad, col="black", size=1, linetype="dashed") +
  geom_line(y=C1_med + C1_mad, col="black", size=1, linetype="dashed") +
  geom_line(y=C1_med - C1_mad, col="black", size=1, linetype="dashed") +
  ylim(0,5) +
  ylab("Copy Number") +
  ggtitle(paste(sampname, ", phi ", phi, ", BIC ", BIC, sep="")) +
  scale_x_continuous(breaks=c(label.df$positions), labels=label.df$label, expand = c(0,0))
ggsave(outPlot, plot = h, height = 4, width = 12) 




# ###WORK IN PROGRESS
# #Borrowed Ginkgo code for CN multiplier to find minimum Sum of Squares Error
# Ginkgo -- Aboukhalil and Garvin, CSHL
# Map breakpoints to kth sample
loc <- cov.data.final[,1:2]
frag <- C1_seg$output[,2:3]
len = dim(frag)[1]
bps = array(0, len)
for (j in 1:len){
  bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
}
bps = sort(bps)
#bps[(len=len+1)] = l
bps[(len=len+1)] = nrow(cov.data.final)

# Track global breakpoint locations
breaks = matrix(0,nrow(cov.data.final),1)
breaks[bps] = 1

# Modify bins to contain median read count/bin within each segment
fixed = matrix(0,nrow(cov.data.final),1)
#fixed[,1][1:bps[2]] = median(normal[,k][1:bps[2]])
fixed[,1][1:bps[2]] = median(longtable$est_copy[1:bps[2]])

for(i in 2:(len-1)){
  fixed[,1][bps[i]:(bps[i+1]-1)] = median(longtable$est_copy[bps[i]:(bps[i+1]-1)])
}
fixed[,1] = fixed[,1]/mean(fixed[,1])

# ----------------------------------------------------------------------------
# -- Determine Copy Number   (SoS Method)
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


abbrev_CN_mult_data <- CN_mult_data %>% arrange(SoS) %>% head(50) %>% select(CN_multiplier, SoS)
abbrev_CN_mult_data <- cbind.data.frame("Cell"=sampname, "Rank"=seq(1,50), abbrev_CN_mult_data)
colnames(abbrev_CN_mult_data)[4] <- "SoS_error"

write.table(abbrev_CN_mult_data, file = outSoSinfo, quote = F, sep="\t",row.names = F, col.names = !file.exists(outSoSinfo), append=T) 



# w=1 # 1 sample
# maxPloidy=6
# CNmult = matrix(0,5,1)
# outerColsums = matrix(0, (20*(maxPloidy-1.5)+1), w)
#
# CNgrid = seq(1.5, maxPloidy, by=0.05)
# outerRaw = unique(longtable$seg.median)/2 %o% CNgrid
# #outerRaw = fixed[,k] %o% CNgrid
# outerRound = round(outerRaw)
# outerDiff = (outerRaw - outerRound) ^ 2
# outerColsums[,1] = colSums(outerDiff, na.rm = FALSE, dims = 1)
# CNmult[,1] = CNgrid[order(outerColsums[,1])[1:5]]
#
#
# outerRaw = sort(unique(longtable$seg.median)/2) %o% CNgrid
# #outerRaw = fixed[,k] %o% CNgrid
# outerRound = round(outerRaw)
# outerDiff = (outerRaw - outerRound) ^ 2
# outerColsums[,1] = colSums(outerDiff, na.rm = FALSE, dims = 1)
# CNmult[,1] = CNgrid[order(outerColsums[,1])[1:5]]
