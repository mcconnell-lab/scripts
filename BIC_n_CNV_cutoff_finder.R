# Contributors: Will Chronister, Inusah Diallo
# 2019-09-06

#This code will identify a cutoff for good quality cells and robust CNVs
#and output lists of cells and CNVs according to whether they passed cutoffs
#Will also plot data to visualize BIC scores and CNVs
#May require manual tweaking of mixed models to find best fit and cutoffs for data

#In the past, this script would incorporate previous results to fill out BIC/CNV distributions but
#this is no longer desirable (similar to decision made re: outlier bin analysis in Bad_bin_finder.R)



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
                install.packages(x,repos = "http://cran.us.r-project.org")
                if(!require(x,character.only=TRUE)) stop("Package not found")
        }
}

pkgTest("dplyr")
pkgTest("ggplot2")
pkgTest("tidyr")
pkgTest("mixtools")

library(dplyr)
library(mixtools)
library(ggplot2)

options(stringsAsFactors = F)

#When Rscript is run, add number of gaussian dists. to use for BIC cutoff
args <- commandArgs(trailing=T)
#directory <- args[1]
fin <- args[1]
num_gaussians <- as.numeric(args[2])
dir_names <- args[3]

#setwd("~/r_script_tests/BIC_tables/")
#directory <- "/nv/vol97/mcconnell_lab/cnvpipe/vandenBos2016_strandseq_human_AD07"
#rundir <- basename(directory)
#run <- strsplit(rundir,split="_")[[1]][1]
#wga <- strsplit(rundir,split="_")[[1]][2]
#species <- strsplit(rundir,split="_")[[1]][3]
dirtable <- read.table(dir_names)

dirtable_vector <- dirtable[,1]
indiv_names <- c()
for(j in dirtable_vector){
  indiv_names <- c(indiv_names, strsplit(basename(j), split="_")[[1]][4])
}
indiv_names <- unique(indiv_names)

bicdirs <- list.dirs(paste(fin,"/..",sep=""), recursive = F)
bicdirs_final <- c()
#Actually, we don't want to look at other directories, so make bicdirs empty and skip this
bicdirs <- c()
for(u in 1:length(bicdirs)){
  #print(paste(bicdirs[u],"/bin_stuff",sep=""))
  add_bicdir <- list.dirs(path=paste(bicdirs[u],"/bic_outputs",sep=""))
  if(length(add_bicdir) > 0){
    bicdirs_final <- c(bicdirs_final, add_bicdir)
  }
}

#bicdirs_final <- grep(glob2rx(trim.head = F, trim.tail = F), bicdirs_final,value = T) 

# bicdirs_final <- bicdirs_final <- bicdirs_final[grep(wga,bicdirs_final)]
# bicdirs_final <- bicdirs_final <- bicdirs_final[grep(species,bicdirs_final)]

#print(bicdirs_final)

#allbic <- data.frame()
allbic <- read.table(paste(fin,"/bic_outputs/BICs_and_other_info_badbinsOUT.txt",sep=''), header=T, sep="\t")
bicdirs_final <- c()


# UNCOMMENT IF NEEDED

#for(i in 1:length(bicdirs_final)){
 # bicfilename <- list.files(bicdirs_final[i],pattern="BICs_and_other_info_badbinsOUT.txt",full.names = T) #this location will need to be changed
  #if(length(bicfilename) > 0) {
  #  bicfile <- read.table(bicfilename, header=T, sep="\t")
   # allbic <- rbind(allbic, bicfile)
  #}
#}

# allbic <- data.frame()
# for(i in 1:length(bictables)){
#   t <- read.table(bictables[i],header = T,sep="\t")
#   allbic <- rbind(allbic, t)
# }

# ggplot(allbic,aes(BIC, fill=Run_dir)) +
#   geom_histogram() +
#   facet_wrap(~Run_dir, ncol=1)
#
# ggplot(allbic,aes(BIC, fill=Run_dir)) +
#   geom_histogram()
if(nrow(allbic) <3 ){
  bic_cutoff <- -2
  print("Making BIC plot")
  bicplot <- ggplot(data=allbic, mapping=aes(as.numeric(BIC), fill=Individual)) +
    geom_histogram(binwidth = 0.05) +
    xlab("BIC score") + ylab("Cells") +
    ggtitle(paste("BIC cutoff <",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
    #  scale_fill_manual(values=unique(all$colors), name="Age") +
    #guides(fill=F) +
    guides(color=F) +geom_vline(xintercept = bic_cutoff,col="purple4", size=1.5, linetype="dashed")
#  ggsave(filename = "BIC_cutoff_plot.png", bicplot, height=8, width=8)
  ggsave(filename = paste(fin,"/bic_outputs/BIC_cutoff_plot_default_cutoff.png",sep=""), bicplot, height=8, width=12)

  print("Making individual's BIC plot")
  for(i in 1:length(indiv_names)){
  indivname <- indiv_names[i]
  bicplot_indiv <- ggplot(data=allbic %>% filter(Individual==indivname),
                          mapping=aes(BIC)) +
    geom_histogram(binwidth = 0.05,fill="gray60") +
    xlab("BIC score") + ylab("Cells") +
    ggtitle(paste("BIC cutoff <",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
    #  scale_fill_manual(values=unique(all$colors), name="Age") +
    #guides(fill=F) +
    #guides(color=F) +
    #Set axes limits to previous
    xlim(layer_scales(bicplot)$x$range$range[1],layer_scales(bicplot)$x$range$range[2]) +
    ylim(layer_scales(bicplot)$y$range$range[1],layer_scales(bicplot)$y$range$range[2]) +
    geom_vline(xintercept = bic_cutoff,
               col="purple4", size=1.5, linetype="dashed")
  #  ggsave(filename = "BIC_cutoff_plot_indiv.png", bicplot_indiv, height=8, width=12)
  ggsave(filename = paste(fin,"/bic_outputs/BIC_cutoff_plot_",indivname,"_default_cutoff.png",sep=""), bicplot_indiv, height=8, width=12) 
  indivbic <- allbic %>% filter(Individual==indivname)
  write.table(indivbic, paste(fin,"/bic_outputs/",indivname,"_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''), sep="\t",quote=F, row.names=F)
  }

  allbic$BIC_cutoff <- sprintf("%.3f",round(bic_cutoff,3))
  allbic$BIC_passfail <- allbic$BIC < bic_cutoff
  allbic$BIC_passfail[which(allbic$BIC_passfail == TRUE)] <- "Pass_BIC"
  allbic$BIC_passfail[which(allbic$BIC_passfail == FALSE)] <- "Fail_BIC"

  write.table(allbic, paste(fin,"/bic_outputs/all_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''), sep="\t",quote=F, row.names=F)
}else if(num_gaussians > 1){
 ad_BIC <- normalmixEM(as.numeric(allbic$BIC), k=num_gaussians)
  #                       mu = c(-2.22,-1.68), lambda=c(0.9,0.1))

  plot(ad_BIC, whichplots=2, breaks = 100, cex.axis=1.4, cex.lab=1.4,
       cex.main=1.8, main2="BIC score distributions", xlab2="BIC score")

  big_gaussian <- which(ad_BIC$lambda==max(ad_BIC$lambda))

  bic_cutoff <- qnorm(0.95, ad_BIC$mu[big_gaussian], ad_BIC$sigma[big_gaussian])
  bic_cutoff <- round(bic_cutoff,3)

  calc.components <- function(x, mix, comp.number) {
    mix$lambda[comp.number]*dnorm(x, mean=mix$mu[comp.number],mix$sigma[comp.number])*length(allbic$BIC)*0.05
  }
  distrib_colors <- c("black","red1")

  numComponents=2
  distComps <- lapply(seq(numComponents), function(i)
    stat_function(fun=calc.components,
                  args=list(mix=ad_BIC, comp.number=i),
                  geom="line",
                  size=1.5,
                  color=distrib_colors[i]))

  print("Making BIC plot")
  bicplot <- ggplot(data=allbic, mapping=aes(as.numeric(BIC), fill=Individual)) +
    geom_histogram(binwidth = 0.05) +
    xlab("BIC score") + ylab("Cells") +
    ggtitle(paste("BIC cutoff <",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
    #  scale_fill_manual(values=unique(all$colors), name="Age") +
    #guides(fill=F) +
    guides(color=F) +
    distComps +
    geom_vline(xintercept = bic_cutoff,
               col="purple4", size=1.5, linetype="dashed")
#  ggsave(filename = "BIC_cutoff_plot.png", bicplot, height=8, width=8)
  ggsave(filename = paste(fin,"/bic_outputs/BIC_cutoff_plot_",num_gaussians,"gaussians.png",sep=""), bicplot, height=8, width=12)

  print("Making individual's BIC plot")
  for(i in 1:length(indiv_names)){
  indivname <- indiv_names[i]
  bicplot_indiv <- ggplot(data=allbic %>% filter(Individual==indivname),
                          mapping=aes(BIC)) +
    geom_histogram(binwidth = 0.05,fill="gray60") +
    xlab("BIC score") + ylab("Cells") +
    ggtitle(paste("BIC cutoff <",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
    #  scale_fill_manual(values=unique(all$colors), name="Age") +
    #guides(fill=F) +
    #guides(color=F) +
    distComps +
    #Set axes limits to previous
    xlim(layer_scales(bicplot)$x$range$range[1],layer_scales(bicplot)$x$range$range[2]) +
    ylim(layer_scales(bicplot)$y$range$range[1],layer_scales(bicplot)$y$range$range[2]) +
    geom_vline(xintercept = bic_cutoff,
               col="purple4", size=1.5, linetype="dashed")
  #  ggsave(filename = "BIC_cutoff_plot_indiv.png", bicplot_indiv, height=8, width=12)
  ggsave(filename = paste(fin,"/bic_outputs/BIC_cutoff_plot_",indivname,"_",num_gaussians,"gaussians.png",sep=""), bicplot_indiv, height=8, width=12)
  indivbic <- allbic %>% filter(Individual==indivname)
  write.table(indivbic, paste(fin,"/bic_outputs/",indivname,"_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''), sep="\t",quote=F, row.names=F)
  }

  allbic$BIC_cutoff <- sprintf("%.3f",round(bic_cutoff,3))
  allbic$BIC_passfail <- allbic$BIC < bic_cutoff
  allbic$BIC_passfail[which(allbic$BIC_passfail == TRUE)] <- "Pass_BIC"
  allbic$BIC_passfail[which(allbic$BIC_passfail == FALSE)] <- "Fail_BIC"

  write.table(allbic, paste(fin,"/bic_outputs/all_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''), sep="\t",quote=F, row.names=F, append=T)
  

} else {
  bic_cutoff <- qnorm(0.95, mean(allbic$BIC), sd(allbic$BIC))
  print("Making BIC plot")
  bicplot <- ggplot(data=allbic, mapping=aes(BIC, fill=Individual)) +
    geom_histogram(binwidth = 0.05) +
    xlab("BIC score") + ylab("Cells") +
    ggtitle(paste("BIC cutoff <",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
    stat_function(fun=function(x,mean,sd,n_obs,binwidth){
      dnorm(x=x, mean=mean, sd=sd) * n_obs * binwidth
    },  args=c(mean=mean(allbic$BIC),
               sd=sd(allbic$BIC),
               n_obs=length(allbic$BIC),
               binwidth=0.05), col="black", size=1.5) +
    geom_vline(xintercept = bic_cutoff,
#                 qnorm(0.95, ad_BIC$mu[1],
#                                  ad_BIC$sigma[1]),
               col="purple4", size=1.5, linetype="dashed")
#  ggsave(filename = "BIC_cutoff_plot_1gaussian.png", bicplot, height=8, width=12)
  ggsave(filename = paste(fin,"/bic_outputs/BIC_cutoff_plot_1gaussian.png",sep=""), bicplot, height=8, width=8)

 for(i in 1:length(indiv_names)){
  print("Making individual's BIC plot")
  indivname <- indiv_names[i]
  bicplot_indiv <- ggplot(data=allbic  %>% filter(Individual==indivname),
                          mapping=aes(BIC)) +
    geom_histogram(binwidth = 0.05, fill="gray60") +
    xlab("BIC score") + ylab("Cells") +
    ggtitle(paste("BIC cutoff <",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
    #  scale_fill_manual(values=unique(all$colors), name="Age") +
    guides(color=F) +
    stat_function(fun=function(x,mean,sd,n_obs,binwidth){
      dnorm(x=x, mean=mean, sd=sd) * n_obs * binwidth
    },  args=c(mean=mean(allbic$BIC),
               sd=sd(allbic$BIC),
               n_obs=length(allbic$BIC),
               binwidth=0.05), col="black", size=1.5) +
    #Set axes limits to previous
    xlim(layer_scales(bicplot)$x$range$range[1],layer_scales(bicplot)$x$range$range[2]) +
    ylim(layer_scales(bicplot)$y$range$range[1],layer_scales(bicplot)$y$range$range[2]) +
    geom_vline(xintercept = bic_cutoff,
               #                 qnorm(0.95, ad_BIC$mu[1],
               #                                  ad_BIC$sigma[1]),
               col="purple4", size=1.5, linetype="dashed")
  #  ggsave(filename = "BIC_cutoff_plot_1gaussian.png", bicplot, height=8, width=12)
  ggsave(filename = paste(fin,"/bic_outputs/BIC_cutoff_plot_",indivname,"_1gaussian.png",sep=""), bicplot_indiv, height=8, width=8) 
  indivbic <- allbic %>% filter(Individual==indivname)
  write.table(indivbic, paste(fin,"/bic_outputs/",indivname,"_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''), sep="\t",quote=F, row.names=F) 
 }

  allbic$BIC_cutoff <- sprintf("%.3f",round(bic_cutoff,3))
  allbic$BIC_passfail <- allbic$BIC < bic_cutoff
  allbic$BIC_passfail[which(allbic$BIC_passfail == TRUE)] <- "Pass_BIC"
  allbic$BIC_passfail[which(allbic$BIC_passfail == FALSE)] <- "Fail_BIC"

  write.table(allbic,paste(fin,"/bic_outputs/all_BICs_and_other_info_badbinsOUT_passfail.txt",sep=''), sep="\t",quote=F, row.names=F,append=T)
  

}

#Save new versions of bic tables with pass/fail column added
#One for specific run, one for whole wga

#Now
#Get CNV cutoffs

#Determine size of smallest chromosome (in bins)
example_longfile <- list.files(paste(fin,"/full_longs_badbinsOUT",sep=''), full.names=T)[1] 
lftable <- read.table(example_longfile, header=T, sep="\t")
lftable_noNA <- lftable %>% filter(!(is.na(est_copy)))
chr_bins <- lftable_noNA %>% group_by(Chrom) %>% summarize(bins=n())
smallest_autosome <- min(chr_bins$bins[which(chr_bins$Chrom < (max(chr_bins$Chrom)-1))])
max_subchrom_seg <- smallest_autosome - 1
segdirs <- list.dirs(paste(fin,"/..",sep=""), recursive = F)
segdirs_final <- c()
segdirs <- c()
for(u in 1:length(segdirs)){
  #print(paste(segdirs[u],"/segment_data",sep=""))
  add_segdir <- list.dirs(path=paste(segdirs[u],"/segment_data",sep=""))
  if(length(add_segdir) > 0){
    segdirs_final <- c(segdirs_final, add_segdir)
  }
}

#segdirs_final <- grep(glob2rx(trim.head = F, trim.tail = F), segdirs_final,value = T)

# segdirs_final <- segdirs_final <- segdirs_final[grep(wga,segdirs_final)]
# segdirs_final <- segdirs_final <- segdirs_final[grep(species,segdirs_final)]

#print(segdirs_final)

allseg <- data.frame()
allseg <-read.table(paste(fin,"/segment_data/all_segment_data_badbinsOUT.txt",sep=''), header=T, sep="\t")
segdirs_final <-c()
#UNCOMMENT IF NEEDED
#for(i in 1:length(segdirs_final)){
 # segfilename <- list.files(segdirs_final[i],pattern="all_segment_data_badbinsOUT.txt",full.names = T)
 # if(length(segfilename) > 0) {
 #   segfile <- read.table(segfilename, header=T, sep="\t")
  #  allseg <- rbind(allseg, segfile)
 # }
#}

#for(i in 1:4){
#  tb <- read.table(list.files("../Segment_data","all_segment",full.names = T)[i], header=T, sep="\t")
#  allseg <- rbind(allseg, tb)
#}

allseg$BIC_cutoff <- sprintf("%.3f",round(bic_cutoff,3))
allseg$BIC_passfail <- allseg$BIC < bic_cutoff
allseg$BIC_passfail[which(allseg$BIC_passfail == TRUE)] <- "Pass_BIC"
allseg$BIC_passfail[which(allseg$BIC_passfail == FALSE)] <- "Fail_BIC"


#Exclude segments from cells that do not pass BIC cutoff (determined above)
#Or are less than 5 bins long
#Or are on chrX or chrY
#allseg_eligible <- allseg %>% filter(BIC < bic_cutoff,  num.mark > 4,  chrom < max(allseg$chrom)-1)
#ggplot(allseg_eligible, aes(x=num.mark)) +
#  geom_histogram(binwidth = 1)


#When defining CNV thresholds, only examine segments smaller than the smallest chromosome
#allseg <- allseg_eligible %>% filter(num.mark <= max_subchrom_seg)

# ggplot(allseg, aes(x=seg.median, fill=Individual)) +
#   geom_histogram(binwidth = 0.05) +
#   theme_bw(base_size = 18) +
#   xlim(0,4)

#mixseg <- normalmixEM(allseg$seg.median, mu = c(1,2,3), lambda=c(0.1,0.8,0.1), maxit = 5000)
#plot(mixseg, whichplots = 2, breaks=50, cex.axis=1.4, cex.lab=1.4,
#     cex.main=1.8, main2 = "Density curves for segments of 5-40 bins",
#     xlab2="Segment CN value")

#qnorm(0.95, mean(allseg$seg.median),
#      sd(allseg$seg.median))

del_cutoff <- qnorm(0.025, mean(allseg$seg.median),sd(allseg$seg.median))
dup_cutoff <- qnorm(0.975, mean(allseg$seg.median),sd(allseg$seg.median))

#Previously determined CNV cutoffs:
del_cutoff_lenient <- 1.34
dup_cutoff_lenient <- 2.60
del_cutoff_stringent <- 1.14
dup_cutoff_stringent <- 2.80

print("Making segment plot")
segplot <- ggplot(data=allseg, mapping=aes(seg.median, fill=Individual)) +
  geom_histogram(binwidth = 0.05) +
  xlab("Copy number") + ylab("Segments") +
  ggtitle(paste("New CNV Cutoffs <",sprintf("%.2f",round(del_cutoff,2))," & >",
                sprintf("%.2f",round(dup_cutoff,2)),
                ", BIC < ",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
  theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
  #  scale_fill_manual(values=unique(all$colors), name="Age") +
  guides(color=F) +
  stat_function(fun=function(x,mean,sd,n_obs,binwidth){
    dnorm(x=x, mean=mean, sd=sd) * n_obs * binwidth
  },  args=c(mean=mean(allseg$seg.median),
             sd=sd(allseg$seg.median),
             n_obs=length(allseg$seg.median),
             binwidth=0.05), col="black", size=1.5) +
  xlim(0,ceiling(max(allseg$seg.median))) +
  geom_vline(xintercept = dup_cutoff, col="red3", size=1.5, linetype="dashed") +
  geom_vline(xintercept = del_cutoff, col="red3", size=1.5, linetype="dashed") +
  geom_vline(xintercept = del_cutoff_lenient,
             col="navyblue", size=1.5, linetype="dotted") +
  geom_vline(xintercept = dup_cutoff_lenient,
             col="navyblue", size=1.5, linetype="dotted") +
  geom_vline(xintercept = del_cutoff_stringent,
             col="green4", size=1.5, linetype="dotted") +
  geom_vline(xintercept = dup_cutoff_stringent,
             col="green4", size=1.5, linetype="dotted")
#ggsave(filename = paste("./CNV_cutoffs_plot_",wga,"_CN2gaussian.png",sep=""), segplot, height=8, width=12)

ggsave(filename = paste(fin,"/cnv_results/CNV_cutoffs_plot_CN2gaussian.png",sep=""), segplot, height=8, width=12) 

#Output new tables with CNV segment info for new cutoffs
#One for specific run, one for whole wga
allcnvs <- allseg %>% filter(num.mark > 4, chrom < max(allseg$chrom)-1, seg.median > dup_cutoff | seg.median < del_cutoff)
allcnvs$cnv_type <- allcnvs$seg.median > dup_cutoff
allcnvs$cnv_type[which(allcnvs$cnv_type == TRUE)] <- "dup"
allcnvs$cnv_type[which(allcnvs$cnv_type == FALSE)] <- "del"
allcnvs$cutoffs_employed <- "new"
allcnvs$del_cutoff <- del_cutoff
allcnvs$dup_cutoff <- dup_cutoff
#Save table of all CNVs unfiltered by BIC

#Output new tables with CNV segment info for lenient cutoffs
#One for specific run, one for whole wga
allcnvs_lenient <- allseg %>% filter(num.mark > 4, chrom < max(allseg$chrom)-1, seg.median > dup_cutoff_lenient | seg.median < del_cutoff_lenient)
allcnvs_lenient$cnv_type <- allcnvs_lenient$seg.median > dup_cutoff_lenient
allcnvs_lenient$cnv_type[which(allcnvs_lenient$cnv_type == TRUE)] <- "dup"
allcnvs_lenient$cnv_type[which(allcnvs_lenient$cnv_type == FALSE)] <- "del"
allcnvs_lenient$cutoffs_employed <- "lenient"
allcnvs_lenient$del_cutoff <- del_cutoff_lenient
allcnvs_lenient$dup_cutoff <- dup_cutoff_lenient
#Save table of all CNVs unfiltered by BIC

#Output new tables with CNV segment info for stringent cutoffs
#One for specific run, one for whole wga
allcnvs_stringent <- allseg %>% filter(num.mark > 4, chrom < max(allseg$chrom)-1, seg.median > dup_cutoff_stringent | seg.median < del_cutoff_stringent)
allcnvs_stringent$cnv_type <- allcnvs_stringent$seg.median > dup_cutoff_stringent
allcnvs_stringent$cnv_type[which(allcnvs_stringent$cnv_type == TRUE)] <- "dup"
allcnvs_stringent$cnv_type[which(allcnvs_stringent$cnv_type == FALSE)] <- "del"
allcnvs_stringent$cutoffs_employed <- "stringent"
allcnvs$del_cutoff <- del_cutoff_stringent
allcnvs$dup_cutoff <- dup_cutoff_stringent
#Save table of all CNVs unfiltered by BIC


print("Making individual segment plot")
for(i in 1:length(indiv_names)){
indivname <- indiv_names[i]
segplot_indiv <- ggplot(data=allseg %>% filter(Individual==indivname), mapping=aes(seg.median)) +
  geom_histogram(binwidth = 0.05, fill="gray60") +
  xlab("Copy number") + ylab("Segments") +
  ggtitle(paste("New CNV Cutoffs <",sprintf("%.2f",round(del_cutoff,2))," & >",
                sprintf("%.2f",round(dup_cutoff,2))," ",indivname,
                ", BIC < ",sprintf("%.3f",round(bic_cutoff,3)),sep="")) +
  theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) +
  #  scale_fill_manual(values=unique(all$colors), name="Age") +
  guides(color=F) +
  stat_function(fun=function(x,mean,sd,n_obs,binwidth){
    dnorm(x=x, mean=mean, sd=sd) * n_obs * binwidth
  },  args=c(mean=mean(allseg$seg.median),
             sd=sd(allseg$seg.median),
             n_obs=length(allseg$seg.median),
             binwidth=0.05), col="black", size=1.5) +
  xlim(layer_scales(segplot)$x$range$range[1],layer_scales(segplot)$x$range$range[2]) +
  ylim(layer_scales(segplot)$y$range$range[1],layer_scales(segplot)$y$range$range[2]) +
  geom_vline(xintercept = qnorm(0.975, mean(allseg$seg.median),
                                sd(allseg$seg.median)),
             col="red3", size=1.5, linetype="dashed") +
  geom_vline(xintercept = qnorm(0.025, mean(allseg$seg.median),
                                sd(allseg$seg.median)),
             col="red3", size=1.5, linetype="dashed") +
  geom_vline(xintercept = del_cutoff_lenient,
             col="navyblue", size=1.5, linetype="dotted") +
  geom_vline(xintercept = dup_cutoff_lenient,
             col="navyblue", size=1.5, linetype="dotted") +
  geom_vline(xintercept = del_cutoff_stringent,
             col="green4", size=1.5, linetype="dotted") +
  geom_vline(xintercept = dup_cutoff_stringent,
             col="green4", size=1.5, linetype="dotted")
#ggsave(filename = paste("./CNV_cutoffs_plot_",wga,"_",indivname,"_CN2gaussian.png",sep=""), segplot_indiv, height=8, width=12)
ggsave(filename = paste(fin,"/cnv_results/CNV_cutoffs_plot_",indivname,"_CN2gaussian.png",sep=""), segplot_indiv, height=8, width=12)

allcnvs_indiv <- allcnvs %>% filter(Individual==indivname)

write.table(allcnvs_indiv, paste(fin,"/cnv_results/new_thresholds/cnvs_",indivname,"_NObicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)

allcnvs_bicpass_indiv <- allcnvs %>% filter(Individual==indivname, BIC < bic_cutoff)

write.table(allcnvs_bicpass_indiv, paste(fin,"/cnv_results/new_thresholds/cnvs_",indivname,"_bicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)

allcnvs_lenient_indiv <- allcnvs_lenient %>% filter(Individual==indivname)
#Save table of indiv's CNVs unfiltered by BIC
write.table(allcnvs_lenient_indiv, paste(fin,"/cnv_results/lenient_thresholds/cnvs_",indivname,"_NObicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)

allcnvs_lenient_bicpass_indiv <- allcnvs_lenient %>% filter(Individual==indivname, BIC < bic_cutoff)
#Save table of indiv's CNVs filtered by BIC
write.table(allcnvs_lenient_bicpass_indiv, paste(fin,"/cnv_results/lenient_thresholds/cnvs_",indivname,"_bicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)
allcnvs_stringent_bicpass_indiv <- allcnvs_stringent %>% filter(Individual==indivname, BIC < bic_cutoff)
#Save table of indiv's CNVs filtered by BIC
allcnvs_stringent_indiv <- allcnvs_stringent %>% filter(Individual==indivname)
#Save table of indiv's CNVs unfiltered by BIC
write.table(allcnvs_stringent_indiv, paste(fin,"/cnv_results/stringent_thresholds/cnvs_",indivname,"_NObicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)
write.table(allcnvs_stringent_bicpass_indiv, paste(fin,"/cnv_results/stringent_thresholds/cnvs_",indivname,"_bicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)
}



write.table(allcnvs, paste(fin,"/cnv_results/new_thresholds/cnvs_NObicfilter.txt",sep=""), sep="\t",quote=F, row.names=F, append=T)

#Save table of indiv's CNVs unfiltered by BIC
allcnvs_bicpass <- allcnvs %>% filter(BIC < bic_cutoff)
#Save table of all CNVs filtered by BIC
write.table(allcnvs_bicpass, paste(fin,"/cnv_results/new_thresholds/cnvs_bicfilter.txt",sep=""), sep="\t",quote=F, row.names=F, append=T)

#Save table of indiv's CNVs filtered by BIC




write.table(allcnvs_lenient, paste(fin,"/cnv_results/lenient_thresholds/cnvs_NObicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)


allcnvs_lenient_bicpass <- allcnvs_lenient %>% filter(BIC < bic_cutoff)
#Save table of all CNVs filtered by BIC
write.table(allcnvs_lenient_bicpass, paste(fin,"/cnv_results/lenient_thresholds/cnvs_bicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)


write.table(allcnvs_stringent, paste(fin,"/cnv_results/stringent_thresholds/cnvs_NObicfilter.txt",sep=""), sep="\t",quote=F, row.names=F) 

allcnvs_stringent_bicpass <- allcnvs_stringent %>% filter(BIC < bic_cutoff)
#Save table of all CNVs filtered by BIC
write.table(allcnvs_stringent_bicpass, paste(fin,"/cnv_results/stringent_thresholds/cnvs_bicfilter.txt",sep=""), sep="\t",quote=F, row.names=F)

print("Step 3 of the pipeline has been completed")