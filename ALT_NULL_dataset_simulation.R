#Script for simulating data under NULL and ALT models

#read in longtable containing bin values
longtable <- read.table("name_of_longtable")

#read in shorttable containing segment information from DNAcopy (segment means)
shorttable <- read.table("name_of_shorttable")

#Set number of autosomal bins and total bins (autosomal + XY)
number_of_autosomal_bins <- 4272
number_all_bins <- 4505


longtable_arstuff <- ar(longtable$bin_value[1:number_of_autosomal_bins], order.max = 1, aic=F)
longtable_phi <- longtable_arstuff$ar
longtable_var <- longtable_arstuff$var.pred

sim_null_values <- 2 + arima.sim(model = list(ar = longtable_phi), n = number_of_autosomal_bins, sd=sqrt(longtable_var))

#Add XY values at copy number 1 to complete dataset
full_sim_null_values <- c(sim_null_values, rep(1,number_all_bins-number_of_autosomal_bins))
bin_number <- c(1:number_all_bins)


## ALT MODEL (mean subtraction) ##
#Subtract mean from each bin value before computing phi (alternate model)
longtable_arstuff_meansubd <- ar((longtable$bin_value[1:number_of_autosomal_bins]-longtable$seg.mean[1:number_of_autosomal_bins]), order.max = 1, aic=F)
longtable_phi_meansubd <- longtable_arstuff_meansubd$ar
longtable_var_meansubd <- longtable_arstuff_meansubd$var.pred

#How many segments, excluding X and Y (chr 23 and 24)?
autosome_segments <- nrow(shorttable[which(shorttable$Chromosome < 23),])

sim_alt_values <- c()
for (row_num in 1:autosome_segments){
  #Add mean of segment to arima.sim'd values using overall phi (mean subtracted) and overall variance (mean subtracted)
  sim_alt_values_grp <- shorttable$seg.mean[row_num] + 
    arima.sim(model = list(ar = longtable_phi_meansubd), n = shorttable$bins_in_segment[row_num], sd=sqrt(longtable_var_meansubd))
  sim_alt_values <- c(sim_alt_values, sim_alt_values_grp)
}

#Add XY values at copy number 1 to complete dataset
full_sim_alt_values <- c(sim_alt_values, rep(1,number_all_bins-number_of_autosomal_bins))
bin_number <- c(1:number_all_bins)
