# scripts
step1, step2, step3, and step4 are run in succession to carry out CNV analysis of single cell whole genome sequencing data. The other scripts are invoked in the process.

The following is a basic outline summarizing the analysis:

step1_prepare_and_submit_jobs.sh
  CNV_pipeline_FASTQorBAM_thru_BIN_COUNTS.py
  windowMaker.py
  dupeStats_MOD.py
  Normalize_read_counts_and_segment.R
  Bad_bin_finder.R
step2_remove_badbins_and_reanalyze.sh
  Remove_bad_bins_Segment_Plot.R
step3_find_BIC_CNV_cutoffs.sh
  BIC_n_CNV_cutoff_finder.R
step4_make_cnv_sawtooth_get_summary.sh
  organize_CNVs_for_bedtools_genomecov.R
  plot_genomecov_beds.R
  compute_final_stats.R

These scripts do require some setup to run, including but not limited to installing software tools and packages and editing the hardcoded filepaths.
