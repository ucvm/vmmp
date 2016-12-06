#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
rmarkdown::render(input = "QC_report.Rmd", output_file = args[1],
	    params = list(fastqc_data = args[2],
	     	          runID = args[3],
			  merge_log = args[4],
			  filter_log = args[5])
)
