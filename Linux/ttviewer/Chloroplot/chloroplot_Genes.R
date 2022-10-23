#!/usr/bin/env Rscript

options (warn = -1)
library('ape')
library('dplyr')
library('coRdon')
library('circlize')
library('genbankr')
library('magrittr')

source('../Chloroplot/R/color_complement.R')
source('../Chloroplot/R/converse_ssc.R')
source('../Chloroplot/R/detect_ir.R')
source('../Chloroplot/R/GC_count.R')
source('../Chloroplot/R/gene_color.R')
source('../Chloroplot/R/gene_info.R')
source('../Chloroplot/R/parse_gb_file.R')
source('../Chloroplot/R/plot_genome.R')
source('../Chloroplot/R/read_gb_file.R')
source('../Chloroplot/R/test_parameters.R')

Args <- commandArgs(trailingOnly=TRUE)
gb_file = Args[1]
plot_table <- PlotTab(gb_file, local.file = TRUE)
PlotPlastidGenome(plot_table, file.name = Args[2])
print("---  Pdf file has been generated, please check the folder under the current file!  ---")
