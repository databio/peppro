#! /usr/bin/env Rscript
###############################################################################
#6/10/19
#Author: Jason Smith
#PEPPRO_summarizer.R
#
#This program is meant to plot multiple library complexity curves on the 
#same plot when called by looper summarize
#
#NOTES:
#usage: Rscript tools/PEPPRO_summarizer.R 
#       /path/to/project_config.yaml
#
#requirements: PEPPROr
#
###############################################################################
version <- 0.1
##### Load dependencies #####

required_libraries <- c("PEPPROr")
for (i in required_libraries) {
    loadLibrary <- tryCatch (
        {
            suppressPackageStartupMessages(
                suppressWarnings(library(i, character.only=TRUE)))
        },
        error=function(e) {
            message("Error: Install the \"", i,
                    "\" R package before proceeding.")
            message('i.e. devtools::install_github("databio/peppro",',
                    ' subdir="PEPPROr")')
            return(NULL)
        },
        warning=function(e) {
            message(e)
            return(1)
        }
    )
    if (length(loadLibrary)!=0) {
        suppressWarnings(library(i, character.only=TRUE))
    } else {
        quit()
    }
}

###############################################################################

configFile <- opt_get_verb()
prj        <- suppressWarnings(Project(configFile))
# Produce output directory (if needed)
dir.create(
    suppressMessages(
        file.path(config(prj)$metadata$output_dir, "summary")),
    showWarnings = FALSE)

# Plot combined library complexity curves for all samples in project
cc <- paste(config(prj)$metadata$output_dir,
            "results_pipeline",
            samples(prj)$sample_name,
            paste0("QC_", samples(prj)$genome),
            paste0(samples(prj)$sample_name, "_preseq_yield.txt"),
            sep="/")
rc <- paste(config(prj)$metadata$output_dir,
            "results_pipeline",
            samples(prj)$sample_name,
            paste0("QC_", samples(prj)$genome),
            paste0(samples(prj)$sample_name, "_preseq_counts.txt"),
            sep="/")
output_name = paste(config(prj)$metadata$output_dir, "summary",
                    paste0(config(prj)$name, "_libComplexity"), sep="/")
p <- plotComplexityCurves(ccurves = cc, coverage = 0, read_length = 0,
                          real_counts_path = rc, ignore_unique = FALSE)

pdf(file = paste0(tools::file_path_sans_ext(output_name), ".pdf"),
            width= 10, height = 7, useDingbats=F)
print(p)
invisible(dev.off())
png(filename = paste0(tools::file_path_sans_ext(output_name), ".png"),
    width = 686, height = 480)
print(p)
invisible(dev.off())

###############################################################################
