#! /usr/bin/env Rscript
###############################################################################
#1/20/2020
#Author: Jason Smith
#PEPPRO_complexity_curves.R
#
#This program is meant to plot multiple library complexity curves on the 
#same plot when called by looper summarize
#
#NOTES:
#usage: Rscript tools/PEPPRO_complexity_curves.R 
#       /path/to/project_config.yaml
#
#requirements: PEPPROr
#
###############################################################################
version <- 0.2
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
        file.path(config(prj)$looper$output_dir, "summary")),
    showWarnings = FALSE)

# Plot combined library complexity curves for all samples in project
cc <- paste(suppressMessages(config(prj)$looper$output_dir),
            "results_pipeline",
            suppressMessages(sampleTable(prj)$sample_name),
            paste0("QC_", suppressMessages(sampleTable(prj)$genome)),
            paste0(suppressMessages(sampleTable(prj)$sample_name),
                   "_preseq_yield.txt"),
            sep="/")
rc <- paste(suppressMessages(config(prj)$looper$output_dir),
            "results_pipeline",
            suppressMessages(sampleTable(prj)$sample_name),
            paste0("QC_", suppressMessages(sampleTable(prj)$genome)),
            paste0(suppressMessages(sampleTable(prj)$sample_name),
                   "_preseq_counts.txt"),
            sep="/")

hasBoth <- file.exists(cc) & file.exists(rc)

ccSub <- cc[hasBoth]
rcSub <- rc[hasBoth]
#message(ccSub, rcSub)
message(paste0(length(ccSub), " of ", length(cc), " files available"))

output_name <- paste(config(prj)$looper$output_dir, "summary",
                     paste0(config(prj)$name, "_libComplexity"), sep="/")

if (sum(hasBoth) > 0){

    p <- plotComplexityCurves(ccurves = ccSub, coverage = 0, read_length = 0,
                              real_counts_path = rcSub, ignore_unique = FALSE)

    pdf(file = paste0(tools::file_path_sans_ext(output_name), ".pdf"),
        width= 10, height = 7, useDingbats=F)
    print(p)
    invisible(dev.off())
    png(filename = paste0(tools::file_path_sans_ext(output_name), ".png"),
        width = 686, height = 480)
    print(p)
    invisible(dev.off())
} else {
    message("No samples have available library complexity files.")
}

###############################################################################
