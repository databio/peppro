#! /usr/bin/env Rscript
#
# PEPPRO_summarizer.R
#
# Interface to produce project level summary files and reports for nascent RNA
# profiling output when called using `looper`.
#
# Authors@R: as.person(c(
#    "Jason Smith <jasonsmith@virginia.edu> [aut, cre]" ))
#
# Created: 01/20/2020 Last updated: 04/27/2020
#
# usage: Rscript /path/to/Rscript/PEPPRO_summarizer.R
#        /path/to/project_config.yaml Depends: R (>= 3.5.1) Imports: PEPPROr, argparser
###############################################################################
##### LOAD ARGUMENTPARSER #####
loadLibrary <- tryCatch (
    {
        suppressWarnings(suppressPackageStartupMessages(library(argparser)))
    },
    error=function(e) {
        message("Error: Install the \"argparser\"",
                " library before proceeding.")
        return(NULL)
    },
    warning=function(e) {
        message(e)
        return(TRUE)
    }
)
if (length(loadLibrary)!=0) {
    suppressWarnings(library(argparser))
} else {
    quit()
}
# Create a parser
p <- arg_parser("Produce nascent RNA profiling Summary Reports, Files, and Plots")
# Add command line arguments
p <- add_argument(p, "config", help="PEPPRO project_config.yaml")
# Parse the command line arguments
argv <- parse_args(p)

###############################################################################
##### LOAD DEPENDENCIES #####
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
##### MAIN #####

# Identify the project configuration file
pep <- argv$config
prj <- invisible(suppressWarnings(pepr::Project(pep)))

# Project genomes
genomes <- invisible(suppressWarnings(pepr::sampleTable(prj)$genome))

# Produce output directory (if needed)
output_dir <- suppressMessages(
    file.path(pepr::config(prj)$looper$output_dir, "summary"))
#output_dir <- system(paste0("echo ", output_dir), intern = TRUE)
dir.create(output_dir, showWarnings = FALSE)

# Create project assets summary
assets <- PEPPROr::createAssetsSummary(pep)

# Produce library complexity summary plots
complexity_path <- PEPPROr::buildFilePath("_libComplexity.pdf", prj)
#complexity_path <- system(paste0("echo ", complexity_path), intern = TRUE)
if (!file.exists(complexity_path)) {
    cc <- paste(suppressMessages(pepr::config(prj)$looper$output_dir),
                "results_pipeline",
                suppressMessages(pepr::sampleTable(prj)$sample_name),
                paste0("QC_", suppressMessages(pepr::sampleTable(prj)$genome)),
                paste0(suppressMessages(pepr::sampleTable(prj)$sample_name),
                       "_preseq_yield.txt"),
                sep="/")
    #cc <- system(paste0("echo ", cc), intern = TRUE)
    rc <- paste(suppressMessages(pepr::config(prj)$looper$output_dir),
                "results_pipeline",
                suppressMessages(pepr::sampleTable(prj)$sample_name),
                paste0("QC_", suppressMessages(pepr::sampleTable(prj)$genome)),
                paste0(suppressMessages(pepr::sampleTable(prj)$sample_name),
                       "_preseq_counts.txt"),
                sep="/")
    #rc <- system(paste0("echo ", rc), intern = TRUE)
    hasBoth <- file.exists(cc) & file.exists(rc)
    ccSub <- cc[hasBoth]
    rcSub <- rc[hasBoth]
    message(paste0(length(ccSub), " of ", length(cc), " files available"))
    if (sum(hasBoth) > 0){
        p <- PEPPROr::plotComplexityCurves(ccurves = ccSub, coverage = 0,
                                            read_length = 0,
                                            real_counts_path = rcSub,
                                            ignore_unique = FALSE)
        output_file <- PEPPROr::buildFilePath("_libComplexity.pdf", prj)
        #output_file <- system(paste0("echo ", output_file), intern = TRUE)
        pdf(file = output_file, width= 10, height = 7, useDingbats=F)
        suppressWarnings(print(p))
        invisible(dev.off())
        output_file <- PEPPROr::buildFilePath("_libComplexity.png", prj)
        #output_file <- system(paste0("echo ", output_file), intern = TRUE)
        png(filename = output_file, width = 686, height = 480)
        suppressWarnings(print(p))
        invisible(dev.off())
    } else {
        complexity_path <- NULL
        complexity_flag <- FALSE
        message("No samples have available library complexity files.")
    }
    if (!is.null(complexity_path) && file.exists(complexity_path)) {
        complexity_flag <- TRUE
    }
} else {
    complexity_flag <- TRUE
}

counts_path <- PEPPROr::buildFilePath("_countData.csv", prj)
#counts_path <- system(paste0("echo ", counts_path), intern = TRUE)
if (!file.exists(counts_path)) {
    counts_file <- PEPPROr::calcCountsTable(prj)
    if (!is.null(counts_file)) {
        data.table::fwrite(counts_file, file=counts_path)
        # Export as csv with rownames as the first column 
        # Will require modification upon loading into R to convert to matrix
        message(paste0("Counts table: ", file.path(counts_path), "\n"))
    } else {
        message("No samples have available gene count data.")
    }
}
