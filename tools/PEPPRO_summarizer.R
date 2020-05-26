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
p <- add_argument(p, arg="config", short="-c",
                  help="PEPPRO project_config.yaml")
p <- add_argument(p, arg="output", short="-o",
                  help="Project parent output directory path")
p <- add_argument(p, arg="results", short="-r",
                  help="Project results output subdirectory path")
p <- add_argument(p, arg="--new-start", short="-N", flag=TRUE,
                  help=paste0("New start mode. This flag will tell the ",
                       "summarizer to start over, and run every command, even ",
                       "if its target output already exists."))
# Parse the command line arguments
argv <- parse_args(p)
#print(argv)  # DEBUG

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

# Set the project configuration file
pep <- argv$config
# Load the project
prj <- invisible(suppressWarnings(pepr::Project(pep)))
# Convenience
project_name <- config(prj)$name

# Set the output directory
summary_dir <- suppressMessages(file.path(argv$output, "summary"))
# Produce output directory (if needed)
dir.create(summary_dir, showWarnings = FALSE)

# Set the results subdirectory
if (dir.exists(argv$results)) {
    results_subdir <- suppressMessages(file.path(argv$results))
} else {
    warning(paste0("The project results subdirectory (", argv$results,
            ") does not exist."))
    quit()
}

# Get project genomes
genomes <- invisible(suppressWarnings(pepr::sampleTable(prj)$genome))

# Create project assets summary
assets  <- PEPPROr::createAssetsSummary(prj, argv$output, results_subdir)

# Produce library complexity summary plots
complexity_path <- file.path(summary_dir,
                             paste0(project_name, "_libComplexity.pdf"))
if (!file.exists(complexity_path) || argv$new_start) {
    cc <- paste(results_subdir,
                suppressMessages(pepr::sampleTable(prj)$sample_name),
                paste0("QC_", suppressMessages(pepr::sampleTable(prj)$genome)),
                paste0(suppressMessages(pepr::sampleTable(prj)$sample_name),
                       "_preseq_yield.txt"),
                sep="/")
    rc <- paste(results_subdir,
                suppressMessages(pepr::sampleTable(prj)$sample_name),
                paste0("QC_", suppressMessages(pepr::sampleTable(prj)$genome)),
                paste0(suppressMessages(pepr::sampleTable(prj)$sample_name),
                       "_preseq_counts.txt"),
                sep="/")
    hasBoth <- file.exists(cc) & file.exists(rc)
    ccSub <- cc[hasBoth]
    rcSub <- rc[hasBoth]
    message(paste0(length(ccSub), " of ", length(cc),
            " library complexity files available."))
    if (sum(hasBoth) > 0){
        p <- PEPPROr::plotComplexityCurves(ccurves = ccSub, coverage = 0,
                                            read_length = 0,
                                            real_counts_path = rcSub,
                                            ignore_unique = FALSE)
        output_file <- PEPPROr::buildFilePath("_libComplexity.pdf", prj)
        pdf(file = output_file, width= 10, height = 7, useDingbats=F)
        suppressWarnings(print(p))
        invisible(dev.off())
        output_file <- PEPPROr::buildFilePath("_libComplexity.png", prj)
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

counts_path <- file.path(summary_dir, paste0(project_name, "_countData.csv"))
if (!file.exists(counts_path) || argv$new_start) {
    counts_file <- PEPPROr::calcCountsTable(prj, results_subdir)
    if (!is.null(counts_file)) {
        data.table::fwrite(counts_file, file=counts_path)
        # Export as csv with rownames as the first column 
        # Will require modification upon loading into R to convert to matrix
        message(paste0("Counts table: ", file.path(counts_path), "\n"))
    } else {
        message("No samples have available gene count data.")
    }
}
