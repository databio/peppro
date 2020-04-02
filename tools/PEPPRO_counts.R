#! /usr/bin/env Rscript
###############################################################################
#1/20/2020
#Author: Jason Smith
#PEPPRO_counts.R
#
#For each sample in a project, generate a counts object of gene counts
#for each gene for each sample in the 
#format [gene],[sample1],[sample...],[sampleN]
#
#NOTES:
#usage: Rscript tools/PEPPRO_counts.R 
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

readPepproGeneCounts = function(project) {
    cwd            <- getwd()
    project_dir    <- pepr::config(project)$looper$output_dir
    sample_names   <- pepr::sampleTable(project)$sample_name
    genomes        <- as.list(pepr::sampleTable(project)$genome)
    names(genomes) <- sample_names
    paths          <- vector("list", length(sample_names))
    names(paths)   <- sample_names

    for (sample in sample_names) {
        paths[[sample]] <- paste(project_dir, 'results_pipeline', sample,
                                 paste0('signal_', genomes[[sample]]),
                                 paste0(sample, "_gene_coverage.bed"), sep="/")
    }

    result <- lapply(paths, function(x){
        if (file.exists(x)) {
            df <- fread(x)
            colnames(df) <- c('chr', 'start', 'end', 'geneName',
                              'score', 'strand', 'count')
            gr <- GenomicRanges::GRanges(df) 
        } else {
            gr <- GenomicRanges::GRanges() 
        }
    })

    setwd(cwd)
    return(GenomicRanges::GRangesList(Filter(length, result)))
}

###############################################################################

configFile <- opt_get_verb()
prj        <- suppressWarnings(Project(configFile))

message("Creating counts table...")

# Produce output directory (if needed)
dir.create(
    suppressMessages(
        file.path(pepr::config(prj)$looper$output_dir, "summary")),
    showWarnings = FALSE)

# Load gene counts files
totalSamples  <- length(suppressMessages(pepr::sampleTable(prj)$sample_name))
countsGR      <- readPepproGeneCounts(prj)
actualSamples <- length(countsGR)

message(paste0(actualSamples, " of ", totalSamples, " files available"))

# Generate output file name
output_name <- paste(pepr::config(prj)$looper$output_dir, "summary",
                     paste0(pepr::config(prj)$name, "_countData.csv"), sep="/")

if (length(countsGR) > 0) {
    # Create gene name data table
    i <- 1
    while (length(countsGR[[i]]) == 0) {
      i <- i + 1
    }
    count_dt <- data.table(geneName = countsGR[[i]]$geneName,
                           seqnames=as.character(seqnames(countsGR[[i]])),
                           start=start(countsGR[[i]]),
                           end=end(countsGR[[i]]),
                           width=width(countsGR[[i]]),
                           strand=as.character(strand(countsGR[[i]])))

    # Populate count table
    while (i < length(names(countsGR))) {
        dt1   <- as.data.table(countsGR[[names(countsGR)[i]]])
        name1 <- paste0(".",names(countsGR)[i])
        i     <- i + 1
        dt2   <- as.data.table(countsGR[[names(countsGR)[i]]])
        name2 <- paste0(".",names(countsGR)[i])

        dt    <- merge(dt1[,-"score"], dt2[,-"score"],
                       by=c("seqnames", "start", "end",
                            "width", "strand", "geneName"),
                       sort=TRUE, suffix=c(name1, name2))

        count_dt <- merge(count_dt, dt, sort=TRUE,
                          by=c("geneName", "seqnames", "start", "end",
                               "width", "strand"))
        i     <- i + 1
    }

    colnames(count_dt) <- sub(".*\\.","",colnames(count_dt))
    counts             <- count_dt[,-c("seqnames", "start", "end",
                                       "width", "strand")]
    #rownames(counts)   <- count_dt$geneName
    # Export as csv with rownames as the first column
    # Will require modification upon loading into R to convert to matrix
    fwrite(counts, file=output_name)
    message(paste0("Counts table: ", file.path(output_name), "\n"))
} else {
    message("No samples have available gene count data.")
}

###############################################################################
