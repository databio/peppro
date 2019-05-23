#! /usr/bin/env Rscript
###############################################################################
#8/13/18
#Last updated: 05/22/19
#Author: Jason Smith
#frif.R
#
#This program is meant to plot the fraction of reads in regions of interest
#
#NOTES:
#usage: Rscript /path/to/Rscript/frif.R
#               sample name
#               /path/to/*_coverage.bed file(s)
#               total_mapped_reads 
#               /path/to/outputFile
#
#requirements: argparser, GenomicDistributions, ggplot2
#
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
argP <- arg_parser("Produce Fraction of Reads in Features (FRiF) plot(s)")

# Add command line arguments
argP <- add_argument(argP, "name",
                     help="Sample name")
argP <- add_argument(argP, "reads", 
                     help="Number of mapped reads")
argP <- add_argument(argP, "output", 
                     help="Output file")
argP <- add_argument(argP, "--bed", nargs=Inf,
                     help="Coverage file(s)")


# Parse the command line arguments
argv <- parse_args(argP)

###############################################################################

##### LOAD DEPENDENCIES #####
required_libraries <- c("GenomicDistributions", "ggplot2")
for (i in required_libraries) {
    loadLibrary <- tryCatch (
        {
            suppressPackageStartupMessages(
                suppressWarnings(library(i, character.only=TRUE)))
        },
        error=function(e) {
            message("Error: Install the \"", i,
                    "\" library before proceeding.")
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

#### FUNCTIONS ####
calcFRiF <- function(bedFile) {
    colnames(bedFile) <- c("chromosome","start","end","count")
    grObj   <- makeGRangesFromDataFrame(bedFile)
    grObj   <- reduce(grObj)
    redBed  <- data.frame(chromosome=seqnames(grObj),
                          start=start(grObj), end=end(grObj))
    bedFile <- merge(redBed, bedFile, by=c("chromosome","start","end"))
    bedFile <- cbind(bedFile, size=(bedFile$end-bedFile$start))
    bedFile <- bedFile[order(-bedFile$count),]
    bedFile <- bedFile[apply(bedFile != 0, 1, all),]
    bedFile <- cbind(bedFile, cumsum=cumsum(bedFile$count))
    bedFile <- cbind(bedFile, cumSize=cumsum(bedFile$size))
    bedFile <- cbind(bedFile, frip=bedFile$cumsum/as.numeric(argv$reads))
    bedFile <- cbind(bedFile, numfeats=as.numeric(1:nrow(bedFile)))
    return(bedFile)
}

###############################################################################

labels  <- data.frame(xPos=numeric(), yPos=numeric(), name=character(),
                      val=numeric(), color=character(), stringsAsFactors=FALSE)
palette <- colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                              "#B97BC8", "#009E73", "#C92404", "#E3E550",
                              "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                              "#897779", "#6114F8", "#19C42B", "#56B4E9"))
plotColors <- palette(length(argv$bed))
info    <- file.info(file.path(argv$bed[1]))

if (file.exists(file.path(argv$bed[1])) && info$size != 0) {
    bed        <- read.table(file.path(argv$bed[1]))
    bedCov     <- calcFRiF(bed)
    name       <- basename(tools::file_path_sans_ext(argv$bed[1]))
    name       <- gsub(argv$name, "", name)
    name       <- gsub("^.*?_", "", name)
    numFields  <- 2
    for(i in 1:numFields) name <- gsub("_[^_]*$", "", name)  
    labels[1,] <- c(0.95*max(log10(bedCov$cumSize)), max(bedCov$frip)+0.001,
                    name, round(max(bedCov$frip),2), "#FF0703")
    bedCov$feature <- name
}  else {
    if (info$size == 0) {
        message(paste0(name, " coverage file is empty"))
    } else {
        message(paste0(name, " coverage file is missing"))
    }
}

if (exists("bedCov")) {
    covDF <- bedCov
}

if (length(argv$bed) > 1) {
    for (i in 2:length(argv$bed)) {
        name       <- basename(tools::file_path_sans_ext(argv$bed[i]))
        name       <- gsub(argv$name, "", name)
        name       <- gsub("^.*?_", "", name)
        numFields  <- 2
        for(j in 1:numFields) name <- gsub("_[^_]*$", "", name) 

        if (file.exists(file.path(argv$bed[i])) && info$size != 0) {
            bed     <- read.table(file.path(argv$bed[i]))
        }  else {
            outFile <- file.path(argv$output)
            system2(paste("touch"), outFile)
            quit()
        }
        if (max(bed[,4] > 0)) {
            if (exists("covDF")) {
                covFile <- calcFRiF(bed)
                covFile$feature <- name
                covDF   <- rbind(covDF, covFile)
                labels  <- rbind(labels, c(0.95*max(log10(covFile$cumSize)),
                                           max(covFile$frip)+0.001,
                                           name, round(max(covFile$frip),2),
                                           plotColors[i]))           
            } else {
                covDF         <- calcFRiF(bed)
                covDF$feature <- name
                labels        <- rbind(labels, c(0.95*max(log10(covDF$cumSize)),
                                                 max(covDF$frip)+0.001,
                                                 name, round(max(covDF$frip),2),
                                                 plotColors[i]))
            }
        }
    }
}

# Reorder by labels
if (exists("covDF")) {
    covDF$feature <- factor(covDF$feature, levels=(labels$name))
}

if (!is.null(argv$bed)) {
    # Produce plot with bed files
    p <- ggplot(covDF, aes(x=log10(cumSize), y=frip,
                           group=feature, color=feature)) +
        geom_line() +
        labs(x="log(number of bases)", y="FRiF") +
        theme_classic() +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    
    # Recolor and reposition legend
    p <- p + scale_color_manual(labels=paste0(labels$name, ": ", labels$val),
                                values=labels$color) +
             theme(legend.position=c(0.05,0.95),
                   legend.justification=c(0.1,0.9))
} else {
    write("Unable to produce FRiF plot!\n", stdout())
}

if (!exists("p")) {
    p <- ggplot()
}

pdf(file = paste0(tools::file_path_sans_ext(argv$output), ".pdf"),
    width= 7, height = 7, useDingbats=F)
p
invisible(dev.off())
png(filename = paste0(tools::file_path_sans_ext(argv$output), ".png"),
    width = 480, height = 480)
p
invisible(dev.off())

if (exists("p")) {
    write("Cumulative FRiF plot completed!\n", stdout())
} else {
    write("Unable to produce FRiF plot!\n", stdout())
}

###############################################################################