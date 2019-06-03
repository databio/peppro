# PACKAGE DOCUMENTATION
#' PEPPRO R tools (PEPPROr)
#'
#'
#' PEPPROr is a collection of functions to generate QC and analysis plots
#' for PRO-(or GRO-)seq data analysis.
#'
#' @docType package
#' @name PEPPROr
#' @author Jason Smith
#'
#' @references \url{http://github.com/databio}
NULL

################################################################################
# FUNCTIONS

#' Plot library complexity curves
#'
#' This function plots library complexity curves using data from
#' preseq and produces pdf/png output files.
#' Adapted from preseq_complexity_curves.py from
#' https://github.com/ewels/ngi_visualizations/tree/master/ngi_visualizations/preseq_complexity_curves
#'
#' @param ccurves A single preseq output file from one sample.
#' @param more_ccurves A list of additional preseq output files from more samples
#' @param coverage Use coverage on axes instead of read counts. Enter the number
#'                 of base pairs of your reference genome.
#' @param read_length Sequence read length, for use in coverage calculations.
#' @param real_counts_path File name for a file with three columns -
#'                         preseq filename, total number of reads, number of
#'                         unique reads
#'                         (unique is optional, file is whitespace delimited)
#' @param use_unique If FALSE, ignore information about unique read counts
#'                   found in real_counts_path file.
#' @param output_name Desired output name (produces both .pdf and .png files).
#' @param x_min Lower x-limit (default 0)
#' @param x_max Upper x-limit (default 500 million)
#' @keywords preseq library complexity
#' @export
#' @examples
#' data("ccurve")
#' data("counts")
#' plot_complexity_curves(ccurve, coverage=3099922541, read_length=30,
#'                        real_counts_path=counts)
plot_complexity_curves <- function(ccurves,
                                   coverage=0, read_length=0,
                                   real_counts_path=FALSE, use_unique=TRUE,
                                   output_name='complexity_curves',
                                   x_min=0, x_max=500000000) {

    if (x_min < 0 || x_max <= x_min) {
        message(paste0("problem with x-min or x-max (", x_min, " ", x_max,
                       "). x-min must be >= 0 and < x-max"))
        quit()
    }

    # Convert limit counts to coverage
    if (coverage > 0) {
        if (read_length == 0) {
            message("Error: --coverage specified but not --read_length")
            quit()
        } else {
            coverage <- as.numeric(coverage) / as.numeric(read_length)
        }
        x_max    <- as.numeric(x_max) / coverage
        x_min    <- as.numeric(x_min) / coverage
    }

    # Get the real counts if we have them
    real_counts_total  <- character()
    real_counts_unique <- character()
    real_counts_name   <- character()
    if ("data.frame" %in% class(real_counts_path)) {
        rc_file            <- real_counts_path
        real_counts_name   <- basename(rc_file$V1)
        real_counts_total  <- as.integer(rc_file$V2)
        if (ncol(rc_file) == 3 && use_unique) {
            real_counts_unique <- as.integer(rc_file$V3)
        }
    } else if (!is.na(real_counts_path)) {
        info <- file.info(file.path(real_counts_path))
        if (file.exists(real_counts_path) && info$size != 0) {
            rc_file            <- fread(real_counts_path)
            real_counts_name   <- basename(rc_file$V1)
            real_counts_total  <- as.integer(rc_file$V2)
            if (ncol(rc_file) == 3 && use_unique) {
                real_counts_unique <- as.integer(rc_file$V3)
            }
        } else {
            message(paste0("Error loading real counts file: ", real_counts_path))
            if (!file.exists(real_counts_path)) {
                message("File could not be found.")
            } else if (info$size == 0) {
                message("File is empty.")
            }
            quit()
        }
    }

    # Convert real counts to coverage
    if (coverage > 0) {
        real_counts_total  <- as.numeric(real_counts_total)  / coverage
        real_counts_unique <- as.numeric(real_counts_unique) / coverage
    }

    # Set up plot params
    global_x_max_ccurve_limit <- 0
    global_y_max_ccurve_limit <- 0
    fig <- ggplot()
    max_label_length <- 0

    # Each ccurve will get a different color
    palette <- colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))

    clist <- data.table(TOTAL_READS = list(), EXPECTED_DISTINCT = list())
    # Go through inputs and plot line
    # if (!is.na(more_ccurves)) {
    #     ccurves <- list(ccurves, more_ccurves)
    # } else {
    #     ccurves <- list(ccurves)
    # }
    ccurves  <- list(ccurves)
    colormap <- palette(length(ccurves))
    for (c in 1:length(ccurves)) {
        message(paste0("c: ", c))
        name        <- basename(tools::file_path_sans_ext(ccurves[[c]]))
        numFields   <- 2
        for(j in 1:numFields) name <- gsub("_[^_]*$", "", name)
        sample_name <- name
        message(paste0("Processing ", sample_name))
        ctable <- fread(ccurves[[c]])
        if (coverage > 0) {
            if ("TOTAL_READS" %in% colnames(ctable)) {
                ccurve_TOTAL_READS <- as.numeric(ctable$TOTAL_READS) / coverage
            } else if ("total_reads" %in% colnames(ctable)) {
                ccurve_TOTAL_READS <- as.numeric(ctable$total_reads) / coverage
            }
            if ("EXPECTED_DISTINCT" %in% colnames(ctable)) {
                ccurve_EXPECTED_DISTINCT <- (
                    as.numeric(ctable$EXPECTED_DISTINCT) / coverage)
            } else if ("distinct_reads" %in% colnames(ctable)) {
                ccurve_EXPECTED_DISTINCT <- (
                    as.numeric(ctable$distinct_reads) / coverage)
            } else {
                messsage(paste0("Error, table ", c, " is not in the expected "))
                message("format... has it been generated with preseq?")
                quit()
            }
        } else {
            if ("TOTAL_READS" %in% colnames(ctable)) {
                ccurve_TOTAL_READS <- ctable$TOTAL_READS
            } else if ("total_reads" %in% colnames(ctable)) {
                ccurve_TOTAL_READS <- ctable$total_reads
            }
            if ("EXPECTED_DISTINCT" %in% colnames(ctable)) {
                ccurve_EXPECTED_DISTINCT <- ctable$EXPECTED_DISTINCT
            } else if ("distinct_reads" %in% colnames(ctable)) {
                ccurve_EXPECTED_DISTINCT <- ctable$distinct_reads
            } else {
                messsage(paste0("Error, table ", c, " is not in the expected "))
                message("format... has it been generated with preseq?")
                quit()
            }
        }
        if (c == 1) {
            clist <- data.table(
                TOTAL_READS = list(ccurve_TOTAL_READS),
                EXPECTED_DISTINCT = list(ccurve_EXPECTED_DISTINCT)
            )
        } else {
            clist <- rbindlist(list(clist,
                                    data.table(TOTAL_READS = list(ccurve_TOTAL_READS),
                                               EXPECTED_DISTINCT = list(ccurve_EXPECTED_DISTINCT))),
                               use.names=TRUE)
        }
        # ggplot(ctable, aes(rep(1:10000), EXPECTED_DISTINCT)) +
        #     geom_point() +
        #     geom_point(aes(rep(1:10000), TOTAL_READS, col='red')) +
        #     scale_y_continuous(limits=c(0,75000000)) +
        #     scale_x_continuous(limits=c(0,1000))

        # ggplot(ctable, aes(rep(1:10000), EXPECTED_DISTINCT/coverage)) +
        #     geom_point() +
        #     geom_point(aes(rep(1:10000), TOTAL_READS/coverage, col='red')) +
        #     scale_y_continuous(limits=c(0,75000000/coverage)) +
        #     scale_x_continuous(limits=c(1,501))

        x_min_ccurve_limit <- computeLimit(x_min, ccurve_TOTAL_READS)
        x_max_ccurve_limit <- computeLimit(x_max, ccurve_TOTAL_READS)
        if (x_max_ccurve_limit > global_x_max_ccurve_limit) {
            global_x_max_ccurve_limit <- x_max_ccurve_limit
        }
        if (ccurve_EXPECTED_DISTINCT[x_max_ccurve_limit] > global_y_max_ccurve_limit) {
            if (x_max_ccurve_limit <= length(ccurve_EXPECTED_DISTINCT)) {
                global_y_max_ccurve_limit <- ccurve_EXPECTED_DISTINCT[x_max_ccurve_limit]
            } else {
                x_max_ccurve_limit <- length(ccurve_EXPECTED_DISTINCT)
            }
        }
        # Add a few points to be sure
        x_max_ccurve_limit <- x_max_ccurve_limit + 3

        # Plot the curve
        fig <- ggplot() +
            geom_line(aes(
                ccurve_TOTAL_READS[x_min_ccurve_limit:x_max_ccurve_limit],
                ccurve_EXPECTED_DISTINCT[x_min_ccurve_limit:x_max_ccurve_limit],
                col=colormap[c])
            )

        # Plot the real data if we have it
        if (length(real_counts_total) > 0 && length(real_counts_unique) > 0) {
            fig <- fig +
                geom_point(aes(real_counts_total,
                               real_counts_unique,
                               col=colormap[c]),
                           shape=23, size=3)
            message(paste0("INFO: Found real counts for ", sample_name,
                           " - Total: ", real_counts_total, " Unique: ",
                           real_counts_unique))
        } else if (length(real_counts_total) > 0) {
            ggp <- ggplot_build(fig)
            xvalues <- ggp$layout$panel_scales_x[[1]]$range$range
            yvalues <- ggp$layout$panel_scales_y[[1]]$range$range
            if (real_counts_total > max(xvalues)) {
                message(paste0("WARNING: Total reads for ",  sample_name,
                               "(", real_counts_total,
                               ") > max preseq value (", max(xvalues),
                               ") - skipping this point..."))
            }
            else {
                interp <- approx(xvalues, yvalues, real_counts_total)$y
                fig <- fig + geom_point(aes(real_counts_total,
                                            interp, col=colormap[c]))
                message(paste0("INFO: Found real count for ",
                               sample_name, " - Total: ",
                               real_counts_total, "(preseq unique reads: ",
                               interp, ")"))
            }
        } else {
            message(paste0("INFO: No real counts found for ", sample_name))
        }

        # plot perfect library as dashed line
        fig <- fig + geom_segment(aes(x = 0, xend=x_max, y=0, yend=x_max),
                                  linetype=2)

        # Set the axis limits
        max_total <- 0
        if (length(real_counts_total) > 0) {
            max_total <- as.integer(max(real_counts_total))
        }

        if (x_max < max_total) {
            message(paste0("WARNING: x-max value ", x_max,
                           " is less than max real data ", max_total))
        }

        max_unique <- 0
        if (length(real_counts_unique) > 0) {
            max_unique <- as.integer(max(real_counts_unique))
            max_unique <- max_unique + (max_unique * 0.1)
        }
        preseq_ymax <- global_y_max_ccurve_limit
        preseq_ymax <- preseq_ymax + (global_y_max_ccurve_limit * 0.1)

        default_ylim <- 100000
        if (coverage > 0) {
            default_ylim <- as.numeric(default_ylim) / coverage
        }
        fig <- fig +
            coord_cartesian(xlim=c(x_min, x_max),
                            ylim = c(default_ylim,
                                     max(preseq_ymax, max_unique)))

        if (preseq_ymax < max_unique) {
            message(paste0("WARNING: y-max value changed from default ",
                           int(preseq_ymax), " to the max real data ",
                           max_unique))
        }

        # label the axis
        # Change labels if we're using coverage
        if (coverage > 0) {
            if (length(real_counts_unique) > 0) {
                fig <- fig +
                    xlab(paste0("Total Coverage (incl. duplicates)\n",
                                "Points show read count versus deduplicated ",
                                "read counts (externally calculated)"))
            } else if (length(real_counts_total) > 0) {
                fig <- fig +
                    xlab(paste0("Total Coverage (incl. duplicates)\n",
                                "Points show externally calculated read ",
                                "counts on the curves"))
            }
            fig <- fig +
                ylab("Unique Coverage") +
                ggtitle("Complexity Curve: preseq")
        } else {
            if (length(real_counts_unique) > 0) {
                fig <- fig +
                    xlab(paste0("Total Molecules (incl. duplicates)\n",
                                "Points show read count versus deduplicated ",
                                "read counts (externally calculated)"))
            } else if (length(real_counts_total) > 0) {
                fig <- fig +
                    xlab(paste0("Total Molecules (incl. duplicates)\n",
                                "Points show externally calculated read ",
                                "counts on the curves"))
            } else {
                fig <- fig +
                    xlab(paste0("Total Molecules (incl. duplicates)"))
            }
            fig <- fig +
                ylab("Unique Molecules") +
                ggtitle("Complexity Curve: preseq")
        }

        fig <- fig +
            labs(col = "") +
            scale_color_discrete(labels=c(sample_name)) +
            theme_classic(base_size=14) +
            theme(axis.line = element_line(size = 0.5)) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black",
                                              fill=NA, size=0.5)) +
            theme(plot.title = element_text(hjust = 0.5))

        # now save the plot
        pdf(file = paste0(tools::file_path_sans_ext(output_name), ".pdf"),
            width= 9, height = 7, useDingbats=F)
        print(fig)
        invisible(dev.off())
        png(filename = paste0(tools::file_path_sans_ext(output_name), ".png"),
            width = 617, height = 480)
        print(fig)
        invisible(dev.off())
    }
}

#' Compute the axis value limit
#'
#' This function returns the index of ccurve_TOTAL_READS containing the
#' closest value to x_max
#' @param value An axis limit value.
#' @param ccurve_TOTAL_READS A vector of read counts from a sample.
#' @keywords preseq limit
#' @examples
#' computeLimit()
computeLimit <- function(value, ccurve_TOTAL_READS) {
    # This function returns the index of ccurve_TOTAL_READS containing the
    # closest value to x_max

    if (max(ccurve_TOTAL_READS) < value) {
        message(paste0("WARNING: ", value, " is set higher than the highest ",
                       "extrapolated point by preseq (value=",
                       max(ccurve_TOTAL_READS)))
    }
    first_point  <- 0
    middle_point <- 0
    last_point   <- length(ccurve_TOTAL_READS)
    iterations   <- 0
    while (first_point != last_point) {
        middle_point <- (first_point + last_point)/2
        middle_value <- ccurve_TOTAL_READS[middle_point]
        if (middle_value == value || iterations >= 10000) {
            return(middle_point)
        } else if (middle_value >= value) {
            last_point  <- middle_point - 1
        } else {
            first_point <- middle_point + 1
        }
        iterations <- iterations + 1
    }
    return(first_point)
}

#' Calculate the Fraction of Reads in Features (FRiF)
#'
#' This function calculates the fraction of reads in a feature and returns
#' a modified BED file with the cumulative sum of reads, cumulative size
#' of covered features, the fraction of reads in those features, and the
#' number of total features.
#'
#' @param bedFile A BED format file
#' @param reads Number of aligned reads
#' @keywords FRiF
#' @examples
#' calcFRiF()
calcFRiF <- function(bedFile, reads) {
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
    bedFile <- cbind(bedFile, frip=bedFile$cumsum/as.numeric(reads))
    bedFile <- cbind(bedFile, numfeats=as.numeric(1:nrow(bedFile)))
    return(bedFile)
}


#' Plot Fraction of Reads in Features (FRiF)
#'
#' This function plots the fraction of reads in a set of features
#' and produces pdf/png output files.
#'
#' @param sample_name Name of sample
#' @param num_reads Number of aligned reads in sample
#' @param output_name Output file name
#' @param bedFile A BED format file
#' @keywords FRiP FRiF BED
#' @export
#' @examples
#' data("promoter")
#' data("promoter_flanking")
#' data("exon")
#' data("intron")
#' data("utr3")
#' data("utr5")
#' plotFRiF(sample_name="example", num_reads=87520,
#'          output_name="example_frif.pdf",
#'          bedFile = c("promoter", "promoter_flanking", "exon",
#'                      "intron", "utr3", "utr5"))
#' @export
plotFRiF <- function(sample_name, num_reads, output_name, bedFile) {
    labels  <- data.frame(xPos=numeric(), yPos=numeric(), name=character(),
                          val=numeric(), color=character(),
                          stringsAsFactors=FALSE)
    palette <- colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))
    plotColors <- palette(length(bedFile))
    if (!exists(bedFile[1])) {
        info   <- file.info(file.path(bedFile[1]))
    }
    name       <- sample_name

    if (exists(bedFile[1])) {
        bed        <- get(bedFile[1])
        bedCov     <- calcFRiF(bed, num_reads)
        name       <- bedFile[1]
        labels[1,] <- c(0.95*max(log10(bedCov$cumSize)), max(bedCov$frip)+0.001,
                        name, round(max(bedCov$frip),2), "#FF0703")
        bedCov$feature <- name
    } else if (file.exists(file.path(bedFile[1])) && info$size != 0) {
        bed        <- read.table(file.path(bedFile[1]))
        bedCov     <- calcFRiF(bed, num_reads)
        name       <- basename(tools::file_path_sans_ext(bedFile[1]))
        name       <- gsub(sample_name, "", name)
        name       <- gsub("^.*?_", "", name)
        numFields  <- 2
        for(i in 1:numFields) name <- gsub("_[^_]*$", "", name)
        labels[1,] <- c(0.95*max(log10(bedCov$cumSize)), max(bedCov$frip)+0.001,
                        name, round(max(bedCov$frip),2), "#FF0703")
        bedCov$feature <- name
    }  else {
        if (is.na(info[1])) {
            message(paste0(name, " coverage file is missing"))
        } else if (info$size == 0) {
            message(paste0(name, " coverage file is empty"))
        } else {
            message(paste0(name, " coverage file is missing"))
        }
    }

    if (exists("bedCov")) {
        covDF <- bedCov
    }

    if (length(bedFile) > 1) {
        for (i in 2:length(bedFile)) {
            if (exists(bedFile[i])) {
                name       <- bedFile[i]
            } else {
                info       <- file.info(file.path(bedFile[i]))
                name       <- basename(tools::file_path_sans_ext(bedFile[i]))
                name       <- gsub(sample_name, "", name)
                name       <- gsub("^.*?_", "", name)
                numFields  <- 2
                for(j in 1:numFields) name <- gsub("_[^_]*$", "", name)
            }

            if (exists(bedFile[i])) {
                bed     <- get(bedFile[i])
            } else if (file.exists(file.path(bedFile[i])) && info$size != 0) {
                bed     <- read.table(file.path(bedFile[i]))
            } else {
                outFile <- file.path(output_name)
                system2(paste("touch"), outFile)
                quit()
            }

            if (max(bed[,4] > 0)) {
                if (exists("covDF")) {
                    covFile <- calcFRiF(bed, num_reads)
                    covFile$feature <- name
                    covDF   <- rbind(covDF, covFile)
                    labels  <- rbind(labels, c(0.95*max(log10(covFile$cumSize)),
                                               max(covFile$frip)+0.001,
                                               name, round(max(covFile$frip),2),
                                               plotColors[i]))
                } else {
                    covDF         <- calcFRiF(bed, num_reads)
                    covDF$feature <- name
                    labels        <- rbind(labels,
                                           c(0.95*max(log10(covDF$cumSize)),
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

    if (!is.null(bedFile)) {
        # Produce plot with bed files
        p <- ggplot(covDF, aes(x=log10(cumSize), y=frip,
                               group=feature, color=feature)) +
            geom_line() +
            labs(x="log(number of bases)", y="FRiF") +
            theme_classic() +
            theme(panel.border = element_rect(colour = "black",
                                              fill=NA,
                                              size=0.5))

        # Recolor and reposition legend
        p <- p + scale_color_manual(labels=paste0(labels$name, ": ",
                                                  labels$val),
                                    values=labels$color) +
            theme(legend.position=c(0.05,0.95),
                  legend.justification=c(0.1,0.9))
    } else {
        write("Unable to produce FRiF plot!\n", stdout())
    }

    if (!exists("p")) {
        p <- ggplot()
    }

    pdf(file = paste0(tools::file_path_sans_ext(output_name), ".pdf"),
        width= 7, height = 7, useDingbats=F)
    print(p)
    invisible(dev.off())
    png(filename = paste0(tools::file_path_sans_ext(output_name), ".png"),
        width = 480, height = 480)
    print(p)
    invisible(dev.off())

    if (exists("p")) {
        write("Cumulative FRiF plot completed!\n", stdout())
    } else {
        write("Unable to produce FRiF plot!\n", stdout())
    }
}


#' Plot TSS enrichment
#'
#' This function plots the global TSS enrichment and produces pdf/png files.
#'
#' @param TSSfile TSS enrichment file
#' @keywords TSS enrichment
#' @export
#' @examples
#' data("tss")
#' plotTSS(TSSfile = "TSSfile")
#' @export
plotTSS <- function(TSSfile) {
    write(paste("\nGenerating TSS plot with ", TSSfile, sep=""), stdout())

    t1 <- theme(
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.line    = element_blank(),
        axis.text.x  = element_text(face = "plain", color = "black",
                                    size = 20, hjust = 0.5),
        axis.text.y  = element_text(face = "plain", color = "black",
                                    size = 20, hjust = 0.5),
        axis.title.x = element_text(face = "plain", color = "black", size = 22,
                                    hjust = 0.5, vjust=0.5),
        axis.title.y = element_text(face = "plain", color = "black", size = 22,
                                    hjust = 0.5),
        plot.title   = element_text(face="bold", color = "black", size=12,
                                    hjust=0.5),
        legend.position="none",
        axis.ticks.length = unit(2, "mm")
    )

    if (exists(TSSfile)) {
        insertionsMat <- data.frame(get(TSSfile))
    } else {
        insertionsMat <- read.table(TSSfile, header=FALSE, row.names=NULL,
                                    as.is=TRUE, check.names=FALSE)
    }

    normTSS <- insertionsMat / mean(insertionsMat[1:200,])
    colnames(normTSS) <- c("score")
    TSSscore <- round(mean(normTSS[1950:2050,]),1)
    if (is.nan(TSSscore)) {
        message(paste("\nNaN produced.  Check ", TSSfile, "\n", sep=""))
        quit()
    }
    lineColor <- "red2"
    if (TSSscore > TSS_CUTOFF)
    {
        lineColor <- "springgreen4"
    }


    png(filename = paste(tools::file_path_sans_ext(TSSfile), ".png", sep=""),
        width = 480, height = 480)
    pre <- ggplot(normTSS, aes(x=(as.numeric(rownames(normTSS))-2000), y=score,
                               group=1, colour="black")) +
        geom_hline(yintercept = 6, linetype = 2,
                   color = "grey", size = 0.25) +
        geom_smooth(method="loess", span=0.02,
                    se=FALSE, colour=lineColor) +
        labs(x = "Distance from TSS (bp)", y = "TSS Enrichment Score")
    p <- pre + t1 + scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        coord_cartesian(xlim=c(-2300,2300), ylim=c(0,32)) +
        #ylim(0,30) + xlim(-2100,2100) +
        annotate("rect", xmin=1200, xmax=2300, ymin=27, ymax=32,
                 fill="gray95", size = 0.5) +
        annotate("text", x=1750, y=31, label="TSS Score", fontface = 1,
                 size=6, hjust=0.5) +
        annotate("text", x=1750, y=29, label=TSSscore, fontface = 2,
                 size=10, hjust=0.5)
    print(p)
    invisible(dev.off())

    pdf(file = paste(tools::file_path_sans_ext(TSSfile), ".pdf", sep=""),
        width= 7, height = 7, useDingbats=F)
    pre <- ggplot(normTSS, aes(x=(as.numeric(rownames(normTSS))-2000), y=score,
                               group=1, colour="black")) +
        geom_hline(yintercept = 6, linetype = 2,
                   color = "grey", size = 0.25) +
        geom_smooth(method="loess", span=0.02,
                    se=FALSE, colour=lineColor) +
        labs(x = "Distance from TSS (bp)", y = "TSS Enrichment Score")
    p <- pre + t1 + scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        coord_cartesian(xlim=c(-2300,2300), ylim=c(0,32)) +
        #ylim(0,30) + xlim(-2100,2100) +
        annotate("rect", xmin=1200, xmax=2300, ymin=27, ymax=32,
                 fill="gray95", size = 0.5) +
        annotate("text", x=1750, y=31, label="TSS Score", fontface = 1,
                 size=6, hjust=0.5) +
        annotate("text", x=1750, y=29, label=TSSscore, fontface = 2,
                 size=10, hjust=0.5)
    print(p)
    invisible(dev.off())

    write("Completed TSS enrichment plot!\n", stdout())
}

################################################################################
