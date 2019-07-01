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
#' @param ignore_unique If FALSE, ignore information about unique read counts
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
                                   real_counts_path=NA, ignore_unique=FALSE,
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
    rcDT <- data.table(name = character(length(real_counts_path)),
                       total = integer(length(real_counts_path)),
                       unique = integer(length(real_counts_path)))
    if ("data.frame" %in% class(real_counts_path[1])) {
        rc_file    <- real_counts_path[1]
        rcDT$name  <- basename(rc_file$V1)
        rcDT$total <- as.integer(rc_file$V2)
        if (ncol(rc_file) == 3 && !ignore_unique) {
            rcDT$unique <- as.integer(rc_file$V3)
        } else {
            rcDT$unique <- NA
        }

        if (length(real_counts_path) > 1) {
            for (rc in 2:length(real_counts_path)) {
                rc_file <- real_counts_path[rc]
                rcDT$name[rc]  <- basename(rc_file$V1)
                rcDT$total[rc] <- as.integer(rc_file$V2)
                if (ncol(rc_file) == 3 && !ignore_unique) {
                    rcDT$unique[rc] <- as.integer(rc_file$V3)
                } else {
                    rcDT$unique[rc] <- NA
                }
            }
        }
    } else if (!is.na(real_counts_path[1])) {
        for (rc in 1:length(real_counts_path)) {
            info <- file.info(file.path(real_counts_path[rc]))
            if (file.exists(real_counts_path[rc]) && info$size != 0) {
                rc_file        <- fread(real_counts_path[rc])
                rcDT$name[rc]  <- basename(rc_file$V1)
                rcDT$total[rc] <- as.integer(rc_file$V2)
                if (ncol(rc_file) == 3 && !ignore_unique) {
                    rcDT$unique[rc] <- as.integer(rc_file$V3)
                } else {
                    rcDT$unique[rc] <- NA
                }
            } else {
                message(paste0("Error loading real counts file: ",
                               real_counts_path[rc]))
                if (!file.exists(real_counts_path[rc])) {
                    message("File could not be found.")
                } else if (info$size == 0) {
                    message("File is empty.")
                }
                quit()
            }
        }
    }

    # Convert real counts to coverage
    if (coverage > 0) {
        rcDT[,total  := as.numeric(total)  / coverage]
        rcDT[,unique := as.numeric(unique) / coverage]
    }

    # Set up plot params
    global_x_max_ccurve_limit <- 0
    global_y_max_ccurve_limit <- 0
    max_label_length          <- 0

    # Each ccurve will get a different color
    palette <- colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))

    clist <- data.table(TOTAL_READS = list(), EXPECTED_DISTINCT = list())
    ccurves  <- as.list(ccurves)
    colormap <- palette(length(ccurves))
    for (c in 1:length(ccurves)) {
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
                SAMPLE_NAME = list(rep(sample_name,
                                   length(ccurve_TOTAL_READS))),
                TOTAL_READS = list(ccurve_TOTAL_READS),
                EXPECTED_DISTINCT = list(ccurve_EXPECTED_DISTINCT),
                COLOR = list(rep(colormap[c], length(ccurve_TOTAL_READS)))
            )
        } else {
            clist <- rbindlist(list(clist,
                data.table(SAMPLE_NAME = list(rep(sample_name,
                                              length(ccurve_TOTAL_READS))),
                           TOTAL_READS = list(ccurve_TOTAL_READS),
                           EXPECTED_DISTINCT = list(ccurve_EXPECTED_DISTINCT),
                           COLOR = list(rep(colormap[c],
                                        length(ccurve_TOTAL_READS))))),
                use.names=TRUE)
        }
    }

    x_min_ccurve_limit <- computeLimit(x_min, unlist(clist$TOTAL_READS))
    x_max_ccurve_limit <- computeLimit(x_max, unlist(clist$TOTAL_READS))
    if (x_max_ccurve_limit > global_x_max_ccurve_limit) {
        global_x_max_ccurve_limit <- x_max_ccurve_limit
    }
    if (unlist(clist$EXPECTED_DISTINCT)[x_max_ccurve_limit] > global_y_max_ccurve_limit) {
        if (x_max_ccurve_limit <= length(unlist(clist$EXPECTED_DISTINCT))) {
            global_y_max_ccurve_limit <- unlist(clist$EXPECTED_DISTINCT)[x_max_ccurve_limit]
        } else {
            x_max_ccurve_limit <- length(unlist(clist$EXPECTED_DISTINCT))
        }
    }
    # Add a few points to be sure
    x_max_ccurve_limit <- x_max_ccurve_limit + 3

    sn <- clist$SAMPLE_NAME[[1]][x_min_ccurve_limit:x_max_ccurve_limit]
    tr <- clist$TOTAL_READS[[1]][x_min_ccurve_limit:x_max_ccurve_limit]
    ed <- clist$EXPECTED_DISTINCT[[1]][x_min_ccurve_limit:x_max_ccurve_limit]
    co <- clist$COLOR[[1]][x_min_ccurve_limit:x_max_ccurve_limit]
    df <- data.frame(sample_name = sn,
                     total_reads = tr,
                     expected_distinct = ed,
                     color = co)
    if (nrow(clist) > 1) {
        for (i in 2:nrow(clist)) {
            sn <- clist$SAMPLE_NAME[[i]][x_min_ccurve_limit:x_max_ccurve_limit]
            tr <- clist$TOTAL_READS[[i]][x_min_ccurve_limit:x_max_ccurve_limit]
            ed <- clist$EXPECTED_DISTINCT[[i]][x_min_ccurve_limit:x_max_ccurve_limit]
            co <- clist$COLOR[[i]][x_min_ccurve_limit:x_max_ccurve_limit]
            df <- rbind(df, data.frame(sample_name = sn,
                                       total_reads = tr,
                                       expected_distinct = ed,
                                       color = co))
        }
    }

    # Plot the curve
    fig <- ggplot(df, aes(total_reads,
                          expected_distinct,
                          group = sample_name,
                          col=color)) + geom_line()

    # Plot the real data if we have it
    numFields   <- 2
    for(j in 1:numFields) rcDT$name <- gsub("_[^_]*$", "", rcDT$name)
    rcDT$color <- colormap

    if (any(rcDT$total > 0) && !any(is.na(rcDT$unique)) && !ignore_unique) {
        fig <- fig + geom_point(data=rcDT,
                                aes(total, unique, col=color),
                                shape=23, size=3)
        message(paste0("INFO: Found real counts for ",
                       paste(rcDT$name, sep=","), " - Total: ",
                       rcDT$total, " Unique: ",
                       rcDT$unique, "\n"))
    } else if (any(rcDT$total > 0)) {
        if (max(rcDT$total) > max(df$total_reads)) {
            message(paste0("WARNING: Max total reads (", max(rcDT$total),
                           ") > ", "max preseq value (", max(df$total_reads),
                           ") - skipping..."))
        } else {
            interp <- approx(df$total_reads, df$expected_distinct, rcDT$total)$y
            fig <- fig + geom_point(data=rcDT,
                                    aes(total, interp, col=color),
                                    shape=23, size=3)
            message(paste0("INFO: Found real counts for ",
                           paste(rcDT$name, sep=","), " - Total: ",
                           rcDT$total, " (preseq unique reads: ",
                           interp, ")\n"))
        }
    } else {
        message(paste0("INFO: No real counts provided."))
    }

    # plot perfect library as dashed line
    fig <- fig + geom_segment(aes(x = 0, xend=x_max, y=0, yend=x_max),
                              linetype=2, col ='black')

    # Set the axis limits
    max_total <- 0
    if (any(rcDT$total > 0)) {
        max_total <- as.numeric(max(rcDT$total))
    }

    if (x_max < max_total) {
        message(paste0("WARNING: x-max value ", x_max,
                       " is less than max real data ", max_total))
    }

    max_unique <- 0
    if (!any(is.na(rcDT$unique)) && any(rcDT$unique > 0)) {
        max_unique <- as.numeric(max(rcDT$unique))
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
                       preseq_ymax, " to the max real data ",
                       max_unique))
    }

    # label the axis
    # Change labels if we're using coverage
    if (coverage > 0) {
        if (!any(is.na(rcDT$unique)) && any(rcDT$unique > 0)) {
            fig <- fig +
                xlab(paste0("Total Coverage (incl. duplicates)\n",
                            "Points show read count versus deduplicated ",
                            "read counts (externally calculated)"))
        } else if (any(rcDT$total > 0)) {
            fig <- fig +
                xlab(paste0("Total Coverage (incl. duplicates)\n",
                            "Points show read count versus projected unique ",
                            "read counts on the curves"))
        } else {
            fig <- fig +
                xlab(paste0("Total Coverage (incl. duplicates)"))
        }
        fig <- fig +
            ylab("Unique Coverage") +
            ggtitle("Complexity Curve: preseq")
    } else {
        if (!any(is.na(rcDT$unique)) && any(rcDT$unique > 0)) {
            fig <- fig +
                xlab(paste0("Total Reads (incl. duplicates)\n",
                            "Points show read count versus deduplicated ",
                            "read counts (externally calculated)"))
        } else if (any(rcDT$total > 0)) {
            fig <- fig +
                xlab(paste0("Total Reads (incl. duplicates)\n",
                            "Points show externally calculated read ",
                            "counts on the curves"))
        } else {
            fig <- fig +
                xlab(paste0("Total Reads (incl. duplicates)"))
        }
        fig <- fig +
            ylab("Unique Reads") +
            ggtitle("Complexity Curve: preseq")
    }

    fig <- fig +
        labs(col = "") +
        scale_color_discrete(labels=c(clist$SAMPLE_NAME)) +
        theme_classic(base_size=14) +
        theme(axis.line = element_line(size = 0.5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio = 1,
              panel.border = element_rect(colour = "black",
                                          fill=NA, size=0.5)) +
        theme(plot.title = element_text(hjust = 0.5))

    # inset zoom plot
    zoom_theme <- theme(legend.position = "none",
                        axis.line = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        aspect.ratio = 1,
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(color='black'),
                        plot.margin = unit(c(0,0,-6,-6),"mm"))

    if (!any(is.na(rcDT$unique)) && any(rcDT$unique > 0)) {
        zoom_fig <- ggplot(df, aes(total_reads,
                       expected_distinct,
                       group = sample_name,
                       col=color)) +
            geom_line() +
            geom_abline(intercept = 0, slope = 1, linetype="dashed") +
            coord_cartesian(xlim = c(0,max(rcDT$unique)*2),
                            ylim = c(0,max(rcDT$unique)*2)) +
            #geom_hline(data=rcDT, aes(yintercept=unique, col=color),
            #           linetype="dotted") +
            #geom_vline(data=rcDT, aes(xintercept=unique, col=color),
            #           linetype="dotted") +
            geom_point(data=rcDT,
                       aes(total, unique, col=color),
                       shape=23, size=3) +
            theme_classic(base_size=14) +
            zoom_theme
        g   <- ggplotGrob(zoom_fig)
        fig <- fig +
            annotation_custom(grob = g,
                              xmin = x_max / 2,
                              xmax = x_max,
                              ymin = 0,
                              ymax = max(preseq_ymax, max_unique)/2)
    } else if (any(rcDT$total > 0)) {
        interp <- approx(df$total_reads, df$expected_distinct, rcDT$total)$y
        zoom_fig <- ggplot(df, aes(total_reads,
                                   expected_distinct,
                                   group = sample_name,
                                   col=color)) +
            geom_line() +
            geom_abline(intercept = 0, slope = 1, linetype="dashed") +
            coord_cartesian(xlim = c(0,max(rcDT$total)*2),
                            ylim = c(0,max(rcDT$total)*2)) +
            # geom_hline(data=rcDT, aes(yintercept=total, col=color),
            #            linetype="dotted") +
            # geom_vline(data=rcDT, aes(xintercept=total, col=color),
            #            linetype="dotted") +
            geom_point(data=rcDT,
                       aes(total, interp, col=color),
                       shape=23, size=3) +
            theme_classic(base_size=14) +
            zoom_theme
        g   <- ggplotGrob(zoom_fig)
        fig <- fig +
            annotation_custom(grob = g,
                              xmin = x_max / 2,
                              xmax = x_max,
                              ymin = 0,
                              ymax = max(preseq_ymax, max_unique)/2)
    }

    # now save the plot
    pdf(file = paste0(tools::file_path_sans_ext(output_name), ".pdf"),
        width= 10, height = 7, useDingbats=F)
    print(fig)
    invisible(dev.off())
    png(filename = paste0(tools::file_path_sans_ext(output_name), ".png"),
        width = 686, height = 480)
    print(fig)
    invisible(dev.off())
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
    if (length(TSSfile) == 1) {
        write(paste0("\nGenerating TSS plot with ", TSSfile), stdout())
    } else {
        if (length(TSSfile) == 2) {
            write(paste0("\nGenerating TSS plot with ",
                         paste(TSSfile, collapse=" and ")),
                  stdout())
        } else {
            write(paste0("\nNot sure how to merge the following: ",
                         paste(TSSfile, collapse=", ")),
                  stdout())
            write(paste0("Did you mean to pass more than 2 files?"), stdout())
            quit()
        }
    }

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

    iMat <- data.table(V1 = numeric())
    if (length(TSSfile) == 1) {
        if (exists(TSSfile[i])) {
            iMat <- data.table(get(TSSfile))
        } else {
            iMat <- fread(TSSfile)
        }
    } else if (length(TSSfile) == 2) {
        for (i in 1:length(TSSfile)) {
            if (exists(TSSfile[i])) {
                if (i == 1) {
                    iMat <- data.table(get(TSSfile[i]))
                } else {
                    iMat <- list(iMat, data.table(get(TSSfile[i])))
                }
            } else {
                if (i == 1) {
                    iMat <- fread(TSSfile[i])
                } else {
                    iMat <- list(iMat, fread(TSSfile[i]))
                }
            }
        }
    } else {
        write(paste0("\nNot sure how to merge the following: ",
                     paste(TSSfile, collapse=", ")),
              stdout())
        write(paste0("Did you mean to pass more than 2 files?"), stdout())
        quit()
    }

    if (length(TSSfile) == 1) {
        plusMinus <- iMat
    } else {
        plus      <- iMat[[1]]
        minus     <- iMat[[2]]
    }

    if (exists("plusMinus")) {
        val      <- 0.05*nrow(plusMinus)
        #normTSS  <- (plusMinus / mean(plusMinus[c(1:val,
        #            (nrow(plusMinus)-val):nrow(plusMinus)), V1]))
        normTSS           <- plusMinus / mean(plusMinus[c(1:val), V1])
        colnames(normTSS) <- c("score")
        peakPos  <- which.max(normTSS$score)
        TSSscore <- round(mean(normTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                                       peakPos+50)), score]),1)
        if (is.nan(TSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile, "\n"))
            quit()
        }
    } else {
        val      <- 0.05*nrow(plus)
        #normTSS  <- (plus / mean(plus[c(1:val,
        #            (nrow(plus)-val):nrow(plus)), V1]))
        normTSS           <- plus / mean(plus[c(1:val), V1])
        colnames(normTSS) <- c("score")
        peakPos  <- which.max(normTSS$score)
        TSSscore <- round(mean(normTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                                       peakPos+50)), score]),1)
        if (is.nan(TSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile[1], "\n"))
            quit()
        }
    }
    
    lineColor <- "red2"
    if (TSSscore > TSS_CUTOFF)
    {
        lineColor <- "springgreen4"
    }

    name        <- basename(tools::file_path_sans_ext(TSSfile[1]))
    numFields   <- 2
    for(j in 1:numFields) name <- gsub("_[^_]*$", "", name)
    sample_name <- paste(dirname(TSSfile[1]), name, sep="/")

    pre <- ggplot(normTSS, aes(x=(as.numeric(rownames(normTSS))-
                                 (nrow(normTSS)/2)),
                               y=score, group=1, colour="black")) +
        geom_hline(yintercept = 6, linetype = 2,
                   color = "grey", size = 0.25) +
        geom_smooth(method="loess", span=0.02,
                    se=FALSE, colour=lineColor) +
        labs(x = "Distance from TSS (bp)", y = "TSS Enrichment Score")
    p <- pre + t1 +
         scale_x_continuous(expand=c(0,0)) +
         scale_y_continuous(expand=c(0,0)) +
         coord_cartesian(xlim=c(-2300,2300), ylim=c(0,32))
    if (exists("minus")) {
        val      <- 0.025*nrow(minus)
        # normTSS  <- (minus / mean(minus[c(1:val,
        #             (nrow(minus)-val):nrow(minus)), V1]))
        minusNormTSS           <- minus / mean(minus[c(1:val), V1])
        colnames(minusNormTSS) <- c("score")
        peakPos       <- which.max(minusNormTSS$score)
        minusTSSscore <- round(
            mean(minusNormTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                               peakPos+50)), score]),1)
        if (is.nan(minusTSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile[2], "\n"))
            quit()
        }
        p <- p + geom_smooth(data=minusNormTSS,
                             aes(x=(as.numeric(rownames(minusNormTSS))-
                                   (nrow(normTSS)/2)),
                                 y=score, group=1, colour="black"),
                             method="loess", span=0.02,
                             se=FALSE, colour="blue") +
                 annotate("rect", xmin=1200, xmax=2300, ymin=25, ymax=32,
                          fill="gray95", size = 0.5) +
                 annotate("text", x=1750, y=31, label="TSS Score", fontface = 1,
                          size=6, hjust=0.5) +
                 annotate("text", x=1500, y=29, label="+", fontface = 2,
                          size=8, hjust=0.5, color=lineColor) +
                 annotate("text", x=1500, y=27, label=TSSscore, fontface = 2,
                          size=8, hjust=0.5, color=lineColor) +
                 annotate("text", x=2000, y=29, label="-",
                          fontface = 2, size=8, hjust=0.5, color="blue") +
                 annotate("text", x=2000, y=27, label=minusTSSscore,
                          fontface = 2, size=8, hjust=0.5, color="blue")
    } else {
        p <- p + annotate("rect", xmin=1200, xmax=2300, ymin=27, ymax=32,
                          fill="gray95", size = 0.5) +
                 annotate("text", x=1750, y=31, label="TSS Score",
                          fontface = 1, size=6, hjust=0.5) +
                 annotate("text", x=1750, y=29, label=TSSscore, fontface = 2,
                          size=10, hjust=0.5)
    }

    png(filename = paste0(sample_name, "_TSSenrichment.png"),
        width = 480, height = 480)
    print(p)
    invisible(dev.off())

    pdf(file = paste0(sample_name, "_TSSenrichment.pdf"),
        width= 7, height = 7, useDingbats=F)
    print(p)
    invisible(dev.off())

    write("Completed TSS enrichment plot!\n", stdout())
}


#' Plot fragment length distribution
#'
#' This function plots the fragment length distribution of a paired-end sample
#' and produces pdf/png files.
#'
#' @param fragL infile containing single column of fragment lengths
#' @param fragL_count counts of each fragment length identified
#' @param fragL_dis1 pdf filename
#' @param ragL_dis2 fragment length distribution stats file
#' @keywords fragment distribution
#' @export
#' @examples
#' data("fragL")
#' data("fragL_count")
#' plotFLD(fragL = "fragL", fragL_count = "fragL_count",
#'         fragL_dis1 = "fragLenDistribution_example.pdf",
#'         fragL_dis2 = "fragLenDistribution_example.txt")
#' @export
plotFLD <- function(fragL, fragL_count,
                    fragL_dis1="fragLenDistribution.pdf",
                    fragL_dis2="fragLenDistribution.txt") {

    outfile_png <- gsub('pdf', 'png', fragL_dis1)

    dat  <- fread(fragL_count)
    dat1 <- dat[dat$V2<=600,]
    tmp  <- seq(1:as.numeric(dat1[1,2]-1))
    dat0 <- data.table(V1=rep(0,length(tmp)),V2=tmp)
    dat2 <- rbind(dat0, dat1)

    t1 = theme_classic(base_size=14) +
         theme(axis.line = element_line(size = 0.5)) +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.position = "none",
               aspect.ratio = 1,
               panel.border = element_rect(colour = "black",
                                          fill=NA, size=0.5)) +
         theme(plot.title = element_text(hjust = 0.5))

    p <- ggplot(dat1, aes(x=V2, y=V1)) +
             geom_line(aes(color='red')) +
             xlab("Read length") + 
             ylab("Read counts") +
             ggtitle("Insert size distribution") +
             t1

    # Save plot to pdf file
    pdf(file=fragL_dis1, width= 7, height = 7, useDingbats=F)
    print(p)
    invisible(dev.off())
         
    # Save plot to png file
    png(filename=outfile_png, width = 480, height = 480)
    print(p)
    invisible(dev.off())

    dat  <- fread(fragL)
    summ <- data.table(Min=min(dat$V1), Max=max(dat$V1), Median=median(dat$V1),
                       Mean=mean(dat$V1), Stdev=sd(dat$V1))
    # Write summary table to stats file
    fwrite(summ, file=fragL_dis2, row.names=F, quote=F, sep="\t")
}

#' Return the count of the plot at 95% of the upper limit.
#'
#' From:
#' DeCicco, L. (2018, August 10). Exploring ggplot2 boxplots - 
#'  Defining limits and adjusting style [Blog post]. 
#'  Retrieved from https://waterdata.usgs.gov/blog/boxplots/
#'
#' @param x A vector of numbers.
#' @param lim Upper limit.
#' @keywords limits
#' @examples
#' n_fun()
n_fun <- function(x, lim=50) {
    return(data.frame(y = 0.95*log10(lim),
                      label = length(x)))
}

#' This function forces the y-axis breaks to be on every 10^x.
#'
#' From:
#' DeCicco, L. (2018, August 10). Exploring ggplot2 boxplots - 
#'  Defining limits and adjusting style [Blog post]. 
#'  Retrieved from https://waterdata.usgs.gov/blog/boxplots/
#'
#' @param x A vector of numbers.
#' @keywords limits
#' @examples
#' prettyLogs()
prettyLogs <- function(x){
    pretty_range    <- range(x[x > 0])
    pretty_logs     <- 10^(-10:10)
    log_index       <- which(pretty_logs < pretty_range[2] & 
                             pretty_logs > pretty_range[1])
    log_index       <- c(log_index[1]-1,log_index,
                         log_index[length(log_index)]+1)
    pretty_logs_new <-  pretty_logs[log_index] 
    return(pretty_logs_new)
}

#' Custom formatting function for the log axis.
#'
#' From:
#' DeCicco, L. (2018, August 10). Exploring ggplot2 boxplots - 
#'  Defining limits and adjusting style [Blog post]. 
#'  Retrieved from https://waterdata.usgs.gov/blog/boxplots/
#'
#' @param n A vector of numbers.
#' @keywords limits
#' @examples
#' fancyNumbers()
fancyNumbers <- function(n){
    nNoNA     <- n[!is.na(n)]
    x         <- gsub(pattern = "1e",replacement = "10^",
                      x = format(nNoNA, scientific = TRUE))
    exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
    base      <- ifelse(exponents == 0, "1",
                 ifelse(exponents == 1, "10","10^"))
    exponents[base == "1" | base == "10"] <- ""
     textNums            <- rep(NA, length(n))  
     textNums[!is.na(n)] <- paste0(base,exponents)
    textReturn          <- parse(text=textNums)
     return(textReturn)
}

#' Plot the distribution of genic exonRPKM/intronRPKM ratios
#'
#' This function plots the distribution of by gene exon RPKM divided by
#' intron RPKM ratios. Can produce raw or log10 distributions, but reports
#' both median values.
#'
#' @param rpkm A three column TSV format file containing
#'             "gene", "intron RPKM", "exon RPKM" columns.
#' @param raw Plot raw distribution
#' @keywords mRNA contamination
#' @export
#' @examples
#' data("rpkm")
#' mRNAcontamination(rpkm = "rpkm")
#' @export
mRNAcontamination <- function(rpkm, raw=FALSE) {
    if (exists(rpkm)) {
        RPKM <- data.table(get(rpkm))
    } else {
        RPKM <- fread(rpkm)
    }
    colnames(RPKM) <- c("gene","intron","exon")

    name           <- basename(tools::file_path_sans_ext(rpkm))
    numFields      <- 2
    for(j in 1:numFields) name <- gsub("_[^_]*$", "", name)
    sample_name <- paste(dirname(rpkm), name, sep="/")

    finite_rpkm <- RPKM[is.finite(RPKM$exon/RPKM$intron),]

    if (raw) {
        q <- ggplot(data = finite_rpkm, 
                    aes(x="", y=(exon/intron))) +
                stat_boxplot(geom ='errorbar', width = 0.25) +
                geom_boxplot(width = 0.25,
                             outlier.color='red',
                             outlier.shape=1) +
                stat_summary(fun.y = "mean", geom = "point",
                             shape = 1, size = 2) +
                labs(x=name,
                     y=expression((over(exon[RPKM], intron[RPKM]))~X~Gene)) +
                ylim(c(0, ceiling(summary(finite_rpkm$exon/finite_rpkm$intron)[5]))) +
                theme_classic(base_size=14) +
                theme(axis.line = element_line(size = 0.5)) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      aspect.ratio = 1,
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, size=0.5))
    } else {
        q <- ggplot(data = finite_rpkm, 
                    aes(x="", y=log10(exon/intron))) +
                stat_boxplot(geom ='errorbar', width = 0.25) +
                geom_boxplot(width = 0.25,
                             outlier.color='red',
                             outlier.shape=1) +
                stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
                stat_summary(fun.y = "mean", geom = "point",
                             shape = 1, size = 2) +
                scale_y_log10(limits = c(0.001, 50),
                              expand = expand_scale(mult = c(0, 0)),
                              labels=fancyNumbers,
                              breaks=prettyLogs) +
                annotation_logticks(sides = c("rl")) +
                labs(x=name,
                     y=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene)) +
                theme_classic(base_size=14) +
                theme(axis.line = element_line(size = 0.5)) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      aspect.ratio = 1,
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, size=0.5))
    }

    label1 <- c(paste("'median'[log[10]]", ":~",
                round(median(log10((finite_rpkm$exon/finite_rpkm$intron))), 2)),
                paste("'median'[raw]", ":",
                round(median(finite_rpkm$exon/finite_rpkm$intron), 2)))

    max_y  <- layer_scales(q)$y$range$range[2]

    if (raw) {
        q <- q + annotate("text", x = 0.5, y = c(max_y, 0.95*max_y),
                      hjust=0, vjust=1, label = label1, parse=TRUE)
    } else {
        q <- q + annotate("text", x = 0.5, y = c(10^max_y, 10^max_y-10),
                      hjust=0, vjust=1, label = label1, parse=TRUE)
    }

    # Save plot to pdf file
    pdf(file=paste0(sample_name, "_mRNA_contamination.pdf"),
        width= 7, height = 7, useDingbats=F)
    print(q)
    invisible(dev.off())
         
    # Save plot to png file
    png(filename = paste0(sample_name, "_mRNA_contamination.png"),
        width = 480, height = 480)
    print(q)
    invisible(dev.off())
}

################################################################################
