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
#' A standardized ggplot theme for PEPPRO plots
#'
#' @keywords ggplot2 theme
#' @examples
#' theme_PEPPRO()
theme_PEPPRO <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 14, ...) +
  theme(
    axis.line = element_line(size = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    aspect.ratio = 1,
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )
}


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
#' @param x_min Lower x-limit (default 0)
#' @param x_max Upper x-limit (default 500 million)
#' @keywords preseq library complexity
#' @export
#' @examples
#' data("ccurve")
#' data("counts")
#' plotComplexityCurves(ccurves = "ccurve", coverage=3099922541, read_length=30,
#'                      real_counts_path="counts")
plotComplexityCurves <- function(ccurves,
                                 coverage=0, read_length=0,
                                 real_counts_path=NA, ignore_unique=FALSE,
                                 x_min=0, x_max=500000000) {

    if (x_min < 0 || x_max <= x_min) {
        message(paste0("problem with x-min or x-max (", x_min, " ", x_max,
                       "). x-min must be >= 0 and < x-max"))
        quit(save = "no", status = 1, runLast = FALSE)
    }

    # Convert limit counts to coverage
    if (coverage > 0) {
        if (read_length == 0) {
            message("Error: --coverage specified but not --read_length")
            quit(save = "no", status = 1, runLast = FALSE)
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

    if (exists(real_counts_path[1])) {
        rc_file    <- data.table(get(real_counts_path[1]))
        rcDT$name  <- basename(rc_file$V1)
        rcDT$total <- as.integer(rc_file$V2)
        if (ncol(rc_file) == 3 && !ignore_unique) {
            rcDT$unique <- as.integer(rc_file$V2)
        } else {
            rcDT$unique <- NA
        }

        if (length(real_counts_path) > 1) {
            for (rc in 2:length(real_counts_path)) {
                rc_file <- real_counts_path[rc]
                rcDT$name[rc]  <- basename(rc_file$V1)
                rcDT$total[rc] <- as.integer(rc_file$V3)
                if (ncol(rc_file) == 3 && !ignore_unique) {
                    rcDT$unique[rc] <- as.integer(rc_file$V2)
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
                rcDT$total[rc] <- as.integer(rc_file$V3)
                if (ncol(rc_file) == 3 && !ignore_unique) {
                    rcDT$unique[rc] <- as.integer(rc_file$V2)
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
                quit(save = "no", status = 1, runLast = FALSE)
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

        if (exists(ccurves[[c]])) {
            ctable <- data.table(get(ccurves[[c]]))
        } else if (file.exists(ccurves[[c]])) {
            ctable <- fread(ccurves[[c]])
        } else {
            stop(paste0("FileExistsError: ", ccurves[[c]],
                        " could not be found."))
            quit(save = "no", status = 1, runLast = FALSE)
        }

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
                quit(save = "no", status = 1, runLast = FALSE)
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
                quit(save = "no", status = 1, runLast = FALSE)
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

    # Plot by millions of reads
    plottingFactor <- 1000000

    df$total_reads       <- df$total_reads/plottingFactor
    df$expected_distinct <- df$expected_distinct/plottingFactor
    rcDT$total  <- rcDT$total/plottingFactor
    rcDT$unique <- rcDT$unique/plottingFactor
    x_max       <- x_max/plottingFactor

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
                       paste(rcDT$name, sep=","), " - Total (M): ",
                       rcDT$total, " Unique (M): ",
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
                           paste(rcDT$name, sep=","), " - Total (M): ",
                           rcDT$total, " (preseq unique reads (M): ",
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

    # Adjust limits by plottingFactor
    default_ylim <- default_ylim/plottingFactor
    preseq_ymax  <- preseq_ymax/plottingFactor

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
                labs(x = paste0("total coverage (incl. duplicates)"),
                     caption = paste0("Points show read count versus ",
                                      "deduplicated read counts ",
                                      "(externally calculated)"))
        } else if (any(rcDT$total > 0)) {
            fig <- fig +
                labs(x = "total coverage (incl. duplicates)",
                     caption = paste0("Points show read count versus projected ",
                                      "unique read counts on the curves"))
        } else {
            fig <- fig +
                labs(x = "total coverage (incl. duplicates)")
        }
        fig <- fig +
            labs = (y = "unique coverage")
            #ggtitle("Complexity Curve: preseq")
    } else {
        if (!any(is.na(rcDT$unique)) && any(rcDT$unique > 0)) {
            fig <- fig +
                labs(x = "total reads (M) (incl. duplicates)",
                     caption = paste0("Points show read count versus deduplicated ",
                                      "read counts (externally calculated)"))
        } else if (any(rcDT$total > 0)) {
            fig <- fig +
                labs(x = "total reads (M) (incl. duplicates)",
                     caption = paste0("Points show externally calculated read ",
                                      "counts on the curves"))
        } else {
            fig <- fig +
                labs(x = "total reads (M) (incl. duplicates)")
        }
        fig <- fig +
            labs(y = "unique reads (M)")
            #ggtitle("Complexity Curve: preseq")
    }

    fig <- fig +
        labs(col = "") +
        scale_color_discrete(labels=c(clist$SAMPLE_NAME)) +
        theme_PEPPRO() +
        theme(legend.position = "right",
              plot.caption = element_text(size = 8, face = "italic"))

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
                        plot.margin = unit(c(0.1,0.1,-6,-6),"mm"))

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
                              ymin = 10,
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
                              ymin = 10,
                              ymax = max(preseq_ymax, max_unique)/2)
    }

    return(fig)
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
        middle_point <- as.numeric((first_point + last_point)/2)
        middle_value <- as.numeric(ccurve_TOTAL_READS[middle_point])
        if (length(middle_value)==0) {
            return(middle_point)
        } else if (middle_value == value || iterations >= 10000) {
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
    colnames(bedFile) <- c("chromosome", "start", "end", "count")
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
#' @param genome_size Size of genome in bp
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
plotFRiF <- function(sample_name, num_reads, genome_size,
                     type = c("frif", "prif", "both"),
                     output_name, bedFile) {
    # TODO: Get total genome size value as an input
    #genome_size <- 3099922541
    labels  <- data.frame(xPos=numeric(), yPos=numeric(), name=character(),
                          val=numeric(), color=character(),
                          stringsAsFactors=FALSE)
    feature_dist  <- data.frame(feature=character(), numfeats=numeric(),
                                numbases=numeric(), expected=numeric(),
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
        feature_dist[1,] <- c(name, nrow(bed),
                              as.numeric(sum(abs(bed$V3-bed$V2))),
                              as.numeric((sum(abs(bed$V3-bed$V2))/genome_size)))
        bedCov$feature <- name
    } else if (file.exists(file.path(bedFile[1])) && info$size != 0) {
        bed <- read.table(file.path(bedFile[1]))
        if (nrow(bed[which(bed$V4 != 0),]) == 0) {
            message(paste0(name, "  has no covered features"))
        } else {
            bedCov     <- calcFRiF(bed, num_reads)
            name       <- basename(tools::file_path_sans_ext(bedFile[1]))
            name       <- gsub(sample_name, "", name)
            name       <- gsub("^.*?_", "", name)
            numFields  <- 2
            for(i in 1:numFields) name <- gsub("_[^_]*$", "", name)
            labels[1,] <- c(0.95*max(log10(bedCov$cumSize)), max(bedCov$frip)+0.001,
                            name, round(max(bedCov$frip),2), "#FF0703")
            feature_dist[1,] <- c(name, nrow(bed),
                                  as.numeric(sum(abs(bed$V3-bed$V2))),
                                  as.numeric((sum(abs(bed$V3-bed$V2))/genome_size)))
            bedCov$feature <- name
        }
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
    } else {
        return(ggplot())
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
                quit(save = "no", status = 1, runLast = FALSE)
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
                    feature_dist <- rbind(feature_dist,
                        c(name, nrow(bed), as.numeric(sum(abs(bed$V3-bed$V2))),
                          as.numeric((sum(abs(bed$V3-bed$V2))/genome_size))))
                } else {
                    covDF         <- calcFRiF(bed, num_reads)
                    covDF$feature <- name
                    labels        <- rbind(labels,
                                           c(0.95*max(log10(covDF$cumSize)),
                                             max(covDF$frip)+0.001,
                                             name, round(max(covDF$frip),2),
                                             plotColors[i]))
                    feature_dist <- rbind(feature_dist,
                        c(name, nrow(bed), as.numeric(sum(abs(bed$V3-bed$V2))),
                          as.numeric((sum(abs(bed$V3-bed$V2))/genome_size))))
                }
            }
        }
    }

    # Reorder by labels
    if (exists("covDF")) {
        covDF$feature <- factor(covDF$feature, levels=(labels$name))
    }

    feature_dist$numbases <- as.numeric(feature_dist$numbases)
    feature_dist$expected <- as.numeric(feature_dist$expected)
    feature_dist$observed <- as.numeric(labels$val)
    feature_dist$logOE <- log10(feature_dist$observed/feature_dist$expected)
    feature_dist <- merge(feature_dist, labels, by.x="feature", by.y="name")
    #feature_dist <- feature_dist[order(feature_dist$logOE, decreasing=TRUE),]
    feature_dist <- feature_dist[order(feature_dist$logOE),]
    rownames(feature_dist) <- NULL
    feature_dist$feature <- factor(feature_dist$feature,
                                   levels=feature_dist$feature)
    feature_dist$color <- factor(feature_dist$color,
                                 levels=feature_dist$color)

    if (!is.null(bedFile)) {

        if (tolower(type) == "both") {
            # Produce plot with bed files
            p <- ggplot(covDF[which(covDF$frip > min(density(covDF$frip)$y)),],
                        aes(x=log10(cumSize), y=frip,
                            group=feature, color=feature)) +
                #geom_line(aes(linetype=feature), size=2, alpha=0.5) +
                geom_line(size=2, alpha=0.5) +
                guides(linetype = FALSE) +
                labs(x="log10(number of bases)", y="FRiF") +
                theme_PEPPRO()

            # Recolor and reposition legend
            p <- p + scale_color_manual(labels=paste0(labels$name, ": ",
                                                      labels$val),
                                        values=labels$color) +
                labs(color="FRiF") +
                theme(legend.position="right",
                      legend.justification=c(0.1,0.9),
                      legend.background=element_blank(),
                      legend.text = element_text(size = rel(0.65)),
                      legend.key = element_blank(),
                      axis.text.x = element_text(angle = 0, hjust = 1,
                                                 vjust=0.5))

            p2 <- ggplot(feature_dist, aes(x = feature, y = logOE)) +
                geom_bar(stat="identity", fill=labels$color, alpha=0.5) + 
                geom_hline(aes(yintercept=0), linetype="dotted") +
                xlab('') +
                ylab('log10(Obs/Exp)') +
                coord_flip() +
                scale_x_discrete(position="top") +
                theme_PEPPRO() +
                theme(plot.background = element_rect(fill = "transparent",
                                                     color = NA,),
                      panel.background = element_rect(fill = "transparent"),
                      rect = element_rect(fill = "transparent"),
                      plot.margin = unit(c(0,0,-6.5,-6.5),"mm"))

            g   <- ggplotGrob(p2)
            min_x <- min(layer_scales(p)$x$range$range)
            max_x <- max(layer_scales(p)$x$range$range)
            min_y <- min(layer_scales(p)$y$range$range)
            max_y <- max(layer_scales(p)$y$range$range)

            p <- p + annotation_custom(grob = g, xmin = 1.05*min_x,
                                       xmax=min_x*2.05, ymin=max_y/2,
                                       ymax=max_y)
        } else if (tolower(type) == "frif") {
            p <- ggplot(covDF[which(covDF$frip > min(density(covDF$frip)$y)),],
                        aes(x=log10(cumSize), y=frip,
                            group=feature, color=feature)) +
                #geom_line(aes(linetype=feature), size=2, alpha=0.5) +
                geom_line(size=2, alpha=0.5) +
                guides(linetype = FALSE) +
                labs(x="log10(number of bases)", y="FRiF") +
                theme_PEPPRO()

            # Recolor and reposition legend
            p <- p + scale_color_manual(labels=paste0(labels$name, ": ",
                                                      labels$val),
                                        values=labels$color) +
                labs(color="FRiF") +
                theme(legend.position=c(0.075,0.975),
                      legend.justification=c(0.1,0.9),
                      legend.title = element_blank(),
                      legend.text = element_text(size = rel(0.65)), 
                      legend.background=element_blank(),
                      legend.key = element_blank(),
                      axis.text.x = element_text(angle = 0, hjust = 1,
                                                 vjust=0.5))
        } else if (tolower(type) == "prif") {
            p <- ggplot(feature_dist, aes(x = feature, y = logOE)) +
                geom_bar(stat="identity",
                         fill = feature_dist$color,
                         alpha = 0.5) + 
                geom_hline(aes(yintercept=0), linetype="dotted") +
                xlab('') +
                ylab('log10(Obs/Exp)') +
                coord_flip() +
                theme_PEPPRO()
        } else {
            #default to both
            # Produce plot with bed files
            p <- ggplot(covDF[which(covDF$frip > min(density(covDF$frip)$y)),],
                        aes(x=log10(cumSize), y=frip,
                            group=feature, color=feature)) +
                geom_line(aes(linetype=feature), size=2, alpha=0.5) +
                guides(linetype = FALSE) +
                labs(x="log10(number of bases)", y="FRiF") +
                theme_PEPPRO()

            # Recolor and reposition legend
            p <- p + scale_color_manual(labels=paste0(labels$name, ": ",
                                                      labels$val),
                                        values=labels$color) +
                labs(color="FRiF") +
                theme(legend.position="right",
                      legend.justification=c(0.1,0.9),
                      legend.background=element_blank(),
                      legend.key = element_blank(),
                      axis.text.x = element_text(angle = 0, hjust = 1,
                                                 vjust=0.5))

            p2 <- ggplot(feature_dist, aes(x = feature, y = logOE)) +
                geom_bar(stat="identity", fill=labels$color, alpha=0.5) + 
                geom_hline(aes(yintercept=0), linetype="dotted") +
                xlab('') +
                ylab('log10(Obs/Exp)') +
                coord_flip() +
                scale_x_discrete(position="top") +
                theme_PEPPRO() +
                theme(plot.background = element_rect(fill = "transparent",
                                                     color = NA,),
                      panel.background = element_rect(fill = "transparent"),
                      rect = element_rect(fill = "transparent"),
                      plot.margin = unit(c(0,0,-6.5,-6.5),"mm"))

            g   <- ggplotGrob(p2)
            min_x <- min(layer_scales(p)$x$range$range)
            max_x <- max(layer_scales(p)$x$range$range)
            min_y <- min(layer_scales(p)$y$range$range)
            max_y <- max(layer_scales(p)$y$range$range)

            p <- p + annotation_custom(grob = g, xmin = 1.05*min_x,
                                       xmax=min_x*2.05, ymin=max_y/2,
                                       ymax=max_y)
        }

        
    } else {
        write("Unable to produce FRiF plot!\n", stdout())
    }

    if (!exists("p")) {
        p <- ggplot()
    }

    return(p)
}


#' The function rounds the up to the nearest "nice" number.
#'
#' From:
#' https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
roundUpNice <- function(x, nice=c(1,2,3,4,5,6,7,8,9,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}


#' Plot TSS enrichment
#'
#' This function plots the global TSS enrichment.
#'
#' @param TSSfile TSS enrichment file
#' @keywords TSS enrichment
#' @export
#' @examples
#' data("tss")
#' plotTSS(TSSfile = "tss")
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
            quit(save = "no", status = 1, runLast = FALSE)
        }
    }

    t1 <- theme_classic(base_size=14) + 
            theme(plot.background  = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border     = element_rect(colour = "black",
                                                  fill=NA, size=0.5),
                  panel.background = element_blank(),
                  axis.line    = element_blank(),
                  legend.position="none",
                  aspect.ratio = 1,
                  axis.ticks.length = unit(2, "mm"))

    iMat <- data.table(V1 = numeric())
    if (length(TSSfile) == 1) {
        if (exists(TSSfile)) {
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
        quit(save = "no", status = 1, runLast = FALSE)
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
        # check for true peak
        if ((normTSS$score[peakPos]/normTSS$score[peakPos-1]) > 1.5 &
            (normTSS$score[peakPos]/normTSS$score[peakPos+1]) > 1.5) {
            tmpTSS  <- normTSS$score[-peakPos]
            peakPos <- which.max(tmpTSS) + 1
        }
        TSSscore <- round(mean(normTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                                       peakPos+50)), score]),1)
        if (is.nan(TSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile, "\n"))
            quit(save = "no", status = 1, runLast = FALSE)
        }
    } else {
        val      <- 0.05*nrow(plus)
        #normTSS  <- (plus / mean(plus[c(1:val,
        #            (nrow(plus)-val):nrow(plus)), V1]))
        normTSS           <- plus / mean(plus[c(1:val), V1])
        colnames(normTSS) <- c("score")
        peakPos  <- which.max(normTSS$score)
        # check for true peak
        if ((normTSS$score[peakPos]/normTSS$score[peakPos-1]) > 1.5 &
            (normTSS$score[peakPos]/normTSS$score[peakPos+1]) > 1.5) {
            tmpTSS  <- normTSS$score[-peakPos]
            peakPos <- which.max(tmpTSS) + 1
        }
        TSSscore <- round(mean(normTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                                       peakPos+50)), score]),1)
        if (is.nan(TSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile[1], "\n"))
            quit(save = "no", status = 1, runLast = FALSE)
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
        # geom_hline(yintercept = 6, linetype = 2,
        #            color = "grey", size = 0.25) +
        geom_smooth(method="loess", span=0.02,
                    se=FALSE, colour=lineColor) +
        labs(x = "Distance from TSS (bp)", y = "TSS Enrichment Score")
    y_max <- roundUpNice(TSSscore)
    p <- pre + t1 +
         scale_x_continuous(expand=c(0,0)) +
         scale_y_continuous(expand=c(0,0)) +
         coord_cartesian(xlim=c(-2300, 2300), ylim=c(0, 1.1*y_max))
    if (exists("minus")) {
        val      <- 0.025*nrow(minus)
        # normTSS  <- (minus / mean(minus[c(1:val,
        #             (nrow(minus)-val):nrow(minus)), V1]))
        minusNormTSS           <- minus / mean(minus[c(1:val), V1])
        colnames(minusNormTSS) <- c("score")
        peakPos       <- which.max(minusNormTSS$score)
        # check for true peak
        if ((minusNormTSS$score[peakPos]/minusNormTSS$score[peakPos-1]) > 1.5 &
            (minusNormTSS$score[peakPos]/minusNormTSS$score[peakPos+1]) > 1.5) {
            tmpTSS  <- minusNormTSS$score[-peakPos]
            peakPos <- which.max(tmpTSS) + 1
        }

        minusTSSscore <- round(
            mean(minusNormTSS[(max(0, peakPos-50)):(min(nrow(minusNormTSS),
                               peakPos+50)), score]),1)
        if (is.nan(minusTSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile[2], "\n"))
            quit(save = "no", status = 1, runLast = FALSE)
        }
        p <- p + geom_smooth(data=minusNormTSS,
                             aes(x=(as.numeric(rownames(minusNormTSS))-
                                   (nrow(minusNormTSS)/2)),
                                 y=score, group=1, colour="black"),
                             method="loess", span=0.02,
                             se=FALSE, colour="blue") +
                annotate("rect", xmin=1200, xmax=2300, ymin=0.9*y_max,
                         ymax=1.1*y_max, fill="gray95") +
                annotate("text", x=1750, y=1.05*y_max, label="TSS Score",
                         fontface = 1, hjust=0.5) +
                annotate("text", x=1500, y=y_max, label="+", fontface = 2,
                          hjust=0.5, color=lineColor) +
                annotate("text", x=1500, y=0.95*y_max, label=TSSscore,
                         fontface = 2,  hjust=0.5, color=lineColor) +
                annotate("text", x=2000, y=y_max, label="-",
                         fontface = 2,  hjust=0.5, color="blue") +
                annotate("text", x=2000, y=0.95*y_max, label=minusTSSscore,
                         fontface = 2,  hjust=0.5, color="blue")
    } else {
        p <- p + annotate("rect", xmin=1200, xmax=2300, ymin=0.9*y_max,
                          ymax=1.1*y_max, fill="gray95") +
                 annotate("text", x=1750, y=1.05*y_max, label="TSS Score",
                          fontface = 1, hjust=0.5) +
                 annotate("text", x=1750, y=0.95*y_max, label=TSSscore,
                          fontface = 2, hjust=0.5)
    }

    return(p)
}


#' Derive the sample name from input file and return with full path
#'
#' @param path A path to a file for which you wish to extract the sample name
#' @param num_fields An integer representing the number of fields to strip
#' @param delim A delimiter for the fields splitting a path or string
sampleName <- function(path, num_fields=2, delim='_') {
    name <- basename(tools::file_path_sans_ext(path))
    if(num_fields == 0) {return(name)}
    for(n in 1:num_fields) name <- gsub(paste0(delim, "[^", delim, "]*$"), "", name)
    return(paste(dirname(path), name, sep="/"))
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
#' data("frag_len")
#' data("frag_len_count")
#' plotFLD(fragL = "frag_len", fragL_count = "frag_len_count",
#'         fragL_txt = "fragLenDistribution_example.txt", max_fragment=200)
#' @export
plotFLD <- function(fragL,
                    fragL_count,
                    fragL_txt="fragLenDistribution.txt",
                    max_fragment = 200) {

    if (exists(fragL_count)) {
        dat <- data.table(get(fragL_count))
    } else if (file.exists(fragL_count)) {
        dat <- fread(fragL_count)
    } else {
        stop(paste0("FileExistsError: ", fragL_count, " could not be found."))
        quit(save = "no", status = 1, runLast = FALSE)
    }

    if (exists(fragL)) {
        summary_table <- data.table(get(fragL))
    } else if (file.exists(fragL)) {
        summary_table <- fread(fragL)
    } else {
        stop(paste0("FileExistsError: ", fragL, " could not be found."))
        quit(save = "no", status = 1, runLast = FALSE)
    }

    dat1 <- dat[dat$V2<=max_fragment,]
    tmp  <- seq(1:as.numeric(dat1[1,2]-1))
    dat0 <- data.table(V1=rep(0,length(tmp)),V2=tmp)
    dat2 <- rbind(dat0, dat1)

    x_min = which.min(dat1$V1[1:which.max(dat1$V1)])

    p <- ggplot(dat1[x_min:nrow(dat1),], aes(x=V2, y=V1)) +
            geom_point(size=1, alpha=0.25) +
            geom_line(alpha=0.5) +
            annotate("rect", xmin=-Inf, xmax=20, ymin=-Inf, ymax=Inf,
                 alpha=0.1, fill="#ff001e") +
            annotate("text", x=25, y=(max(dat1$V1)/2),
                     size=theme_get()$text[["size"]]/4,
                     label="partial degradation", angle=90, col="#858585") +
            annotate("rect", xmin=-Inf, xmax=30, ymin=-Inf, ymax=Inf,
                     alpha=0.1, fill="#ffee00") + 
            annotate("text", x=7.5, y=(max(dat1$V1)/2),
                     size=theme_get()$text[["size"]]/4,
                     label="high degradation", angle=90, col="#858585") +
            xlab("fragment length") + 
            ylab("number of reads") +
            theme_PEPPRO()

    summ <- data.table(Min=min(summary_table$V1),
                       Max=max(summary_table$V1),
                       Median=median(summary_table$V1),
                       Mean=mean(summary_table$V1),
                       Stdev=sd(summary_table$V1))
    # Write summary table to stats file
    fwrite(summ, file=fragL_txt, row.names=F, quote=F, sep="\t")

    return(p)
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
    textNums[!is.na(n)] <- paste0(base, exponents)
    textReturn          <- parse(text=textNums)
    return(textReturn)
}


#' Determine which quantile to use as cutoff.
#'
#' Modified from: GenomicDistributions (Tessa Danehy)
#'
#' @param vec A vector of numbers.
#' @param baseline A minimum quantile cutoff.
#' @keywords cutoff quantiles
#' @examples
#' calcQuantileCutoff()
calcQuantileCutoff = function(vec, baseline=1){
  n <- length(vec) # number of observations
  if (n > 1000) {n = 1000}
  if (n < 100) {n = 100}
  q <- max(baseline, (11 - round(n/100))) # finding quantiles for flanking bins
  return(q)
}


#' Calculate histogram binning by quantile.
#'
#' Modified from: GenomicDistributions (Tessa Danehy)
#'
#' @param vec A vector of numbers.
#' @param q The quantile to use for cutoffs.
#' @param bins The number of bins to use for a histogram.
#' @param transformed Adjust divisions if transformed
#' @keywords quantiles
#' @examples
#' calcDivisions()
calcDivisions = function(vec, q, bins=NULL, transformed=FALSE){
  q = as.numeric(q)
  if(q > 50){
    message("Quantile should not be larger than 50. Optimal size is under 10.")
    q <- 50
  }
  if(!is.null(bins)){
    b <- bins + 1
  }
  else {
    b <- abs((30-(2*q))/q) # finding the number of bins based on the quantiles
  }
  quant  <- unname(quantile(vec, probs = c((q/100), (1-(q/100)))))
  seq_10 <- seq(quant[1], quant[2], length = b)
  if (transformed) {
    div  <- c(-Inf, round(seq_10, 2), Inf)
  } else {
    div  <- c(-Inf, round(seq_10), Inf)
  }
  
  return(div)
}


#' Create a character vector of bin labels.
#'
#' From: GenomicDistributions (Tessa Danehy)
#'
#' @param breakPoints A vector of numbers.
#' @param digits The number of digits to round to.
#' @param collapse The character to separate paste values by.
#' @param infBins Whether to use infinite bins.
#' @keywords labels
#' @examples
#' labelCuts()
labelCuts = function(breakPoints, digits=1, collapse="-", infBins=FALSE) {
    labels <- 
        apply(round(cbind(breakPoints[-length(breakPoints)],   
                          breakPoints[-1]),digits), 1, paste0,
                          collapse=collapse) 
    
    if (infBins) {
        labels[1] <- paste0("<", breakPoints[2])
        labels[length(labels)] <- paste0(">", breakPoints[length(breakPoints)-1])
    }
    return(labels)
}


#' Create a character vector of bin labels.
#'
#' Modified from: GenomicDistributions (Tessa Danehy)
#'
#' @param vec A vector (or list of vectors) of numbers.
#' @param divisions The break points of the distribution
#' @keywords quantiles
#' @examples
#' cutDists()
cutDists = function(vec, divisions = c(-Inf, -1e6, -1e4, -1000, -100, 0,
                                         100, 1000, 10000, 1e6, Inf)) {
    if (is.list(vec)) {
        x = lapply(vec, cutDists)
        
        # To accommodate multiple lists, we'll need to introduce a new 'name'
        # column to distinguish them.
        nameList = names(vec)
        if(is.null(nameList)) {
            nameList = 1:length(query) # Fallback to sequential numbers
        }
        
        # Append names
        xb = rbindlist(x)
        xb$name = rep(nameList, sapply(x, nrow))
        
        return(xb)
    }
    divisions <- unique(divisions)
    labels <- labelCuts(signif(divisions, 3), collapse=" to ", infBins=TRUE)
    #message(paste0("breaks: ", paste0(divisions, collapse=" ")))
    cuts   <- cut(vec, divisions, labels)
    return(as.data.frame(table(cuts)))
}


#' Plot the distribution of genic exonRPKM/intronRPKM ratios
#'
#' This function plots the distribution of by gene exon RPKM divided by
#' intron RPKM ratios. Can produce raw or log10 distributions, but reports
#' both median values.
#'
#' @param rpkm A three column TSV format file containing
#'             "gene", "intron RPKM", "exon RPKM" columns.
#' @param name X-axis label, typically a sample name
#' @param raw Plot raw distribution
#' @param type Plot format
#' @param annotate Display mean and median values on plot
#' @keywords mRNA contamination
#' @export
#' @examples
#' data("rpkm_ratios")
#' mRNAcontamination(rpkm = "rpkm_ratios")
#' @export
mRNAcontamination <- function(rpkm,
                              name='mRNA contamination ratios',
                              raw=TRUE,
                              type=c("histogram", "boxplot", "violin"),
                              annotate=TRUE) {
    if (exists(rpkm)) {
        RPKM <- data.table(get(rpkm))
    } else if (file.exists(rpkm)) {
        RPKM <- fread(rpkm)
    } else {
        stop(paste0("FileExistsError: ", rpkm, " could not be found."))
        quit(save = "no", status = 1, runLast = FALSE)
    }
    colnames(RPKM) <- c("chr", "start", "end", "gene","ratio","strand")

    finite_rpkm <- RPKM[is.finite(RPKM$ratio),]

    if (raw) {
        div <- calcDivisions(finite_rpkm$ratio,
                             calcQuantileCutoff(finite_rpkm$ratio,
                                                baseline = 3))
    } else {
        div <- calcDivisions(log10(finite_rpkm$ratio),
                             calcQuantileCutoff(log10(finite_rpkm$ratio),
                                                baseline = 3),
                             transformed=TRUE)
    }

    # ensure breaks are not duplicated
    div        <- unique(div)
    quantLabel <- paste(calcQuantileCutoff(finite_rpkm$ratio),"%", sep='')

    if (raw) {
        if (type == "histogram") {
            if (length(div) <= 3) {
                base_plot  <- ggplot(data = finite_rpkm, aes(x=ratio))  
            } else {
                # calculate a frequency table with the specified divisions
                rpkm_table <- cutDists(finite_rpkm$ratio, divisions = div)
                base_plot  <- ggplot(data = rpkm_table,  aes(x=cuts, y=Freq))
            }
            
        } else {
            base_plot  <- ggplot(data = finite_rpkm,  aes(x="", y=(ratio)))
        }
    } else {
        if (type == "histogram") {
            if (length(div) <= 3) {
                base_plot  <- ggplot(data = finite_rpkm, aes(x=log10(ratio)))  
            } else {
                # calculate a frequency table with the specified divisions
                rpkm_table <- cutDists(log10(finite_rpkm$ratio),
                                       divisions = div)
                base_plot  <- ggplot(data = rpkm_table,  aes(x=cuts, y=Freq))
            }
        } else {
            base_plot  <- ggplot(data = finite_rpkm, aes(x="", y=log10(ratio)))
        }
    }

    if (type == "histogram") {
        if (raw) {
            if (length(div) <= 3) {
                plot <- base_plot +
                    geom_histogram(col="black", fill=I("transparent")) +
                    geom_vline(aes(xintercept=median(ratio)),
                               color="gray", linetype="dashed", size=1) +
                    annotate("text", x=median(finite_rpkm$ratio),
                             y=(ceiling(quantile(finite_rpkm$ratio, 0.25))),
                                label="median", angle=90,
                                color="gray", vjust=-0.5) +
                    geom_vline(aes(xintercept=mean(ratio)),
                               color="light gray", linetype="dotted", size=1) +
                    annotate("text", x=mean(finite_rpkm$ratio),
                             y=(ceiling(quantile(finite_rpkm$ratio, 0.25))),
                                label="mean", angle=90,
                                color="light gray", vjust=-0.5) +
                    labs(x=expression((over(exon[RPKM], intron[RPKM]))~X~Gene),
                         y="frequency") +
                    xlim(c(0, ceiling(quantile(finite_rpkm$ratio, 0.90)))) +
                    theme_PEPPRO()
            } else {
                plot <- base_plot +
                    geom_bar(stat="identity",
                             fill = c("maroon",
                                      rep("gray", (length(div)-3)),
                                      "maroon")) + 
                    labs(x=expression((over(exon[RPKM], intron[RPKM]))~X~Gene),
                         y="frequency") +
                    geom_text(aes(label= quantLabel),
                              data=rpkm_table[c(1,length(rpkm_table$Freq)),],
                              vjust=-1)
            }
        } else {
            if (length(div) <= 3) {
                plot = base_plot +
                    geom_histogram(col="black", fill=I("transparent")) +
                    geom_vline(aes(xintercept=median(log10(ratio))),
                               color="gray", linetype="dashed", size=1) +
                    geom_vline(aes(xintercept=mean(log10(ratio))),
                               color="light gray", linetype="dotted", size=1) +
                    labs(x=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene)) +
                    scale_x_log10(limits = c(0.001, 50),
                                  expand = expand_scale(mult = c(0, 0)),
                                  labels=fancyNumbers,
                                  breaks=prettyLogs) +
                    annotation_logticks(sides = c("rl"))
            } else {
                plot <- base_plot +
                    geom_bar(stat="identity",
                             fill = c("maroon",
                                      rep("gray", (length(div)-3)),
                                      "maroon")) + 
                    labs(x=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene),
                         y="frequency") +
                    geom_text(aes(label= quantLabel),
                              data=rpkm_table[c(1,length(rpkm_table$Freq)),],
                              vjust=-1)
            }
        }
    } else if (type == "boxplot") {
        if (raw) {
            plot = base_plot +
                    stat_boxplot(geom ='errorbar', width = 0.25) +
                    geom_boxplot(width = 0.25,
                                 outlier.color='red',
                                 outlier.shape=1) +
                    stat_summary(fun.y = "mean", geom = "point",
                                 shape = 1, size = 2) +
                    labs(x=name,
                         y=expression((over(exon[RPKM], intron[RPKM]))~X~Gene)) +
                    ylim(c(0, ceiling(quantile(finite_rpkm$ratio, 0.90))))
        } else {
            plot <- base_plot +
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
                     y=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene))
        }
    } else if (type == "violin") {
        if (raw) {
            plot = base_plot +
                #stat_boxplot(geom ='errorbar', width = 0.25) +
                geom_violin(width = 0.25, draw_quantiles = c(0.25,0.75),
                            linetype="dashed") +
                geom_violin(width=0.25, fill="transparent",
                            draw_quantiles = 0.5) +
                stat_summary(fun.y = "mean", geom = "point",
                             shape = 1, size = 2) +
                labs(x=name,
                     y=expression((over(exon[RPKM], intron[RPKM]))~X~Gene)) +
                ylim(c(0, ceiling(quantile(finite_rpkm$ratio, 0.90))))
        } else {
            plot <- base_plot +
                #stat_boxplot(geom ='errorbar', width = 0.25) +
                geom_violin(width = 0.25, draw_quantiles = c(0.25,0.75),
                            linetype="dashed") +
                geom_violin(width=0.25, fill="transparent",
                            draw_quantiles = 0.5) +
                stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
                stat_summary(fun.y = "mean", geom = "point",
                             shape = 1, size = 2) +
                scale_y_log10(limits = c(0.001, 50),
                              expand = expand_scale(mult = c(0, 0)),
                              labels=fancyNumbers,
                              breaks=prettyLogs) +
                annotation_logticks(sides = c("rl")) +
                labs(x=name,
                     y=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene))
        }
    } else {
        # Default to histogram
        if (raw) {
            if (length(div) <= 3) {
                plot <- base_plot +
                    geom_histogram(col="black", fill=I("transparent")) +
                    geom_vline(aes(xintercept=median(ratio)),
                               color="gray", linetype="dashed", size=1) +
                    annotate("text", x=median(finite_rpkm$ratio),
                             y=(ceiling(quantile(finite_rpkm$ratio, 0.25))),
                                label="median", angle=90,
                                color="gray", vjust=-0.5) +
                    geom_vline(aes(xintercept=mean(ratio)),
                               color="light gray", linetype="dotted", size=1) +
                    annotate("text", x=mean(finite_rpkm$ratio),
                             y=(ceiling(quantile(finite_rpkm$ratio, 0.25))),
                                label="mean", angle=90,
                                color="light gray", vjust=-0.5) +
                    labs(x=expression((over(exon[RPKM], intron[RPKM]))~X~Gene),
                         y="frequency") +
                    xlim(c(0, ceiling(quantile(finite_rpkm$ratio, 0.90)))) +
                    theme_PEPPRO()
            } else {
                plot <- base_plot +
                    geom_bar(stat="identity",
                             fill = c("maroon",
                                      rep("gray", (length(div)-3)),
                                      "maroon")) + 
                    labs(x=expression((over(exon[RPKM], intron[RPKM]))~X~Gene),
                         y="frequency") +
                    geom_text(aes(label= quantLabel),
                              data=rpkm_table[c(1,length(rpkm_table$Freq)),],
                              vjust=-1)
            }
        } else {
            if (length(div) <= 3) {
                plot = base_plot +
                    geom_histogram(col="black", fill=I("transparent")) +
                    geom_vline(aes(xintercept=median(log10(ratio))),
                               color="gray", linetype="dashed", size=1) +
                    geom_vline(aes(xintercept=mean(log10(ratio))),
                               color="light gray", linetype="dotted", size=1) +
                    labs(x=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene)) +
                    scale_x_log10(limits = c(0.001, 50),
                                  expand = expand_scale(mult = c(0, 0)),
                                  labels=fancyNumbers,
                                  breaks=prettyLogs) +
                    annotation_logticks(sides = c("rl"))
            } else {
                plot <- base_plot +
                    geom_bar(stat="identity",
                             fill = c("maroon",
                                      rep("gray", (length(div)-3)),
                                      "maroon")) + 
                    labs(x=expression(log[10](over(exon[RPKM], intron[RPKM]))~X~Gene),
                         y="frequency") +
                    geom_text(aes(label= quantLabel),
                              data=rpkm_table[c(1,length(rpkm_table$Freq)),],
                              vjust=-1)
            }
        }
    }

    label1 <- paste("'median'[log[10]]", ":~",
                    round(median(log10((finite_rpkm$ratio))), 2))
    label2 <- paste("'median'[raw]", ":", round(median(finite_rpkm$ratio), 2))

    if (type == "histogram") {
        max_x <- length(layer_scales(plot)$x$range$range)
    } else {
        max_x <- suppressMessages(
            suppressWarnings(layer_scales(plot)$x$range$range[2]))
    }
    if (is.na(max_x)) {max_x <- Inf}
    max_y <- suppressMessages(
        suppressWarnings(layer_scales(plot)$y$range$range[2]))

    # TODO: make summary plot of these that IS boxplots
    if (annotate) {
        q <- plot + annotate("text", x = floor(max_x), y = floor(max_y),
                         hjust="right", vjust=1.05, 
                         label = label1, parse=TRUE) +
            annotate("text", x = floor(max_x), y = floor(max_y),
                     hjust="right", vjust=2.05,
                     label = label2, parse=TRUE) +
            theme_PEPPRO()
    } else {
        q <- plot + theme_PEPPRO()
    }
    

    return(q)
}


#' Plot the distribution of highest covered TSS density/gene body density ratios
#'
#' @param pi A single column containing the ratio of TSS densities/gene body
#'           densities for the highest scoring TSSs
#' @param name X-axis label, typically a sample name
#' @param type Plot format
#' @param annotate Display mean and median values on plot
#' @keywords pause index
#' @export
#' @examples
#' data("pidx")
#' plotPI(pi = "pidx")
#' @export
plotPI <- function(pi, name='pause indicies',
                   type=c("histogram", "boxplot", "violin"),
                   annotate=TRUE) {
    # TODO: make summary plot of these that IS boxplots
    if (exists(pi)) {
        PI <- data.table(get(pi))
    } else if (file.exists(pi)) {
        PI <- fread(pi)
    } else {
        stop(paste0("FileExistsError: ", pi, " could not be found."))
        quit(save = "no", status = 1, runLast = FALSE)
    }
    colnames(PI) <- c("chr", "start", "end", "name", "pi", "strand")

    div <- calcDivisions(PI$pi, calcQuantileCutoff(PI$pi, baseline = 3))
    quantLabel <- paste(calcQuantileCutoff(PI$pi, baseline = 3),"%", sep='')

    if (type == "histogram") {
        if (length(div) <= 3) {
            base_plot <- ggplot(data = PI, aes(x=pi))
        } else {
            # calculate a frequency table with the specified divisions
            pi_table  <- cutDists(PI$pi, divisions = div)
            base_plot <- ggplot(data = pi_table,  aes(x=cuts, y=Freq))
        }
        
    } else {
        base_plot <- ggplot(data = PI,  aes(x="", y=pi))
    }

    if (type == "histogram") {
        if (length(div) <= 3) {
            q <- base_plot +
                    geom_histogram(col="black", fill=I("transparent")) +
                    geom_vline(aes(xintercept=median(PI$pi)),
                               color="gray", linetype="dashed", size=1) +
                    annotate("text", x=median(PI$pi),
                             y=(ceiling(quantile(PI$pi, 0.25))),
                                label="median", angle=90,
                                color="gray", vjust=-0.5) +
                    geom_vline(aes(xintercept=mean(PI$pi)),
                               color="light gray", linetype="dotted", size=1) +
                    annotate("text", x=mean(PI$pi),
                             y=(ceiling(quantile(PI$pi, 0.25))),
                                label="mean", angle=90,
                                color="light gray", vjust=-0.5) +
                    labs(x="pause indicies", y="frequency") +
                    xlim(c(0, ceiling(quantile(PI$pi, 0.90)))) +
                    theme_PEPPRO()
        } else {
            q <- base_plot +
                geom_bar(stat="identity",
                         fill = c("maroon",
                                  rep("gray", (length(div)-3)),
                                  "maroon")) + 
                labs(x="pause indicies", y="frequency") +
                geom_text(aes(label= quantLabel),
                          data=pi_table[c(1,length(pi_table$Freq)),],
                          vjust=-1)
        }
    } else if (type == "boxplot") {
        plot <- base_plot +
                stat_boxplot(geom ='errorbar', width = 0.25) +
                geom_boxplot(width = 0.25,
                             outlier.color='red',
                             outlier.shape=1) +
                stat_summary(fun.y = "mean", geom = "point",
                             shape = 1, size = 2) +
                labs(x=name, y="each gene's pause index")
    } else if (type == "violin") {
        plot <- base_plot +
                geom_violin(width = 0.25, draw_quantiles = c(0.25,0.75),
                            linetype="dashed") +
                geom_violin(width=0.25, fill="transparent",
                            draw_quantiles = 0.5) +
                stat_summary(fun.y = "mean", geom = "point",
                             shape = 1, size = 2) +
                labs(x=name, y="each gene's pause index")
    } else {
        # default to histogram
         if (length(div) <= 3) {
            q <- base_plot +
                    geom_histogram(col="black", fill=I("transparent")) +
                    geom_vline(aes(xintercept=median(PI$pi)),
                               color="gray", linetype="dashed", size=1) +
                    annotate("text", x=median(PI$pi),
                             y=(ceiling(quantile(PI$pi, 0.25))),
                                label="median", angle=90,
                                color="gray", vjust=-0.5) +
                    geom_vline(aes(xintercept=mean(PI$pi)),
                               color="light gray", linetype="dotted", size=1) +
                    annotate("text", x=mean(PI$pi),
                             y=(ceiling(quantile(PI$pi, 0.25))),
                                label="mean", angle=90,
                                color="light gray", vjust=-0.5) +
                    labs(x="pause indicies", y="frequency") +
                    xlim(c(0, ceiling(quantile(PI$pi, 0.90)))) +
                    theme_PEPPRO()
        } else {
            q <- base_plot +
                geom_bar(stat="identity",
                         fill = c("maroon",
                                  rep("gray", (length(div)-3)),
                                  "maroon")) + 
                labs(x="pause indicies", y="frequency") +
                geom_text(aes(label= quantLabel),
                          data=pi_table[c(1,length(pi_table$Freq)),],
                          vjust=-1)
        }
    }  

    if (type != "histogram") {
        if (max(PI$pi) > 500) {
            q <- plot + scale_y_continuous(breaks = round(seq(min(PI$pi),
                                                  max(PI$pi),
                                                  by = 50), 0),
                                        limits=c(0, max(PI$pi)))
        } else if (max(PI$pi) > 100 & max(PI$pi) < 500) {
            q <- plot + scale_y_continuous(breaks = round(seq(min(PI$pi),
                                                  max(PI$pi),
                                                  by = 25), 0),
                                        limits=c(0, max(PI$pi)))
        } else {
            q <- plot + scale_y_continuous(breaks = round(seq(min(PI$pi),
                                                  max(PI$pi),
                                                  by = 5), 0),
                                        limits=c(0, max(PI$pi)))
        }
        q <- q + coord_cartesian(ylim=c(0, ceiling(boxplot(PI$pi)$stats[5]))) +
            theme_PEPPRO()
        max_x <- suppressMessages(
            suppressWarnings(layer_scales(q)$x$range$range[2]))
        max_y  <- ceiling(boxplot(PI$pi)$stats[5])
    } else {
        max_x <- length(layer_scales(q)$x$range$range)
        max_y <- suppressMessages(
            suppressWarnings(layer_scales(q)$y$range$range[2]))
    }

    if (is.na(max_x)) {max_x <- Inf}
    
    label1 <- paste("'median'", ":", round(median(PI$pi), 2))
    label2 <- paste("'mean'", ":", round(mean(PI$pi), 2))
    if (annotate) {
        q <- q + annotate("text", x = floor(max_x), y = floor(max_y),
                      hjust="right", vjust=1.05, label = label1, parse=TRUE) +
         annotate("text", x = floor(max_x), y = floor(max_y),
                      hjust="right", vjust=2.15, label = label2, parse=TRUE) +
         theme_PEPPRO()
    } else {
        q <- q + theme_PEPPRO()   
    }

    return(q)
}

#' Calculate mode(s) of data
#'
#' From: https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
#' @param x A vector of numbers or characters
#' @param return_multiple Bool to return multiple modes or first in order
#' @param na.rm Bool Remove NAs
#'
#' @keywords mode
mode <- function(x, return_multiple = TRUE, na.rm = FALSE) {
    if(na.rm){
        x <- na.omit(x)
    }
    ux       <- unique(x)
    freq     <- tabulate(match(x, ux))
    mode_loc <- if(return_multiple) which(freq==max(freq)) else which.max(freq)
    return(ux[mode_loc])
}


#' Determine the appropriate abbreviation for large numbers
#'
#' Modified From: https://stackoverflow.com/questions/28159936/formatting-large-currency-or-dollar-values-to-millions-billions
#' @param vec A vector of numbers
#'
#' @keywords abbreviation
getAbbr <- function(vec) { 
    div <- findInterval(as.numeric(gsub("\\,", "", vec)), c(0, 1e3, 1e6, 1e9, 1e12) )
    return(paste(c("","K","M","B","T")[mode(div)]))
}


#' Determine the appropriate dividing factor for large numbers
#'
#' Modified From: https://stackoverflow.com/questions/28159936/formatting-large-currency-or-dollar-values-to-millions-billions
#' @param vec A vector of numbers
#'
#' @keywords abbreviation
getFactor <- function(vec) { 
    div <- findInterval(as.numeric(gsub("\\,", "", vec)), c(0, 1e3, 1e6, 1e9, 1e12) )
    return(as.numeric(paste(c(1, 1e3, 1e6, 1e9, 1e12)[mode(div)])))
}


#' Plot the cutadapt-based distribution of adapter insertions
#'
#' @param input A cutadapt report
#' @param name A sample name or identifier for the plot title
#' @param umi_len The UMI length
#'
#' @keywords cutadapt
#' @export
#' @examples
#' data("cutadapt")
#' plotCutadapt(input = "cutadapt")
#' @export
plotCutadapt <- function(input, name='cutadapt',
                         umi_len = 0,
                         count_factor = 1000000) {
    if (exists(input)) {
        report <- data.table(get(input))
    } else if (file.exists(input)) {
        report <- fread(input)
    } else {
        stop(paste0("FileExistsError: ", input, " could not be found."))
        quit(save = "no", status = 1, runLast = FALSE)
    }
    
    if (umi_len > 0) {
        report <- report[-(which(report$length == (max(report$length) - umi_len))),]
    }

    # only keep sizes where the expected count represents less than 1% of 
    # the actual count
    report <- report[which(report$expect/report$count < 0.01),]
    
    # inverse length to get ascending order
    report$length <- max(report$length)-report$length

    # don't include size 0 insertions
    report <- report[-nrow(report),]

    abbr <- getAbbr(report$count)
    if (abbr == '') {
        ylabel <- "Number of reads"
    } else {
        ylabel <- paste0("Number of reads (", abbr, ")")
    }

    count_factor <- getFactor(report$count)

    if (20 %in% report$length) {
        degraded_upper <- 20
        degraded_lower <- 10
    } else {
        degraded_upper <- min(report$length) + 10
        degraded_lower <- max(1, degraded_upper - 10)
    }

    if (40 %in% report$length) {
        intact_upper <- 40
        intact_lower <- 30
    } else {
        intact_upper <- max(report$length)
        intact_lower <- max(1, intact_upper - 10)
    }
    
    degradation  <- sum(report[which(report$length >= degraded_lower &
                                     report$length <= degraded_upper),]$count) /
                    max(1, sum(report[which(report$length >= intact_lower &
                               report$length <= intact_upper),]$count))

    q <- ggplot(report, aes(x=max(length)-length, y=count/count_factor)) +
            geom_point() +
            geom_vline(xintercept = 20, linetype = "dotted", alpha=0.25) +
            geom_vline(xintercept = 30, linetype = "longdash", alpha=0.5) +
            labs(title=name, x="Size of insertion", y=ylabel) +
            theme_PEPPRO()
    q <- q + 
        annotate("rect", xmin=-Inf, xmax=20, ymin=-Inf, ymax=Inf,
                 alpha=0.1, fill="#ff001e") +
        annotate("text", x=25, y=(max(report$count/count_factor)/2),
                 size=theme_get()$text[["size"]]/4,
                 label="partial degradation", angle=90, col="#858585") +
        annotate("rect", xmin=-Inf, xmax=30, ymin=-Inf, ymax=Inf,
                 alpha=0.1, fill="#ffee00") + 
        annotate("text", x=7.5, y=(max(report$count/count_factor)/2),
                 size=theme_get()$text[["size"]]/4,
                 label="high degradation", angle=90, col="#858585") +
        annotate("text", x=Inf, y=(max(report$count/count_factor)*0.99),
                 size=theme_get()$text[["size"]]/3, hjust=1.1,
                 label=paste0("degradation ratio: ", round(degradation, 2)))

    return(q)
}


#' Plot the distribution of adapter insertions
#'
#' @param input FLASH histogram output
#' @param name A sample name or identifier for the plot title
#' @param umi_len The UMI length
#'
#' @keywords cutadapt
#' @export
#' @examples
#' data("adapt")
#' plotAdapt(input = "adapt")
#' @export
plotAdapt <- function(input, name='adapt', umi_len = 0) {
    if (exists(input)) {
        report <- data.table(get(input))
    } else if (file.exists(input)) {
        report <- fread(input)
    } else {
        stop(paste0("FileExistsError: ", input, " could not be found."))
        quit(save = "no", status = 1, runLast = FALSE)
    }
    
    colnames(report) <- c("length", "count")

    if (umi_len > 0) {
        report$length <- report$length - umi_len
    }
    
    # don't include size 0 insertions
    report <- report[which(report$length > 0),]

    abbr <- getAbbr(report$count)
    if (abbr == '') {
        ylabel <- "Number of reads"
    } else {
        ylabel <- paste0("Number of reads (", abbr, ")")
    }

    count_factor <- getFactor(report$count)

    if (20 %in% report$length) {
        degraded_upper <- 20
        degraded_lower <- 10
    } else {
        degraded_upper <- min(report$length) + 10
        degraded_lower <- max(1, degraded_upper - 10)
    }

    if (40 %in% report$length) {
        intact_upper <- 40
        intact_lower <- 30
    } else {
        intact_upper <- max(report$length)
        intact_lower <- max(1, intact_upper - 10)
    }
    
    degradation  <- sum(report[which(report$length >= degraded_lower &
                                     report$length <= degraded_upper),]$count) /
                    max(1, sum(report[which(report$length >= intact_lower &
                               report$length <= intact_upper),]$count))

    q <- ggplot(report, aes(x=length, y=count/count_factor)) +
            geom_point() +
            geom_vline(xintercept = 20, linetype = "dotted", alpha=0.25) +
            geom_vline(xintercept = 30, linetype = "longdash", alpha=0.5) +
            labs(title=name, x="Size of insertion", y=ylabel) +
            theme_PEPPRO()
    q <- q + 
        annotate("rect", xmin=-Inf, xmax=20, ymin=-Inf, ymax=Inf,
                 alpha=0.1, fill="#ff001e") +
        annotate("text", x=25, y=(max(report$count/count_factor)/2),
                 size=theme_get()$text[["size"]]/4,
                 label="partial degradation", angle=90, col="#858585") +
        annotate("rect", xmin=-Inf, xmax=30, ymin=-Inf, ymax=Inf,
                 alpha=0.1, fill="#ffee00") + 
        annotate("text", x=7.5, y=(max(report$count/count_factor)/2),
                 size=theme_get()$text[["size"]]/4,
                 label="high degradation", angle=90, col="#858585") +
        annotate("text", x=Inf, y=(max(report$count/count_factor)*0.99),
                 size=theme_get()$text[["size"]]/3, hjust=1.1,
                 label=paste0("degradation ratio: ", round(degradation, 2))) 

    return(q)
}


################################################################################
