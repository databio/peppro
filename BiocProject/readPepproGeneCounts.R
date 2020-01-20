readPepproGeneCounts = function(project) {
    cwd            <- getwd()
    project_dir    <- pepr::config(project)$metadata$output_dir
    sample_names   <- pepr::samples(project)$sample_name
    genomes        <- as.list(pepr::samples(project)$genome)
    names(genomes) <- sample_names
    paths          <- vector("list", length(sample_names))
    names(paths)   <- sample_names

    for (sample in sample_names) {
        paths[[sample]] <- paste(project_dir, 'results_pipeline', sample,
                                 paste0('signal_', genomes[[sample]]),
                                 paste0(sample, "_gene_coverage.bed"), sep="/")
    }

    result <- lapply(paths, function(x){
        #message(paste0("x: ", x))
        if (file.exists(x)) {
            df <- read.table(x)
            colnames(df) <- c('chr', 'start', 'end', 'geneName',
                              'score', 'strand', 'count')
            gr <- GenomicRanges::GRanges(df) 
        } else {
            gr <- GenomicRanges::GRanges() 
        }
    })

    setwd(cwd)
    #names(result) <- sample_names
    return(GenomicRanges::GRangesList(Filter(length, result)))
}


