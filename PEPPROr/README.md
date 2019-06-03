# PEPPROr: Portable Encapsulated Project PRO-seq analysis package in R

**PEPPROr** includes analysis and plotting functions for PRO-seq analysis using the PEPPRO pipeline.

## Installing PEPPROr
PEPPROr may be installed from Github:

```
devtools::install_github("databio/peppro", subdir="PEPPROr")
```

or locally after downloading/cloning the source code:

```
install.packages("path/to/PEPPROr/directory", repos=NULL, type="source")
```

## Usage
```
Usage:   PEPPRO.R [command] {args}
Version: 0.2

Command: preseq	 plot preseq complexity curves
	     frif	 plot fraction of reads in features
	     tss	 plot TSS enrichment
```

To learn more about specific arguments, call the subcommand either by itself (e.g. `PEPPRO.R frif`), or followed by `--help` or `-h` (e.g. `PEPPRO.R frif --help`).
