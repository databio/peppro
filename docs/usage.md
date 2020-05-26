# <img src="../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px">  usage reference

`PEPPRO` command-line usage instructions:



`python pipelines/peppro.py --help`
```{console}
usage: peppro.py [-h] [-R] [-N] [-D] [-F] [-T] [--silent] [--verbosity V]
                 [--logdev] [-C CONFIG_FILE] -O PARENT_OUTPUT_FOLDER
                 [-M MEMORY_LIMIT] [-P NUMBER_OF_CORES] -S SAMPLE_NAME -I
                 INPUT_FILES [INPUT_FILES ...]
                 [-I2 [INPUT_FILES2 [INPUT_FILES2 ...]]] -G GENOME_ASSEMBLY
                 [-Q SINGLE_OR_PAIRED]
                 [--protocol {PRO,pro,PRO-SEQ,PRO-seq,proseq,PROSEQ,GRO,gro,groseq,GROSEQ,GRO-SEQ,GRO-seq}]
                 [--adapter-tool {cutadapt,fastp}]
                 [--dedup-tool {seqkit,fqdedup}]
                 [--trimmer-tool {seqtk,fastx}] [--umi-len UMI_LEN]
                 [--max-len MAX_LEN] [--sob] [--scale]
                 [--prealignments PREALIGNMENTS [PREALIGNMENTS ...]]
                 [--TSS-name TSS_NAME] [--pi-tss ENSEMBL_TSS]
                 [--pi-body ENSEMBL_GENE_BODY] [--pre-name PRE_NAME]
                 [--anno-name ANNO_NAME] [--exon-name EXON_NAME]
                 [--intron-name INTRON_NAME] [--search-file SEARCH_FILE]
                 [--coverage] [--keep] [--noFIFO] [--no-complexity]
                 [--prioritize] [-V]

PEPPRO version 0.9.7

optional arguments:
  -h, --help            show this help message and exit
  -R, --recover         Overwrite locks to recover from previous failed run
  -N, --new-start       Overwrite all results to start a fresh run
  -D, --dirty           Don't auto-delete intermediate files
  -F, --force-follow    Always run 'follow' commands
  -T, --testmode        Only print commands, don't run
  --silent              Silence logging. Overrides verbosity.
  --verbosity V         Set logging level (1-5 or logging module level name)
  --logdev              Expand content of logging message format.
  -C CONFIG_FILE, --config CONFIG_FILE
                        Pipeline configuration file (YAML). Relative paths are
                        with respect to the pipeline script.
  -M MEMORY_LIMIT, --mem MEMORY_LIMIT
                        Memory limit for processes accepting such. Default
                        units are megabytes unless specified using the suffix
                        [K|M|G|T].
  -P NUMBER_OF_CORES, --cores NUMBER_OF_CORES
                        Number of cores for parallelized processes
  -I2 [INPUT_FILES2 [INPUT_FILES2 ...]], --input2 [INPUT_FILES2 [INPUT_FILES2 ...]]
                        Secondary input files, such as read2
  -Q SINGLE_OR_PAIRED, --single-or-paired SINGLE_OR_PAIRED
                        Single- or paired-end sequencing protocol
  --protocol {PRO,pro,PRO-SEQ,PRO-seq,proseq,PROSEQ,GRO,gro,groseq,GROSEQ,GRO-SEQ,GRO-seq}
                        Run on sequencing type.
  --adapter-tool {cutadapt,fastp}
                        Name of adapter removal program. Default: cutadapt
  --dedup-tool {seqkit,fqdedup}
                        Program to use to duplicate reads. Default: seqkit
  --trimmer-tool {seqtk,fastx}
                        Name of read trimming program. Default: seqtk
  --umi-len UMI_LEN     Specify the length of the UMI.If your data does not
                        utilize UMIs, set to 0. Default: 0
  --max-len MAX_LEN     Trim reads to maximum length. Set to -1 to disable
                        length trimming. Default: -1
  --sob                 Use seqOutBias to produce signal tracks and
                        incorporate mappability information.
  --scale               Scale signal tracks: Default is to scale by read
                        count. If using seqOutBias, scales by the
                        expected/observed cut frequency.
  --prealignments PREALIGNMENTS [PREALIGNMENTS ...]
                        Space-delimited list of reference genomes to align to
                        before primary alignment.
  --TSS-name TSS_NAME   file_name of TSS annotation file.
  --pi-tss ENSEMBL_TSS  file_name of pause index TSS annotation file.
  --pi-body ENSEMBL_GENE_BODY
                        file_name of pause index gene body annotation file.
  --pre-name PRE_NAME   file_name of pre-mRNA annotation file.
  --anno-name ANNO_NAME
                        file_name of genomic annotation file.
  --exon-name EXON_NAME
                        file_name of exon annotation file.
  --intron-name INTRON_NAME
                        file_name of intron annotation file.
  --search-file SEARCH_FILE
                        file_name of read length matched gt tallymer index
                        search file
  --coverage            Report library complexity using coverage: reads /
                        (bases in genome / read length)
  --keep                Keep prealignment BAM files
  --noFIFO              Do NOT use named pipes during prealignments.
  --no-complexity       Disable library complexity calculation (faster).
  --prioritize          Plot cFRiF/FRiF using mutually exclusive priority
                        ranked features based on the order of feature
                        appearance in the feature annotation asset.
  -V, --version         show program's version number and exit

required named arguments:
  -O PARENT_OUTPUT_FOLDER, --output-parent PARENT_OUTPUT_FOLDER
                        Parent output directory of project
  -S SAMPLE_NAME, --sample-name SAMPLE_NAME
                        Name for sample to run
  -I INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        One or more primary input files
  -G GENOME_ASSEMBLY, --genome GENOME_ASSEMBLY
                        Identifier for genome assembly
```
