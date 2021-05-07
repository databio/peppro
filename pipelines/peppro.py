#!/usr/bin/env python3
"""
PEPPRO - Run-on sequencing pipeline
"""

__author__ = ["Jason Smith", "Nathan Sheffield", "Mike Guertin"]
__email__ = "jasonsmith@virginia.edu"
__version__ = "0.9.11-pipestat"

from argparse import ArgumentParser
import os
import sys
import re
import tempfile
import tarfile
import pypiper
import errno
from pypiper import build_command
from refgenconf import RefGenConf as RGC, select_genome_config
from pipestat.exceptions import PipestatDatabaseError

TOOLS_FOLDER = "tools"
RUNON_SOURCE_PRO = ["PRO", "pro", "PRO-SEQ", "PRO-seq", "proseq", "PROSEQ"]
RUNON_SOURCE_GRO = ["GRO", "gro", "groseq", "GROSEQ", "GRO-SEQ", "GRO-seq"]
RUNON_SOURCE = RUNON_SOURCE_PRO + RUNON_SOURCE_GRO

ADAPTER_REMOVERS = ["cutadapt", "fastp"]
DEDUPLICATORS = ["seqkit", "fqdedup"]
TRIMMERS = ["seqtk", "fastx"]

DEFAULT_REMOVER = "cutadapt"
DEFAULT_DEDUPLICATOR = "seqkit"
DEFAULT_TRIMMER = "seqtk"

BT2_IDX_KEY = "bowtie2_index"
DEFAULT_UMI_LEN = 0
DEFAULT_MAX_LEN = -1

def safely_retrieve(pipeline_manager, result_id, default=None):
    """
    Safely retrieve a value for a result. 
    
    Fall back to a desired default value in case of failure.

    :param pypiper.PipelineManager pipeline_manager: pm with .pipestat 
        property to retrieve the result from
    :param str result_id: result identifier to retrieve
    :param Any defult: a desired default value
    :return Any: retireved value
    """
    try:
        return pm.pipestat.retrieve(result_identifier=result_id)
    except PipestatDatabaseError:
        return default

def parse_arguments():
    """
    Parse command-line arguments passed to the pipeline.
    """
    # Argument Parsing from yaml file
    ###########################################################################
    parser = ArgumentParser(description='PEPPRO version ' + __version__)
    parser = pypiper.add_pypiper_args(parser, groups=
        ['pypiper', 'looper', 'ngs', 'pipestat'],
        required=["input", "genome", "sample-name", "output-parent"])

    # Pipeline-specific arguments
    parser.add_argument("--protocol", dest="protocol",
                        default="pro", choices=RUNON_SOURCE,
                        help="Run on sequencing type.")

    parser.add_argument("--adapter-tool", dest="adapter",
                        default=DEFAULT_REMOVER, choices=ADAPTER_REMOVERS,
                        help="Name of adapter removal program. "
                             "Default: {}".format(DEFAULT_REMOVER))

    parser.add_argument("--dedup-tool", dest="dedup",
                        default=DEFAULT_DEDUPLICATOR, choices=DEDUPLICATORS,
                        help="Program to use to duplicate reads. "
                             "Default: {}".format(DEFAULT_DEDUPLICATOR))

    parser.add_argument("--trimmer-tool", dest="trimmer",
                        default=DEFAULT_TRIMMER, choices=TRIMMERS,
                        help="Name of read trimming program. "
                             "Default: {}".format(DEFAULT_TRIMMER))

    parser.add_argument("--umi-len", 
                        default=DEFAULT_UMI_LEN, type=int,
                        help="Specify the length of the UMI."
                             "If your data does not utilize UMIs, set to 0. "
                             "Default: {}".format(DEFAULT_UMI_LEN))

    parser.add_argument("--max-len",
                        default=DEFAULT_MAX_LEN,
                        help="Trim reads to maximum length. "
                             "Set to -1 to disable length trimming. "
                             "Default: {}".format(DEFAULT_MAX_LEN))

    parser.add_argument("--sob", action='store_true',
                        dest="sob", default=False,
                        help="Use seqOutBias to produce signal tracks and "
                             "incorporate mappability information.")

    parser.add_argument("--scale", action='store_true',
                        dest="scale", default=False,
                        help="Scale signal tracks: "
                             "Default is to scale by read count.\n"
                             "If using seqOutBias, scales by the expected/"
                             "observed cut frequency.")

    parser.add_argument("--prealignments", default=[], type=str, nargs="+",
                        help="Space-delimited list of reference genomes to "
                             "align to before primary alignment.")

    parser.add_argument("--TSS-name", default=None,
                        dest="TSS_name", type=str,
                        help="file_name of TSS annotation file.")

    parser.add_argument("--pi-tss", default=None,
                        dest="ensembl_tss", type=str,
                        help="file_name of pause index TSS annotation file.")

    parser.add_argument("--pi-body", default=None,
                        dest="ensembl_gene_body", type=str,
                        help="file_name of pause index gene body annotation file.")

    parser.add_argument("--pre-name", default=None,
                        dest="pre_name", type=str,
                        help="file_name of pre-mRNA annotation file.")

    parser.add_argument("--anno-name", default=None,
                        dest="anno_name", type=str,
                        help="file_name of genomic annotation file.")

    parser.add_argument("--exon-name", default=None,
                        dest="exon_name", type=str,
                        help="file_name of exon annotation file.")

    parser.add_argument("--intron-name", default=None,
                        dest="intron_name", type=str,
                        help="file_name of intron annotation file.")

    parser.add_argument("--search-file", default=None,
                        dest="search_file", type=str,
                        help="file_name of read length matched gt tallymer "
                             "index search file")

    parser.add_argument("--coverage", action='store_true', default=False,
                        dest="coverage",
                        help="Report library complexity using coverage: "
                             "reads / (bases in genome / read length)")

    parser.add_argument("--keep", action='store_true', default=False,
                        dest="keep",
                        help="Keep prealignment BAM files")

    parser.add_argument("--noFIFO", action='store_true', default=False,
                        dest="no_fifo",
                        help="Do NOT use named pipes during prealignments.")

    parser.add_argument("--no-complexity", action='store_true', default=False,
                        dest="complexity",
                        help="Disable library complexity calculation (faster).")

    parser.add_argument("--prioritize", action='store_true', default=False,
                        dest="prioritize",
                        help="Plot cFRiF/FRiF using mutually exclusive priority"
                             " ranked features based on the order of feature"
                             " appearance in the feature annotation asset.")

    parser.add_argument("-V", "--version", action="version",
                        version="%(prog)s {v}".format(v=__version__))

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        raise SystemExit

    return args


def _remove_adapters(args, res, tools, read2, fq_file, outfolder):
    """
    A helper function to build a command for adapter removal.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict res: binding between resources and
        value, e.g. for resources used by the pipeline
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools used by the pipeline
    :param bool read2: if True, do not deduplicate and do not retain
        intermediate files
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to remove adapters
    """

    sname = args.sample_name  # for concise code

    cutadapt_folder = os.path.join(outfolder, "cutadapt")
    fastp_folder = os.path.join(outfolder, "fastp")
    fastq_folder = os.path.join(outfolder, "fastq")

    if read2:
        cutadapt_report = os.path.join(cutadapt_folder,
            sname + "_R2_cutadapt.txt")
        noadap_fastq = os.path.join(fastq_folder, sname + "_R2_noadap.fastq")
        short_fastq = os.path.join(fastq_folder, sname + "_R2_short.fastq")
        fastp_pfx = os.path.join(fastp_folder, sname + "_R2_fastp_adapter")
    else:
        cutadapt_report = os.path.join(cutadapt_folder,
            sname + "_R1_cutadapt.txt")
        noadap_fastq = os.path.join(fastq_folder, sname + "_R1_noadap.fastq")
        short_fastq = os.path.join(fastq_folder, sname + "_R1_short.fastq")
        fastp_pfx = os.path.join(fastp_folder, sname + "_R1_fastp_adapter")

    fastp_report_txt = fastp_pfx + ".txt"
    fastp_report_html = fastp_pfx + ".html"
    fastp_report_json = fastp_pfx + ".json"

    if _itsa_file(res.adapters):
        pm.info("Using custom adapter file: {}".format(res.adapters))
        five_prime = pm.checkprint("awk '/5prime/{getline; print}' " +
            res.adapters) or "TGGAATTCTCGGGTGCCAAGG"
        three_prime = pm.checkprint("awk '/3prime/{getline; print}' " +
            res.adapters) or "GATCGTCGGACTGTAGAACTCTGAAC"
    else:
        # Default to the hardcoded values as a fallback
        five_prime = "TGGAATTCTCGGGTGCCAAGG"
        three_prime = "GATCGTCGGACTGTAGAACTCTGAAC"

    # Setup report output folders
    if args.adapter == "cutadapt":
        ngstk.make_dir(cutadapt_folder)
        adapter_report = cutadapt_report
    elif args.adapter == "fastp":
        ngstk.make_dir(fastp_folder)
        adapter_report = fastp_report_txt

    # Create adapter trimming command(s).
    if args.adapter == "fastp":
        adapter_cmd_chunks = [
            ("(" + tools.fastp),
            ("--overrepresentation_analysis"),
            ("--thread", str(pm.cores)),
            ("--in1", fq_file)
        ]
        if read2:
            adapter_cmd_chunks.extend([
                ("--adapter_sequence", three_prime)
            ])
        else:
            adapter_cmd_chunks.extend([
                ("--adapter_sequence", five_prime)
            ])

        adapter_cmd_chunks.extend([
            ("--length_required", (2 + int(float(args.umi_len)))),
            ("--html", fastp_report_html),
            ("--json", fastp_report_json),
            ("--report_title", ("'" + sname + "'"))
        ])
        # If calculating library complexity and this is read 1 or single-end,
        # must produce an intermediate file.
        adapter_cmd_chunks.extend([("-o", noadap_fastq)])
        # Must keep intermediates always now

        adapter_cmd_chunks.extend([
            (") 2>", fastp_report_txt)
        ])

        adapter_cmd = build_command(adapter_cmd_chunks)

    elif args.adapter == "cutadapt":
        # Must keep intermediates always now
        cut_version = float(pm.checkprint("cutadapt --version"))
        adapter_cmd_chunks = ["(" + tools.cutadapt]
        # old versions of cutadapt can not use multiple cores
        if cut_version >= 1.15:
            adapter_cmd_chunks.extend([("-j", str(pm.cores))])
        adapter_cmd_chunks.extend([
            ("-m", (2 + int(float(args.umi_len)))),
            ("-O", 1)
        ])
        if read2:
            adapter_cmd_chunks.extend([("-a", three_prime)])
        else:
            adapter_cmd_chunks.extend([("-a", five_prime)])
        adapter_cmd_chunks.extend([
            fq_file,
            ("-o", noadap_fastq),
            ("--too-short-output", short_fastq),
            ")",
            (">", cutadapt_report)
        ])

        adapter_cmd = build_command(adapter_cmd_chunks)

    else:
        # Default to fastp
        adapter_cmd_chunks = [
            ("(" + tools.fastp),
            ("--overrepresentation_analysis"),
            ("--thread", str(pm.cores)),
            ("--in1", fq_file)
        ]
        if read2:
            adapter_cmd_chunks.extend([
                ("--adapter_sequence", three_prime)
            ])
        else:
            adapter_cmd_chunks.extend([
                ("--adapter_sequence", five_prime)
            ])

        adapter_cmd_chunks.extend([
            ("--length_required", (2 + int(float(args.umi_len)))),
            ("--html", fastp_report_html),
            ("--json", fastp_report_json),
            ("--report_title", ("'" + sname + "'"))
        ])
        # If calculating library complexity and this is read 1 or single-end,
        # must produce an intermediate file.
        adapter_cmd_chunks.extend([("-o", noadap_fastq)])
        # Must keep intermediates always now

        adapter_cmd_chunks.extend([
            (") 2>", fastp_report_txt)
        ])

        adapter_cmd = build_command(adapter_cmd_chunks)

    return adapter_cmd


def _deduplicate(args, tools, fq_file, outfolder):
    """
    A helper function to build a command for deduplication.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to remove adapters
    """
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    noadap_fastq = os.path.join(fastq_folder, sname + "_R1_noadap.fastq")
    dedup_fastq = os.path.join(fastq_folder, sname + "_R1_dedup.fastq")

    # Create deduplication command(s).
    if not args.complexity and int(args.umi_len) > 0:
        if args.dedup == "seqkit":
            dedup_cmd_chunks = [
                (tools.seqkit, "rmdup"),
                ("--threads", str(pm.cores)),
                "--by-seq",
                "--ignore-case",
                "-o"           
            ]
            if not args.complexity:
                dedup_cmd_chunks.extend([
                    (dedup_fastq, noadap_fastq)
                ])
            else:
                dedup_cmd_chunks.extend(["-"])
            dedup_cmd = build_command(dedup_cmd_chunks)
        elif args.dedup == "fqdedup":
            dedup_cmd_chunks = [tools.fqdedup]
            dedup_cmd_chunks.extend([("-i", noadap_fastq)])
            dedup_cmd_chunks.extend([("-o", dedup_fastq)])
            dedup_cmd = build_command(dedup_cmd_chunks)
        else:
            # Default to seqkit
            dedup_cmd_chunks = [
                (tools.seqkit, "rmdup"),
                ("--threads", str(pm.cores)),
                "--by-seq",
                "--ignore-case",
                "-o"           
            ]
            dedup_cmd_chunks.extend([
                (dedup_fastq, noadap_fastq)
            ])
            dedup_cmd = build_command(dedup_cmd_chunks)
    else:
        # Don't deduplicate a read2 file nor deduplicate if there are no UMI's
        dedup_cmd = ""

    return dedup_cmd


def _trim_deduplicated_files(args, tools, fq_file, outfolder):
    """
    A helper function to build a command for read trimming using fastq files
    that have been deduplicated.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to trim adapter trimmed and deduplicated reads
    """

    # Only call this when args.complexity32 and int(args.umi_len) > 0
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    dedup_fastq = os.path.join(fastq_folder, sname + "_R1_dedup.fastq")
    processed_fastq = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")

    fastp_folder = os.path.join(outfolder, "fastp")    
    umi_report = os.path.join(fastp_folder, sname + "_R1_rmUmi.html")
    umi_json = os.path.join(fastp_folder, sname + "_R1_rmUmi.json")

    # Check quality encoding for use with FastX_Tools
    if args.trimmer == "fastx":
        encoding = _guess_encoding(fq_file)

    if args.adapter == "fastp":
        # Remove UMI by specifying location of UMI
        # Location is still read1 because it's being treated as SE data
        trim_cmd_chunks = [
            tools.fastp,
            ("--thread", str(pm.cores)),
            ("-i", dedup_fastq),
            "--stdout",
            "--umi",
            ("--umi_loc", "read1"),
            ("--umi_len", args.umi_len),
            ("--html", umi_report),
            ("--json", umi_json)
        ]

        # Trim to max length if specified
        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "trimfq"),
                ("-L", args.max_len),
                "-"
            ])

        # Remove too short reads
        trim_cmd_chunks.extend([
            "|",
            (tools.seqtk, "seq"),
            #("-L", 5)
            ("-L",(2 + int(float(args.umi_len))))
        ])

        # Reverse complement for PRO-seq
        if args.protocol.lower() in RUNON_SOURCE_PRO:
            trim_cmd_chunks.extend([("-r", "-")])
        else:
            trim_cmd_chunks.extend(["-"])

        trim_cmd_chunks.extend([(">", processed_fastq)])

    elif args.trimmer == "seqtk":
        # Remove UMI by blind trimming
        trim_cmd_chunks = [
            tools.seqtk,
            "trimfq",
            ("-b", str(args.umi_len))
        ]

        # Trim to max length if specified
        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([("-L", str(args.max_len))])

        trim_cmd_chunks.extend([dedup_fastq])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                "-",
                (">", processed_fastq)
            ])                           
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                ("-r", "-"),
                (">", processed_fastq)
            ])
    elif args.trimmer == "fastx":
        trim_tool = tools.fastx + "_trimmer"
        rc_tool = tools.fastx + "_reverse_complement"
        trim_cmd_chunks = [trim_tool]

        if encoding == "Illumina-1.8":
            trim_cmd_chunks.extend([
                ("-Q", str(33))
            ])

        # Remove UMI blindly
        trim_cmd_chunks.extend([("-f", str(int(float(args.umi_len)) + 1))])

        # Trim to max length if specified
        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([
                ("-l", (str(int(float(args.max_len)) +
                 int(float(args.umi_len)))))
            ])

        # Remove too short reads
        trim_cmd_chunks.extend([("-m", (2 + int(float(args.umi_len))))])

        # Add input file
        trim_cmd_chunks.extend([("-i", dedup_fastq)])

        # Do not reverse complement if GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([("-o", processed_fastq)])                        
        else:
            trim_cmd_chunks.extend([("|", rc_tool)])
            if encoding == "Illumina-1.8":
                trim_cmd_chunks.extend([
                    ("-Q", str(33))
                ])
            trim_cmd_chunks.extend([("-o", processed_fastq)])
    else:
        # Default to seqtk
        trim_cmd_chunks = [
            tools.seqtk,
            "trimfq",
            ("-b", str(args.umi_len))
        ]

        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([("-L", str(args.max_len))])

        trim_cmd_chunks.extend([dedup_fastq])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                "-",
                (">", processed_fastq)
            ])                           
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                ("-r", "-"),
                (">", processed_fastq)
            ])

    trim_cmd = build_command(trim_cmd_chunks)
    pm.debug("trim_deduplicated_cmd: {}".format(build_command(trim_cmd_chunks)))

    return trim_cmd


def _trim_adapter_files(args, tools, read2, fq_file, outfolder):
    """
    A helper function to build a command for read trimming using fastq files
    without deduplication.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param bool read2: if True, do not deduplicate and do not retain
        intermediate files
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to trim adapter trimmed files
    """
    # Need undeduplicated results for complexity calculation
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    if read2:
        noadap_fastq = os.path.join(fastq_folder, sname + "_R2_noadap.fastq")
        trimmed_fastq = os.path.join(fastq_folder, sname + "_R2_trimmed.fastq")
    else:
        noadap_fastq = os.path.join(fastq_folder,
            sname + "_R1_noadap.fastq")
        trimmed_fastq = os.path.join(fastq_folder,
            sname + "_R1_processed.fastq")

    fastp_folder = os.path.join(outfolder, "fastp")
    if read2:
        umi_report = os.path.join(fastp_folder, sname + "_R2_rmUmi.html")
        umi_json = os.path.join(fastp_folder, sname + "_R2_rmUmi.json")
    else:
        umi_report = os.path.join(fastp_folder, sname + "_R1_rmUmi.html")
        umi_json = os.path.join(fastp_folder, sname + "_R1_rmUmi.json")

    # Check quality encoding for use with FastX_Tools
    if args.trimmer == "fastx":
        encoding = _guess_encoding(fq_file)
    
    if args.adapter == "fastp":
        # Remove UMI and specify location of UMI
        # Still requires seqtk for reverse complementation
        if int(args.umi_len) > 0:
            trim_cmd_chunks = [
                tools.fastp,
                ("--thread", str(pm.cores)),
                ("-i", noadap_fastq),
                "--stdout",
                "--umi",
                ("--umi_loc", "read1"),
                ("--umi_len", args.umi_len),
                ("--html", umi_report),
                ("--json", umi_json)
            ]

            if int(args.max_len) > 0:
                # Trim to max length if specified
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "trimfq"),
                    ("-L", args.max_len),
                    "-"
                ])

        elif int(args.max_len) > 0:
            trim_cmd_chunks = [
                (tools.seqtk, "trimfq"),
                ("-L", args.max_len),
                noadap_fastq
            ]
        else:
            trim_cmd_chunks = []

        if trim_cmd_chunks:
            # Reverse complement for PRO-seq
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5)
                ("-L",(2 + int(float(args.umi_len))))
            ])
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                trim_cmd_chunks.extend([("-r", "-")])
            else:
                trim_cmd_chunks.extend(["-"])

            trim_cmd_chunks.extend([(">", trimmed_fastq)])

        else:
            # If no UMI removal or read trimming, just reverse complement
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                trim_cmd_chunks.extend([
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    ("-r", noadap_fastq),
                    (">", trimmed_fastq)
                ])
            # Otherwise just make sure to remove too short reads
            else:
                trim_cmd_chunks = [
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    noadap_fastq,
                    (">", trimmed_fastq)
                ]
    elif args.trimmer == "seqtk":
        # Remove UMI blindly by position
        trim_cmd_chunks = [
            tools.seqtk,
            "trimfq",
            ("-b", str(args.umi_len))
        ]

        # Trim tp max length if specified
        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([
                ("-L", str(args.max_len))
            ])
        
        trim_cmd_chunks.extend([noadap_fastq])
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                "-",
                (">", trimmed_fastq)
            ])
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                ("-r", "-"),
                (">", trimmed_fastq)
            ])
    elif args.trimmer == "fastx":
        trim_tool = tools.fastx + "_trimmer"
        rc_tool = tools.fastx + "_reverse_complement"
        trim_cmd_chunks = [trim_tool]
        if encoding == "Illumina-1.8":
            trim_cmd_chunks.extend([
                ("-Q", str(33))
            ])

        # Remove UMI blindly by position only
        trim_cmd_chunks.extend([
            ("-f", str(int(float(args.umi_len)) + 1))
        ])

        # Trim tp max length if specified
        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([
                ("-l", (str(int(float(args.max_len)) +
                 int(float(args.umi_len)))))
            ])

        # Remove too short reads
        trim_cmd_chunks.extend([("-m", (2 + int(float(args.umi_len))))])

        # Need undeduplicated results for complexity calculation
        trim_cmd_chunks.extend([
            ("-i", noadap_fastq)
        ])
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                ("-o", trimmed_fastq)
            ])
        else:
            trim_cmd_chunks.extend([("|", rc_tool)])
            if encoding == "Illumina-1.8":
                trim_cmd_chunks.extend([
                    ("-Q", str(33))
                ])
            trim_cmd_chunks.extend([("-o", trimmed_fastq)])
    else:
        # Default to seqtk
        # Remove UMI blindly by position only
        trim_cmd_chunks = [
            tools.seqtk,
            "trimfq",
            ("-b", str(args.umi_len))
        ]

        # Trim to max length if specified
        if int(args.max_len) > 0:
            trim_cmd_chunks.extend([
                ("-L", str(args.max_len))
            ])

        trim_cmd_chunks.extend([noadap_fastq])
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            # Do not reverse complement
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                "-",
                (">", trimmed_fastq)
            ])
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                ("-r", "-"),
                (">", trimmed_fastq)
            ])

    trim_cmd = build_command(trim_cmd_chunks)
    pm.debug("trim_cmd_nodedup: {}".format(build_command(trim_cmd_chunks)))

    return trim_cmd


def _trim_pipes(args, tools, read2, fq_file, outfolder):
    """
    A helper function to build a command for read trimming using pipes.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param bool read2: if True, do not deduplicate and do not retain
        intermediate files
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to trim adapter trimmed and deduplicated reads
    """

    # Only call this when NOT args.complexity or NOT int(args.umi_len) > 0
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    if read2:
        noadap_fastq = os.path.join(fastq_folder,
            sname + "_R2_noadap.fastq")
        processed_fastq = os.path.join(fastq_folder,
            sname + "_R2_trimmed.fastq")
    else:
        noadap_fastq = os.path.join(fastq_folder,
            sname + "_R1_noadap.fastq")
        processed_fastq = os.path.join(fastq_folder,
            sname + "_R1_processed.fastq")

    fastp_folder = os.path.join(outfolder, "fastp")
    if read2:
        umi_report = os.path.join(fastp_folder, sname + "_R2_rmUmi.html")
        umi_json = os.path.join(fastp_folder, sname + "_R2_rmUmi.json")
    else:
        umi_report = os.path.join(fastp_folder, sname + "_R1_rmUmi.html")
        umi_json = os.path.join(fastp_folder, sname + "_R1_rmUmi.json")

    # Check quality encoding for use with FastX_Tools
    if args.trimmer == "fastx":
        encoding = _guess_encoding(fq_file)

    if args.adapter == "fastp":
        # There are no intermediate files, just pipes
        # Remove UMI
        if int(args.umi_len) > 0:
            trim_cmd_chunks = [
                tools.fastp,
                ("--thread", str(pm.cores)),
                ("--stdin", "--stdout"),
                "--umi",
                ("--umi_loc", "read1"),
                ("--umi_len", args.umi_len),
                ("--html", umi_report),
                ("--json", umi_json),
                "|",
                (tools.seqtk, "trimfq")
            ]

            # Trim to max length if specified
            if int(args.max_len) > 0:
                trim_cmd_chunks.extend([("-L", args.max_len)])

            trim_cmd_chunks.extend(["-"])

            # Do not reverse complement for GRO-seq
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    "-",
                    (">", processed_fastq)
                ])
            else:
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    ("-r", "-"),
                    (">", processed_fastq)
                ])
        elif int(args.max_len) > 0:
            # No UMI, but still trim max length
            trim_cmd_chunks = [
                (tools.seqtk, "trimfq"),
                ("-L", args.max_len),
                "-"
            ]
            # Do not reverse complement for GRO-seq
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    "-",
                    (">", processed_fastq)
                ])
            else:
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    ("-r", "-"),
                    (">", processed_fastq)
                ])
        else:
            # No UMI and no trimming
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                trim_cmd_chunks = [
                    (tools.seqtk, "seq"),
                    #("-L", 5),
                    ("-L",(2 + int(float(args.umi_len)))),
                    ("-r", noadap_fastq),
                    (">", processed_fastq)
                ]
            else:
                trim_cmd_chunks = []
    # if not args.complexity and int(args.umi_len) > 0 retain intermediate files
    elif args.trimmer == "seqtk":
        trim_cmd_chunks = [
            tools.seqtk,
            "trimfq"
        ]

        if read2:
            trim_cmd_chunks.extend([("-e", str(args.umi_len))])
        else:
            trim_cmd_chunks.extend([("-b", str(args.umi_len))])
            if int(args.max_len) > 0:
                trim_cmd_chunks.extend([("-L", str(args.max_len))])

        trim_cmd_chunks.extend(["-"])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                "-",
                (">", processed_fastq)
            ])
        else:
            trim_cmd_chunks.extend([
                "|",
                tools.seqtk,
                ("seq", "-r"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                ("-", ">"),
                processed_fastq
            ])
    elif args.trimmer == "fastx":
        trim_tool = tools.fastx + "_trimmer"
        rc_tool = tools.fastx + "_reverse_complement"
        trim_cmd_chunks = [trim_tool]

        if encoding == "Illumina-1.8":
            trim_cmd_chunks.extend([("-Q", str(33))])

        if read2:
            trim_cmd_chunks.extend([("-t", str(int(float(args.umi_len))))])
        else:
            trim_cmd_chunks.extend([("-f", str(int(float(args.umi_len)) + 1))])
            if int(args.max_len) > 0:
                trim_cmd_chunks.extend([
                    ("-l", (str(int(float(args.max_len)) +
                     int(float(args.umi_len)))))
                ])

        # Remove too short reads
        trim_cmd_chunks.extend([("-m", (2 + int(float(args.umi_len))))])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([("-o", processed_fastq)])                        
        else:
            trim_cmd_chunks.extend([("|", rc_tool)])
            if encoding == "Illumina-1.8":
                trim_cmd_chunks.extend([("-Q", str(33))])
            trim_cmd_chunks.extend([("-o", processed_fastq)])
    else:
        # Default to seqtk
        trim_cmd_chunks = [
            tools.seqtk,
            "trimfq"
        ]

        if read2:
            trim_cmd_chunks.extend([("-e", str(args.umi_len))])
        else:
            trim_cmd_chunks.extend([("-b", str(args.umi_len))])
            if int(args.max_len) > 0:
                trim_cmd_chunks.extend([("-L", str(args.max_len))])

        trim_cmd_chunks.extend(["-"])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                "-",
                (">", processed_fastq)
            ])
        else:
            trim_cmd_chunks.extend([
                "|",
                tools.seqtk,
                ("seq", "-r"),
                #("-L", 5),
                ("-L",(2 + int(float(args.umi_len)))),
                ("-", ">"),
                processed_fastq
            ])

    trim_cmd = build_command(trim_cmd_chunks)
    pm.debug("trim_pipes_cmd: {}".format(build_command(trim_cmd_chunks)))
    pm.debug("trim_pipes_cmd read2 status: {}".format(read2))

    return trim_cmd


def _process_fastq(args, tools, res, read2, fq_file, outfolder):
    """
    A helper function to prepare read files for downstream processing.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict res: binding between resources and
        value, e.g. for resources used by the pipeline
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools used by the pipeline
    :param bool read2: if True, do not deduplicate and do not retain
        intermediate files
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return (str, str): pair (R1, R2) of paths to FASTQ files
    """

    sname = args.sample_name  # for concise code

    # Create names for processed FASTQ files.
    fastq_folder = os.path.join(outfolder, "fastq")
    fastp_folder = os.path.join(outfolder, "fastp")
    fastqc_folder = os.path.join(outfolder, "fastqc")
    fastqc_report = os.path.join(fastqc_folder,
        sname + "_R1_processed_fastqc.html")

    preprocessed_fq1 = os.path.join(fastq_folder, sname + "_R1.fastq")
    preprocessed_fq2 = os.path.join(fastq_folder, sname + "_R2.fastq")
    noadap_fq1 = os.path.join(fastq_folder, sname + "_R1_noadap.fastq")
    noadap_fq2 = os.path.join(fastq_folder, sname + "_R2_noadap.fastq")
    short_fq1 = os.path.join(fastq_folder, sname + "_R1_short.fastq")
    short_fq2 = os.path.join(fastq_folder, sname + "_R2_short.fastq")
    dedup_fq = os.path.join(fastq_folder, sname + "_R1_dedup.fastq")
    trimmed_fq1 = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")
    trimmed_fq2 = os.path.join(fastq_folder, sname + "_R2_trimmed.fastq")
    trimmed_dups_fq2 = os.path.join(fastq_folder,
        sname + "_R2_trimmed_dups.fastq")
    processed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")

    if args.adapter == "cutadapt":
        cutadapt_folder = os.path.join(outfolder, "cutadapt")
        if read2:
            cutadapt_report = os.path.join(cutadapt_folder,
                                         sname + "_R2_cutadapt.txt")
        else:
            cutadapt_report = os.path.join(cutadapt_folder,
                                           sname + "_R1_cutadapt.txt")
        adapter_report = cutadapt_report
    else:
        adapter_report = os.path.join(fastqc_folder,
                                      sname + "_R1_rmAdapter.txt")

    fastp_pfx = os.path.join(fastp_folder, sname + "_R1_fastp_adapter")
    fastp_report_txt = fastp_pfx + ".txt"
    fastp_report_html = fastp_pfx + ".html"

    adapter_command = _remove_adapters(args, res, tools, read2, fq_file, outfolder)
    pm.debug("Adapter command: {}".format(adapter_command))
    pm.debug("Read2 status: {}".format(read2))

    # To plot fragment sizes requires keeping intermediate files
    if not args.complexity and int(args.umi_len) > 0:
        deduplicate_command = _deduplicate(args, tools, fq_file, outfolder)
        pm.debug("Dedup command: {}".format(deduplicate_command))
        trim_command = _trim_adapter_files(args, tools, read2, fq_file, outfolder)
        trim_command2 = _trim_deduplicated_files(args, tools, fq_file, outfolder)
    else:
        trim_command = _trim_adapter_files(args, tools, read2, fq_file, outfolder)

    def report_fastq():
        """
        Report QC metrics on intermediate steps of fastq file preparation
        """
        if args.adapter == "cutadapt":
            report = cutadapt_report
            adapter_term = "Reads with adapters:"
            #too_short_term = "Reads that were too short:"
            total_bases_term = "Total basepairs processed:"

            ac_cmd = ("grep '" + adapter_term + "' " +
                      report + " | awk '{print $(NF-1)}'")
            # cutadapt version < 2.9
            #ts_cmd = ("grep '" + too_short_term + "' " +
            #          report + " | awk '{print $(NF-1)}'")
            ts_cmd = "wc -l "
            if read2:
                ts_cmd += short_fq2
            else:
                ts_cmd += short_fq1
            ts_cmd += " | awk '{print $1}'"
            bases = ("grep '" + total_bases_term + "' " +
                     report + " | awk '{print $(NF-1)}'")
            adapter_bases = ("awk '{sum+=$1*$2} END {printf \"%.0f\", sum}' " +
                             report)

        else:  # default to fastp
            report = fastp_report_txt
            adapter_term = "reads with adapter trimmed:"
            too_short_term = "reads failed due to too short:"
            total_bases_term = "total bases:"

            ac_cmd = ("grep '" + adapter_term + "' " +
                       report + " | head -n 1 | awk '{print $NF}'")
            ts_cmd = ("grep '" + too_short_term + "' " +
                       report + " | head -n 1 | awk '{print $NF}'")
            bases = ("grep '" + total_bases_term + "' " +
                     report + " | head -n 1 | awk '{print $NF}'")
            adapter_bases = ("grep 'bases trimmed due to adapters:' " +
                             report + " | awk '{print $NF}'")

            # pm.report_object("FastP_report", fastp_report_html)
            pm.pipestat.report(values={"FastP_report": {"path": mRNApdf, "title": "FastP_report"}})

        if _itsa_file(report):
            tmp = pm.checkprint(ac_cmd)
            if tmp:
                ac = float(tmp.replace(',',''))
            else:
                ac = 0
            pm.pipestat.report(values={"Reads_with_adapter": ac})    
            # pm.report_result("Reads_with_adapter", ac)
            
            tmp = pm.checkprint(bases)
            if tmp:
                total_bases = float(tmp.replace(',',''))
            else:
                total_bases = 0
                
            tmp = pm.checkprint(adapter_bases)
            if tmp:
                total_adapter = float(tmp.replace(',',''))
            else:
                total_adapter = 0
            
            tmp = pm.checkprint(ts_cmd)
            if tmp:
                ts = float(tmp.replace(',',''))
                if args.adapter == "cutadapt":
                    ts = ts/4
            else:
                ts = 0
            
            pm.pipestat.report(values={"Uninformative_adapter_reads": round(ts, 2)})
            # pm.report_result("Uninformative_adapter_reads", round(ts, 2))

            if _itsa_file(noadap_fq1):
                tr = int(ngstk.count_lines(noadap_fq1).strip())
            else:
                tr = 0

            if _itsa_file(dedup_fq):
                dr = int(ngstk.count_lines(dedup_fq).strip())
                dups = max(0, (float(tr)/4 - float(dr)/4))
                pm.pipestat.report(values={"Duplicate_reads": round(dups, 2)})
                # pm.report_result("Duplicate_reads", round(dups, 2))

            if _itsa_file(preprocessed_fq1):
                pre = int(ngstk.count_lines(preprocessed_fq1).strip())
                pm.pipestat.report(values={"Pct_uninformative_adapter_reads": round(float(100*(ts/(float(pre)/4))), 4)})
                # pm.report_result("Pct_uninformative_adapter_reads", 
                    # round(float(100*(ts/(float(pre)/4))), 4))
        else:
            pm.fail_pipeline("Could not find '{}' to report adapter "
                             "removal statistics.".format(report))


    def plot_fragments(infolder, outfolder):
        """
        Plot adapter insertion distribution (from PE data only)

        :param str infolder: path to fastq containing directory
        :param str outfolder: path to output directory for functions
        """
        # merge short fragment reads
        noadap_fq1 = os.path.join(infolder,
            args.sample_name + "_R1_noadap.fastq")
        noadap_fq2 = os.path.join(infolder,
            args.sample_name + "_R2_noadap.fastq")
        rep_fq1 = os.path.join(infolder,
            args.sample_name + "_R1_noadap.fastq.paired.fq")
        rep_fq2 = os.path.join(infolder,
            args.sample_name + "_R2_noadap.fastq.paired.fq")
        orphan_fq1 = os.path.join(infolder,
            args.sample_name + "_R1_noadap.fastq.single.fq")
        orphan_fq2 = os.path.join(infolder,
            args.sample_name + "_R2_noadap.fastq.single.fq")
        flash_hist = os.path.join(outfolder, args.sample_name + ".hist")
        flash_gram = os.path.join(outfolder, args.sample_name + ".histogram")
        flash_extended = os.path.join(outfolder,
            args.sample_name + ".extendedFrags.fastq.gz")
        flash_notCombined_fq1 = os.path.join(outfolder,
            args.sample_name + ".notCombined_1.fastq.gz")
        flash_notCombined_fq2 = os.path.join(outfolder,
            args.sample_name + ".notCombined_2.fastq.gz")

        tmp = float(pm.pipestat.retrieve(None,"Raw_reads"))
        if tmp:
            rr = float(tmp)
        else:
            rr = 0
        if (rr < 1):
            pm.fail_pipeline(RuntimeError("Raw_reads were not reported. "
                "Check output ({})".format(param.outfolder)))

        cmd1 = (tools.fastqpair + " -t " + str(int(0.9*rr)) + " " + 
                noadap_fq1 + " " + noadap_fq2)
        ziptool = "pigz" if ngstk.check_command("pigz") else "gzip"
        cmd2 = (tools.flash + " -q -t " + str(pm.cores) +
                " --compress-prog=" + ziptool + " --suffix=gz " +
                rep_fq1 + " " + rep_fq2 + " -o " + args.sample_name +
                " -d " + outfolder)
        pm.run([cmd1, cmd2], [flash_hist, flash_gram])
        pm.clean_add(rep_fq1)
        pm.clean_add(rep_fq2)
        pm.clean_add(orphan_fq1)
        pm.clean_add(orphan_fq2)
        pm.clean_add(flash_extended)
        pm.clean_add(flash_notCombined_fq1)
        pm.clean_add(flash_notCombined_fq2)

        pm.timestamp("### Plot adapter insertion distribution")

        degradation_pdf = os.path.join(outfolder,
            args.sample_name + "_adapter_insertion_distribution.pdf")
        degradation_png = os.path.join(outfolder,
            args.sample_name + "_adapter_insertion_distribution.png")
        cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") + 
               " adapt -i " + flash_hist + " -o " + outfolder)
        if int(args.umi_len) > 0:
            cmd += (" -u " + str(args.umi_len))
            umi_len = args.umi_len
        else:
            umi_len = 0
        pm.run(cmd, degradation_pdf, nofail=True)
        # pm.report_object("Adapter insertion distribution", degradation_pdf,
                        #  anchor_image=degradation_png)
        pm.pipestat.report(values={"Adapter_insertion_distribution": {"path": degradation_pdf, "thumbnail_path": degradation_png, "title": "Adapter insertion distribution"}})

        if not safely_retrieve(pm, "Peak_adapter_insertion_size") or args.new_start:
            # Report the peak insertion size
            cmd = ("awk 'NR>2 {print prev} {prev=$0}' " + flash_hist + 
                   " | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max)" +
                   " {max=$2; len=$1}; max_len=$1} END{print len-" +
                   str(umi_len) + "}'")
            adapter_peak = pm.checkprint(cmd)
            if adapter_peak:
                ap = int(adapter_peak)
                # pm.report_result("Peak_adapter_insertion_size", ap)
                pm.pipestat.report(values={"Peak_adapter_insertion_size": ap})

        # Report the degradation ratio
        if not pm.pipestat.retrieve(None,'Degradation_ratio') or args.new_start:
            pm.timestamp("###  Calculating degradation ratio")

            cmd = ("awk '{ if (($1-" + str(umi_len) + ") == 10) {status = 1}} " + 
                   "END {if (status) {print status} else {print 0}}' " +
                   flash_hist)
            degraded_lower = pm.checkprint(cmd)
            cmd = ("awk '{ if (($1-" + str(umi_len) + ") == 20) {status = 1}} " + 
                   "END {if (status) {print status} else {print 0}}' " +
                   flash_hist)
            degraded_upper = pm.checkprint(cmd)
            cmd = ("awk '{ if (($1-" + str(umi_len) + ") == 30) {status = 1}} " + 
                   "END {if (status) {print status} else {print 0}}' " +
                   flash_hist)
            intact_lower = pm.checkprint(cmd)
            cmd = ("awk '{ if (($1-" + str(umi_len) + ") == 40) {status = 1}} " + 
                   "END {if (status) {print status} else {print 0}}' " +
                   flash_hist)
            intact_upper = pm.checkprint(cmd)

            if degraded_lower:
                dl = int(degraded_lower)
            if dl == 1:
                dl = 10
            else:
                cmd = ("awk 'NR==1 {print ($1-" + str(umi_len) + ")}' " + flash_hist)
                degraded_lower = pm.checkprint(cmd)
                dl = int(degraded_lower) if degraded_lower else 1
                if dl < 1:
                    dl = 1

            if degraded_upper:
                du = int(degraded_upper)
            if du == 1:
                du = 20
            else:
                du = int(degraded_lower) + 10

            if intact_upper:
                iu = int(intact_upper)
            if iu == 1:
                iu = 40
            else:
                cmd = ("awk 'END {print ($1-" + str(umi_len) + ")}' " + flash_hist)
                intact_upper = pm.checkprint(cmd)
                dl = int(intact_upper) if intact_upper else 40

            if intact_lower:
                il = int(intact_lower)
            if il == 1:
                il = 30
            else:
                il = int(intact_upper) - 10
                if il < 1:
                    il = 30

            cmd = ("awk '(($1-" + str(umi_len) + " ) <= " + str(du) +
                   " && ($1-" + str(umi_len) + " ) >= " + str(dl) +
                   "){degradedSum += $2}; " +  "(($1-" + str(umi_len) +
                   " ) >= " + str(il) + " && ($1-" + str(umi_len) +
                   ") <= " + str(iu) +
                   "){intactSum += $2}  END {if (intactSum < 1) " +
                   "{intactSum = 1} print degradedSum/intactSum}' "
                   + flash_hist)
            degradation_ratio = pm.checkprint(cmd)
            if degradation_ratio:
                dr = float(degradation_ratio)
                # pm.report_result("Degradation_ratio", round(dr, 4))
                pm.pipestat.report(values={"Degradation_ratio": round(dr, 4)})


    def check_trim(trimmed_fastq, paired_end,
                   trimmed_fastq_R2=None, fastqc_folder=None):
        """
        Evaluate read trimming, and optionally run fastqc.
        
        :param str trimmed_fastq: Path to trimmed reads file.
        :param bool paired_end: Whether the processing is being done with
            paired-end sequencing data.
        :param str trimmed_fastq_R2: Path to read 2 file for paired-end case.
        :param str fastqc_folder: Path to folder within which to place fastqc
            output files; if unspecified, fastqc will not be run.
        :return callable: Function to evaluate read trimming and possibly run
            fastqc.
        """
        n_trim = float(ngstk.count_reads(trimmed_fastq, paired_end))
        # pm.report_result("Trimmed_reads_R1", int(n_trim))
        pm.pipestat.report(values={"Trimmed_reads_R1": int(n_trim)})

        try:
            rr = float(pm.pipestat.retrieve(None,"Raw_reads"))
        except:
            print("Can't calculate trim loss rate without raw read result.")
        else:
            pm.pipestat.report(values={"Trim_loss_rate_R1": round((rr - n_trim) * 100 / rr, 2)})
            # pm.report_result("Trim_loss_rate_R1", round((rr - n_trim) * 100 / rr, 2))

        if paired_end and trimmed_fastq_R2:
            n_trim = float(ngstk.count_reads(trimmed_fastq_R2, paired_end))
            pm.pipestat.report(values={"Trimmed_reads_R2": int(n_trim)})
            # pm.report_result("Trimmed_reads_R2", int(n_trim))

            try:
                rr = float(pm.pipestat.retrieve(None,"Raw_reads"))
            except:
                print("Can't calculate trim loss rate without raw read result.")
            else:
                pm.pipestat.report(values={"Trim_loss_rate_R2":
                                 round((rr - n_trim) * 100 / rr, 2)})
                # pm.report_result("Trim_loss_rate_R2",
                                #  round((rr - n_trim) * 100 / rr, 2))

        # Also run a fastqc (if installed/requested)
        if fastqc_folder:
            if fastqc_folder and os.path.isabs(fastqc_folder):
                try:
                    os.makedirs(fastqc_folder)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
            cmd = ngstk.fastqc(trimmed_fastq, fastqc_folder)
            pm.run(cmd, lock_name="trimmed_fastqc", nofail=True)
            fname, ext = os.path.splitext(os.path.basename(trimmed_fastq))
            fastqc_html = os.path.join(fastqc_folder, fname + "_fastqc.html")
            # pm.report_object("FastQC report R1", fastqc_html)
            pm.pipestat.report(values={"FastQC_report_R1": {"path": fastqc_html, "title": "FastQC report R1"}})
            
            if paired_end and trimmed_fastq_R2:
                    cmd = ngstk.fastqc(trimmed_fastq_R2, fastqc_folder)
                    pm.run(cmd, lock_name="trimmed_fastqc_R2", nofail=True)
                    fname, ext = os.path.splitext(
                        os.path.basename(trimmed_fastq_R2))
                    fastqc_html = os.path.join(
                        fastqc_folder, fname + "_fastqc.html")
                    # pm.report_object("FastQC report R2", fastqc_html)
                    pm.pipestat.report(values={"FastQC report R2": {"path": fastqc_html, "title": "FastQC report R2"}})

    # Put it all together
    paired_end = args.paired_end
    if read2:
        pm.run([adapter_command, trim_command], trimmed_fq2)
        if not _itsa_file(fastqc_report) or args.new_start:
            cmd = ("echo '### Calculated the number of trimmed reads'")
            pm.run(cmd, fastqc_report, 
                   follow=check_trim(processed_fastq, paired_end, trimmed_fq2,
                                     fastqc_folder=fastqc_folder))
        if args.adapter == "cutadapt":
            output_folder = os.path.join(outfolder, "cutadapt")
        else:
            output_folder = os.path.join(outfolder, "fastp")
        cp_cmd = ("cp " + trimmed_fq2 + " " + trimmed_dups_fq2)
        pm.run(cp_cmd, trimmed_dups_fq2,
               follow=plot_fragments(fastq_folder, output_folder))
        pm.clean_add(short_fq2)
        return trimmed_fq2, trimmed_dups_fq2
    elif not args.complexity and int(args.umi_len) > 0:
        # This trim command DOES need the adapter file...
        pm.debug("\ntrim_command1: {} +\n {}\n".format(adapter_command, trim_command))
        pm.run([adapter_command, trim_command], processed_fastq)
        if not _itsa_file(fastqc_report) or args.new_start:
            cmd = ("echo '### Calculated the number of trimmed reads'")
            pm.run(cmd, fastqc_report, 
                   follow=check_trim(processed_fastq, paired_end, None,
                                     fastqc_folder=fastqc_folder))
        # This needs to produce the trimmed_fastq file
        pm.debug("\ntrim_command2: {} +\n {}\n".format(deduplicate_command, trim_command2))
        pm.run([deduplicate_command, trim_command2],
               trimmed_fq1, follow=report_fastq)
        pm.clean_add(noadap_fq1)
        pm.clean_add(short_fq1)
        pm.clean_add(dedup_fq)
        pm.clean_add(trimmed_fq1)
        return processed_fastq, trimmed_fq1
    else:
        pm.debug("\nELSE: trim_command: {} + {}\n".format(adapter_command, trim_command))
        pm.run([adapter_command, trim_command], processed_fastq,
               follow=report_fastq)
        if not _itsa_file(fastqc_report) or args.new_start:
            cmd = ("echo '### Calculate the number of trimmed reads'")
            pm.run(cmd, fastqc_report, 
                   follow=check_trim(processed_fastq, paired_end, None,
                                     fastqc_folder=fastqc_folder))
        pm.clean_add(noadap_fq1)
        pm.clean_add(short_fq1)
        return processed_fastq


def _align_with_bt2(args, tools, paired, useFIFO, unmap_fq1, unmap_fq2,
                    assembly_identifier, assembly_bt2, outfolder,
                    aligndir=None, bt2_opts_txt=None, dups=False):
    """
    A helper function to run alignments in series, so you can run one alignment
    followed by another; this is useful for successive decoy alignments.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param bool paired: if True, use paired-end alignment
    :param bool useFIFO: if True, use named pipe instead of file creation
    :param str unmap_fq1: path to unmapped read1 FASTQ file
    :param str unmap_fq2: path to unmapped read2 FASTQ file
    :param str assembly_identifier: text identifying a genome assembly for the
        pipeline
    :param str assembly_bt2: assembly-specific bowtie2 folder (index, etc.)
    :param str outfolder: path to output directory for the pipeline
    :param str aligndir: name of folder for temporary output
    :param str bt2_opts_txt: command-line text for bowtie2 options
    :param bool dups: if True, produce alternative named output
    :return (str, str): pair (R1, R2) of paths to FASTQ files
    """
    if os.path.exists(os.path.dirname(assembly_bt2)):
        pm.timestamp("### Map to " + assembly_identifier)
        if not aligndir:
            align_subdir = "aligned_{}_{}".format(args.genome_assembly,
                                                  assembly_identifier)
            sub_outdir = os.path.join(outfolder, align_subdir)
        else:
            sub_outdir = os.path.join(outfolder, aligndir)

        ngstk.make_dir(sub_outdir)
        if dups:
            bamname = "{}_{}_dups.bam".format(args.sample_name, assembly_identifier)
            summary_name = "{}_{}_bt_aln_dups_summary.log".format(args.sample_name,
                                                                  assembly_identifier)
        else:
            bamname = "{}_{}.bam".format(args.sample_name, assembly_identifier)
            summary_name = "{}_{}_bt_aln_summary.log".format(args.sample_name,
                                                             assembly_identifier)
        mapped_bam = os.path.join(sub_outdir, bamname)
        summary_file = os.path.join(sub_outdir, summary_name)

        out_fastq_pre = os.path.join(
            sub_outdir, args.sample_name + "_" + assembly_identifier)
        if dups:
            out_fastq_r1 = out_fastq_pre + '_unmap_dups_R1.fq'
            out_fastq_r2 = out_fastq_pre + '_unmap_dups_R2.fq'
        else:
            out_fastq_r1 = out_fastq_pre + '_unmap_R1.fq'
            out_fastq_r2 = out_fastq_pre + '_unmap_R2.fq'

        # Use the undeduplicated even for duplicates on re-runs
        out_fastq_r1_gz = out_fastq_pre + '_unmap_R1.fq.gz'
        out_fastq_r2_gz = out_fastq_pre + '_unmap_R2.fq.gz'

        if useFIFO and paired and not args.keep:
            if dups:
                out_fastq_tmp = os.path.join(sub_outdir,
                                             assembly_identifier + "_dups_bt2")
            else:
                out_fastq_tmp = os.path.join(sub_outdir,
                                             assembly_identifier + "_bt2")
            cmd = "mkfifo " + out_fastq_tmp

            if os.path.exists(out_fastq_tmp):
                os.remove(out_fastq_tmp)
            pm.run(cmd, out_fastq_tmp)
            pm.clean_add(out_fastq_tmp)

        else:
            if dups:
                out_fastq_tmp = out_fastq_pre + '_unmap_dups.fq'
            else:
                out_fastq_tmp = out_fastq_pre + '_unmap.fq'

        out_fastq_tmp_gz = out_fastq_pre + '_unmap.fq.gz'

        filter_pair = build_command([tools.perl,
            tool_path("filter_paired_fq.pl"), out_fastq_tmp,
            unmap_fq1, unmap_fq2, out_fastq_r1, out_fastq_r2])
        # TODO: make filter_paired_fq work with SE data
        # cmd = build_command([tools.perl,
           # tool_path("filter_paired_fq.pl"), out_fastq_tmp,
           # unmap_fq1, out_fastq_r1])
        # For now, revert to old method

        if not bt2_opts_txt:
            # Default options
            bt2_opts_txt = "-k 1"  # Return only 1 alignment
            bt2_opts_txt += " -D 20 -R 3 -N 1 -L 20 -i S,1,0.50"

        # samtools sort needs a temporary directory
        tempdir = tempfile.mkdtemp(dir=sub_outdir)
        os.chmod(tempdir, 0o771)
        pm.clean_add(tempdir)

        # Build bowtie2 command
        cmd = "(" + tools.bowtie2 + " -p " + str(pm.cores)
        cmd += " " + bt2_opts_txt
        cmd += " -x " + assembly_bt2
        cmd += " --rg-id " + args.sample_name
        cmd += " -U " + unmap_fq1
        cmd += " --un " + out_fastq_tmp
        if args.keep: #  or not paired
            #cmd += " --un-gz " + out_fastq_bt2 # TODO drop this for paired... because repairing with filter_paired_fq.pl
            # In this samtools sort command we print to stdout and then use > to
            # redirect instead of  `+ " -o " + mapped_bam` because then samtools
            # uses a random temp file, so it won't choke if the job gets
            # interrupted and restarted at this step.
            cmd += " | " + tools.samtools + " view -bS - -@ 1"  # convert to bam
            cmd += " | " + tools.samtools + " sort - -@ 1"  # sort output
            cmd += " -T " + tempdir
            cmd += " -o " + mapped_bam
            cmd += ") 2>&1"
        else:
            cmd += " 2>&1 > /dev/null)"
        #cmd += ")"
        #cmd += ") 2> " + summary_file

        if paired:
            if args.keep or not useFIFO:
                # checkprint() doesn't know how to handle targets
                # must recreate that effect ourselves
                if not _itsa_file(mapped_bam) or args.new_start:
                    aln_stats = pm.checkprint(cmd)
                pm.run(filter_pair, mapped_bam)
            else:
                pm.wait = False
                pm.run(filter_pair, out_fastq_r2_gz)
                pm.wait = True
                if not _itsa_file(out_fastq_r2_gz) or args.new_start:
                    aln_stats = pm.checkprint(cmd)
        else:
            if args.keep:
                if not _itsa_file(mapped_bam) or args.new_start:
                    aln_stats = pm.checkprint(cmd)
            else:
                # TODO: switch to this once filter_paired_fq works with SE
                #pm.run(cmd2, summary_file)
                #pm.run(cmd1, out_fastq_r1)
                if not _itsa_file(out_fastq_tmp_gz) or args.new_start:
                    aln_stats = pm.checkprint(cmd)
                else:
                    aln_stats = None

        pm.clean_add(out_fastq_tmp)

        if not dups:
            if not safely_retrieve(pm, "Aligned_reads_assembly") or args.new_start:
                if aln_stats:
                    pm.info(aln_stats)  # Log alignment statistics
                    try:
                        align_exact = re.search(r".*aligned exactly 1 time", aln_stats).group().split()[0]
                    except AttributeError:
                        align_exact = None
                    except:
                        err_msg = "Unable to determine alignment statistics for {}."
                        pm.fail_pipeline(RuntimeError(err_msg.format(args.genome_assembly)))
                else:
                    align_exact = None
                # cmd = ("grep 'aligned exactly 1 time' " + summary_file +
                #        " | awk '{print $1}'")
                # align_exact = pm.checkprint(cmd)

                if align_exact:
                    if paired:
                        ar = float(align_exact)*2
                    else:
                        ar = float(align_exact)
                else:
                    ar = 0

                # report aligned reads
                #TODO: the pipestat result identifier cannot be variable
                # pm.report_result("Aligned_reads_" + assembly_identifier, ar)
                pm.pipestat.report(values={"Aligned_reads_assembly": ar})
                try:
                    # wrapped in try block in case Trimmed_reads is not reported 
                    # in this pipeline.
                    tr = float(pm.pipestat.retrieve(None,"Trimmed_reads_R1"))
                except:
                    print("Trimmed reads is not reported.")
                else:
                    res_key = "Alignment_rate_" + assembly_identifier
                    #TODO: the pipestat result identifier cannot be variable
                    res_key = "Alignment_rate_assembly"
                    if float(ar) > 0:
                        # pm.pipestat.report(values={res_key, round(float(ar) * 100 / float(tr), 2))
                        pm.pipestat.report(values={res_key: round(float(ar) * 100 / float(tr), 2)})
                    else:
                        pm.report_result(res_key, 0)
                        pm.pipestat.report(values={res_key: 0})
        
        if paired:
            unmap_fq1 = out_fastq_r1
            unmap_fq2 = out_fastq_r2
        else:
            # Use alternate once filter_paired_fq is working with SE
            #unmap_fq1 = out_fastq_r1
            unmap_fq1 = out_fastq_tmp
            unmap_fq2 = ""

        return unmap_fq1, unmap_fq2
    else:
        msg = "No {} index found in {}; skipping.".format(
            assembly_identifier, os.path.dirname(assembly_bt2))
        print(msg)
        return unmap_fq1, unmap_fq2


def tool_path(tool_name):
    """
    Return the path to a tool used by this pipeline.

    :param str tool_name: name of the tool (e.g., a script file_name)
    :return str: real, absolute path to tool (expansion and symlink resolution)
    """

    return os.path.join(os.path.dirname(os.path.dirname(__file__)),
                        TOOLS_FOLDER, tool_name)


def _guess_encoding(fq):
    """
    Adapted from Brent Pedersen's "_guess_encoding.py"
    https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py
    
    Copyright (c) 2018 Brent Pedersen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    Guess the encoding of a stream of qual lines.
    """
    RANGES = {
        'Sanger': (33, 73),
        'Illumina-1.8': (33, 74),
        'Solexa': (59, 104),
        'Illumina-1.3': (64, 104),
        'Illumina-1.5': (67, 105)
    }
    
    def get_qual_range(qual_str):
        vals = [ord(c) for c in qual_str]
        return min(vals), max(vals)
    
    def get_encodings_in_range(rmin, rmax, ranges=RANGES):
        valid_encodings = []
        for encoding, (emin, emax) in ranges.items():
            if rmin >= emin and rmax <= emax:
                valid_encodings.append(encoding)
        return valid_encodings

    gmin = 99
    gmax = 0
    valid = []

    err_exit = False

    with open(fq) as fastq_file:
        for line_num, line in enumerate(fastq_file):
            # Python starts at 0; need to start at 1 for this step
            if (line_num + 1) % 4 == 0:
                lmin, lmax = get_qual_range(line.rstrip())

                if lmin < gmin or lmax > gmax:
                    gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                    valid = get_encodings_in_range(gmin, gmax)

                    if len(valid) == 0:
                        print("no encodings for range: "
                              "{}".format((gmin, gmax)))
                        err_exit = True
                        continue

                    if len(valid) == 1:
                        err_exit = False
                        break

    if err_exit:
        return("Unknown")
    else:
        return(str(valid[-1]))


def _check_commands(commands, ignore=''):
    """
    Check if command(s) can be called

    :param attributedict commands: dictionary of commands to check
    :param list ignore: list of commands that are optional and can be ignored
    """

    # Use `command` to see if command is callable, store exit code
    is_callable = True
    uncallable = []
    for name, command in commands.items():
        if command not in ignore:
            # if an environment variable is not expanded it means it points to
            # an uncallable command
            if '$' in command:
                # try to expand
                command = os.path.expandvars(os.path.expanduser(command))
                if not os.path.exists(command):
                    uncallable.append(command)

            # if a command is a java file, modify the command
            if '.jar' in command:
                command = "java -jar " + command

            code = os.system("command -v {0} >/dev/null 2>&1 || {{ exit 1; }}".format(command))
            # If exit code is not 0, track which command failed
            #print("{} code {}".format(command, code))  # DEBUG
            if code != 0:
                uncallable.append(command)
                is_callable = False
    if is_callable:
        return True
    else:
        print("The following required tool(s) are not callable: {0}".format(' '.join(uncallable)))
        return False


def calc_frip(bamfile, ftfile, frip_func, pipeline_manager,
              aligned_reads_key="Aligned_reads"):
    """
    Calculate the fraction of reads in feature file.

    Use the given function and data from an aligned reads file and a features 
    file, along with a PipelineManager, to calculate.

    :param str bamfile: path to aligned reads file
    :param str ftfile: path to features of interest to calculate overlap
    :param callable frip_func: how to calculate the fraction of reads in feat;
        this must accept the path to the aligned reads file and the path to
        the called peaks file as arguments.
    :param pypiper.PipelineManager pipeline_manager: the PipelineManager in use
        for the pipeline calling this function
    :param str aligned_reads_key: name of the key from a stats (key-value) file
        to use to fetch the count of aligned reads
    :return float: fraction of reads in peaks
    """
    frip_cmd = frip_func(bamfile, ftfile)
    num_in_reads = pipeline_manager.checkprint(frip_cmd)
    num_aligned_reads = pipeline_manager.get_stat(aligned_reads_key)
    print(num_aligned_reads, num_in_reads)
    return float(num_in_reads) / float(num_aligned_reads)


def _itsa_file(anyfile):
    """
    Helper function to confirm a file exists and is not empty.

    :param str anyfile: path to a file
    """
    return(os.path.isfile(anyfile) and os.stat(anyfile).st_size > 0)


def _itsa_empty_file(anyfile):
    """
    Helper function to confirm a file exists but is empty.

    :param str anyfile: path to a file
    """
    return(os.path.isfile(anyfile) and os.stat(anyfile).st_size == 0)


def is_gzipped(file_name):
    """
    Determine whether indicated file appears to be gzipped.
    :param str file_name: Name/path of file to check as gzipped.
    :return bool: Whether indicated file appears to be in gzipped format.
    """
    _, ext = os.path.splitext(file_name)
    return file_name.endswith(".gz")


def _add_resources(args, res, asset_dict=None):
    """
    Add additional resources needed for pipeline.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param pm.config.resources res: pipeline manager resources list
    :param asset_dict list: list of dictionary of assets to add
    """

    rgc = RGC(select_genome_config(res.get("genome_config")))

    key_errors = []
    exist_errors = []
    required_list = []

    # Check that bowtie2 indicies exist for specified prealignments
    for reference in args.prealignments:
        for asset in [BT2_IDX_KEY]:
            try:
                res[asset] = rgc.seek(reference, asset)
            except KeyError:
                err_msg = "{} for {} is missing from REFGENIE config file."
                pm.fail_pipeline(KeyError(err_msg.format(asset, reference)))
            except:
                err_msg = "{} for {} does not exist."
                pm.fail_pipeline(IOError(err_msg.format(asset, reference)))

    # Check specified assets
    if not asset_dict:
        return res, rgc
    else:
        for item in asset_dict:
            pm.debug("item: {}".format(item))  # DEBUG
            asset = item["asset_name"]
            seek_key = item["seek_key"] or item["asset_name"]
            tag = item["tag_name"] or "default"
            arg = item["arg"]
            user_arg = item["user_arg"]
            req = item["required"]

            if arg and hasattr(args, arg) and getattr(args, arg):
                res[seek_key] = os.path.abspath(getattr(args, arg))
            else:
                try:
                    pm.debug("{} - {}.{}:{}".format(args.genome_assembly,
                                                    asset,
                                                    seek_key,
                                                    tag))  # DEBUG
                    res[seek_key] = rgc.seek(args.genome_assembly,
                                                  asset_name=str(asset),
                                                  tag_name=str(tag),
                                                  seek_key=str(seek_key))
                except KeyError:
                    key_errors.append(item)
                    if req:
                        required_list.append(item)
                except:
                    exist_errors.append(item)
                    if req:
                        required_list.append(item)

        if len(key_errors) > 0 or len(exist_errors) > 0:
            pm.info("Some assets are not found. You can update your REFGENIE "
                    "config file or point directly to the file using the noted "
                    "command-line arguments:")

        if len(key_errors) > 0:
            if required_list:
                err_msg = "Required assets missing from REFGENIE config file: {}"
                pm.fail_pipeline(IOError(err_msg.format(", ".join(["{asset_name}.{seek_key}:{tag_name}".format(**x) for x in required_list]))))
            else:
                warning_msg = "Optional assets missing from REFGENIE config file: {}"
                pm.info(warning_msg.format(", ".join(["{asset_name}.{seek_key}:{tag_name}".format(**x) for x in key_errors])))

        if len(exist_errors) > 0:
            if required_list:
                err_msg = "Required assets not existing: {}"
                pm.fail_pipeline(IOError(err_msg.format(", ".join(["{asset_name}.{seek_key}:{tag_name} (--{user_arg})".format(**x) for x in required_list]))))
            else:
                warning_msg = "Optional assets not existing: {}"
                pm.info(warning_msg.format(", ".join(["{asset_name}.{seek_key}:{tag_name} (--{user_arg})".format(**x) for x in exist_errors])))

        return res, rgc


def report_message(pm, report_file, message, annotation=None):
    """
    Writes a string to provided file in a safe way.
    
    :param PipelineManager pm: a pypiper PipelineManager object
    :param str report_file: name of the output file
    :param str message: string to write to the output file
    :param str annotation: By default, the message will be annotated with the
        pipeline name, so you can tell which pipeline records which stats.
        If you want, you can change this; use annotation='shared' if you
        need the stat to be used by another pipeline (using get_stat()).
    """
    # Default annotation is current pipeline name.
    annotation = str(annotation or pm.name)

    message = str(message).strip()
    
    message = "{message}\t{annotation}".format(
        message=message, annotation=annotation)

    # Just to be extra careful, let's lock the file while we we write
    # in case multiple pipelines write to the same file.
    pm._safe_write_to_file(report_file, message)


###############################################################################
def main():
    """
    Main pipeline process.
    """

    args = parse_arguments()

    args.paired_end = args.single_or_paired.lower() == "paired"

    # Initialize, creating global PipelineManager and NGSTk instance for
    # access in ancillary functions outside of main().
    outfolder = os.path.abspath(
        os.path.join(args.output_parent, args.sample_name))
    global pm
    pm = pypiper.PipelineManager(
        name="PEPPRO", outfolder=outfolder, args=args, version=__version__, 
        pipestat_schema=args.pipestat_schema,
        pipestat_results_file=args.pipestat_results_file,
        pipestat_record_id=args.pipestat_record_id,
        pipestat_namespace=args.pipestat_namespace,
        pipestat_config=args.pipestat_config,)
    global ngstk
    ngstk = pypiper.NGSTk(pm=pm)

    pm.pipestat.report(
        values={
            "log": {"path": pm.pipeline_log_file, "title": "Pipeline log file"},
            "profile": {"path": pm.pipeline_profile_file, "title": "Pipeline profile file"},
            "commands": {"path": pm.pipeline_commands_file, "title": "Pipeline commands file"},
            "version": pm.pl_version
        }
    ) 

    # Convenience alias
    tools = pm.config.tools
    param = pm.config.parameters
    res   = pm.config.resources
    #sstructure = pm.sample_structure  # maybe possible in the future?

    # Check that the required tools are callable by the pipeline
    tool_list = [v for k,v in tools.items()]    # extract tool list
    tool_list = [t.replace('fastx', 'fastx_trimmer') for t in tool_list]
    tool_list = [t.replace('seqoutbias', 'seqOutBias') for t in tool_list]
    opt_tools = ["fqdedup", "fastx_trimmer", "seqOutBias", "fastqc"]
    if args.trimmer == "fastx":  # update tool call
        if 'fastx' in opt_tools: opt_tools.remove('fastx_trimmer')

    if args.dedup == "fqdedup":  # update tool call
        if 'fqdedup' in opt_tools: opt_tools.remove('fqdedup')

    if args.sob:
        if 'seqOutBias' in opt_tools: opt_tools.remove('seqOutBias')

    tool_list = dict((t,t) for t in tool_list)  # convert back to dict

    if not _check_commands(tool_list, opt_tools):
        err_msg = "Missing required tools. See message above."
        pm.fail_pipeline(RuntimeError(err_msg))

    if args.input2 and not args.paired_end:
        err_msg = "Incompatible settings: You specified single-end, but provided --input2."
        pm.fail_pipeline(RuntimeError(err_msg))

    # Set up reference resource according to genome prefix.
    check_list = [
        {"asset_name":"fasta", "seek_key":"chrom_sizes",
         "tag_name":"default", "arg":None, "user_arg":None,
         "required":True},
        {"asset_name":"fasta", "seek_key":None,
         "tag_name":"default", "arg":None, "user_arg":None,
         "required":True},
        {"asset_name":BT2_IDX_KEY, "seek_key":None,
         "tag_name":"default", "arg":None, "user_arg":None,
         "required":True}
    ]
    # If user specifies TSS file, use that instead of the refgenie asset
    if not (args.TSS_name):
        check_list.append(
            {"asset_name":"refgene_anno", "seek_key":"refgene_tss",
             "tag_name":"default", "arg":"TSS_name", "user_arg":"TSS-name",
             "required":False}
        )
    # If user specifies a custom pause index TSS file, use that instead
    if not (args.ensembl_tss):
        check_list.append(
            {"asset_name":"ensembl_gtf", "seek_key":"ensembl_tss",
             "tag_name":"default", "arg":"ensembl_tss", "user_arg":"pi-tss",
             "required":False}
        )
    # If user specifies a custom pause index gene body file, use that instead
    if not (args.ensembl_gene_body):
        check_list.append(
            {"asset_name":"ensembl_gtf", "seek_key":"ensembl_gene_body",
             "tag_name":"default", "arg":"ensembl_gene_body",
             "user_arg":"pi-body", "required":False}
        )
    # If user specifies a custom premature RNA file, use that instead
    if not (args.pre_name):
        check_list.append(
            {"asset_name":"refgene_anno", "seek_key":"refgene_pre_mRNA",
             "tag_name":"default", "arg":"pre_name", "user_arg":"pre-name",
             "required":False}
        )
    # If user specifies feature annotation file,
    # use that instead of the refgenie managed asset
    if not (args.anno_name):
        check_list.append(
            {"asset_name":"feat_annotation", "seek_key":"feat_annotation",
            "tag_name":"default", "arg":"anno_name", "user_arg":"anno-name",
            "required":False}
        )
    # If user specifies a custom exon file, use that instead
    if not (args.exon_name):
        check_list.append(
            {"asset_name":"refgene_anno", "seek_key":"refgene_exon",
             "tag_name":"default", "arg":"exon_name", "user_arg":"exon-name",
             "required":False}
        )
    # If user specifies a custom intron file, use that instead
    if not (args.intron_name):
        check_list.append(
            {"asset_name":"refgene_anno", "seek_key":"refgene_intron",
             "tag_name":"default", "arg":"intron_name",
             "user_arg":"intron-name", "required":False}
        )
    res, rgc = _add_resources(args, res, check_list)

    # If the user specifies optional files, add those to our resources
    if ((args.TSS_name) and os.path.isfile(args.TSS_name) and
            os.stat(args.TSS_name).st_size > 0):
        res.refgene_tss = args.TSS_name
    if ((args.ensembl_tss) and os.path.isfile(args.ensembl_tss) and
            os.stat(args.ensembl_tss).st_size > 0):
        res.ensembl_tss = args.ensembl_tss
    if ((args.ensembl_gene_body) and os.path.isfile(args.ensembl_gene_body) and
            os.stat(args.ensembl_gene_body).st_size > 0):
        res.ensembl_gene_body = args.ensembl_gene_body
    if ((args.pre_name) and os.path.isfile(args.pre_name) and
            os.stat(args.pre_name).st_size > 0):
        res.refgene_pre_mRNA = args.pre_name
    if ((args.anno_name) and os.path.isfile(args.anno_name) and
            os.stat(args.anno_name).st_size > 0):
        res.feat_annotation = args.anno_name
    if ((args.exon_name) and os.path.isfile(args.exon_name) and
            os.stat(args.exon_name).st_size > 0):
        res.refgene_exon = args.exon_name
    if ((args.intron_name) and os.path.isfile(args.intron_name) and
            os.stat(args.intron_name).st_size > 0):
        res.refgene_intron = args.intron_name

    # Adapter file can be set in the config; if left null, we use a default.
    # Expects headers to include >5prime and >3prime
    res.adapters = res.adapters or tool_path("adapter.fa")
    param.outfolder = outfolder
    
    # Report utilized assets
    assets_file = os.path.join(param.outfolder, "assets.tsv")
    for asset in res:
        message = "{}\t{}".format(asset, os.path.expandvars(res[asset]))
        report_message(pm, assets_file, message)
        
    # Report primary genome
    message = "genome\t{}".format(args.genome_assembly)
    report_message(pm, assets_file, message)

    ###########################################################################
    #          Check that the input file(s) exist before continuing           #
    ###########################################################################
    if _itsa_file(args.input[0]):
        print("Local input file: " + args.input[0])
    elif _itsa_empty_file(args.input[0]):
        # The read1 file exists but is empty
        err_msg = "File exists but is empty: {}"
        pm.fail_pipeline(IOError(err_msg.format(args.input[0])))
    else:
        # The read1 file does not exist
        err_msg = "Could not find: {}"
        pm.fail_pipeline(IOError(err_msg.format(args.input[0])))

    if args.input2:
        if _itsa_file(args.input2[0]):
            print("Local input file: " + args.input2[0])
        elif _itsa_empty_file(args.input2[0]):
            # The read1 file exists but is empty
            err_msg = "File exists but is empty: {}"
            pm.fail_pipeline(IOError(err_msg.format(args.input2[0])))
        else:
            # The read1 file does not exist
            err_msg = "Could not find: {}"
            pm.fail_pipeline(IOError(err_msg.format(args.input2[0])))

    container = None # legacy

    ###########################################################################
    #                      Grab and prepare input files                       #
    ###########################################################################
    # pm.report_result(
    #     "File_mb",
    #     round(ngstk.get_file_size([x for x in [args.input, args.input2] if x is not None]), 2)
    # )
    # pm.report_result("Read_type", args.single_or_paired)
    # pm.report_result("Genome", args.genome_assembly)


    pm.pipestat.report(values={
        "File_mb": round(ngstk.get_file_size([x for x in [args.input, args.input2] if x is not None]), 2),
        "Read_type": args.single_or_paired,
        "Genome": args.genome_assembly,
    })

    # PRO-seq pipeline
    if args.protocol.lower() in RUNON_SOURCE_GRO:
        pm.info("Detected GRO input")
    elif args.protocol.lower() in RUNON_SOURCE_PRO:
        pm.info("Detected PRO input")
    else:
        pm.fail_pipeline(RuntimeError("Input protocol must be GRO or PRO."))

    # Each (major) step should have its own subfolder
    raw_folder = os.path.join(param.outfolder, "raw")
    fastq_folder = os.path.join(param.outfolder, "fastq")
    fastqc_folder=os.path.join(param.outfolder, "fastqc")
    ngstk.make_dir(fastqc_folder)

    pm.timestamp("### Merge/link and fastq conversion: ")
    # This command will merge multiple inputs so you can use multiple
    # sequencing lanes in a single pipeline run.
    local_input_files = ngstk.merge_or_link(
        [args.input, args.input2], raw_folder, args.sample_name)
    cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(
        local_input_files, args.sample_name, args.paired_end, fastq_folder)

    #if not pm.pipestat.retrieve(None,"Raw_reads") or args.new_start:
    # TODO: improve the skipping of these steps on recovery runs
    #       issue here is that process_fastq is still trying to run
    #       if we skip this step
    pm.run(cmd, unaligned_fastq,
           follow=ngstk.check_fastq(
               local_input_files, unaligned_fastq, args.paired_end))
    pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)

    pm.info(local_input_files)

    untrimmed_fastq1 = out_fastq_pre + "_R1.fastq"
    untrimmed_fastq2 = out_fastq_pre + "_R2.fastq" if args.paired_end else None

    ############################################################################
    #                          Process read files                              #
    ############################################################################
    pm.timestamp("### FASTQ processing: ")

    unmap_fq1 = out_fastq_pre + '_unmap_R1.fq'
    unmap_fq2 = out_fastq_pre + '_unmap_R2.fq'
    unmap_fq1_dups = out_fastq_pre + '_unmap_dups_R1.fq'
    unmap_fq2_dups = out_fastq_pre + '_unmap_dups_R2.fq'

    cutadapt_folder = os.path.join(outfolder, "cutadapt")
    cutadapt_report = os.path.join(cutadapt_folder,
                                   args.sample_name + "_R1_cutadapt.txt")

    processed_target_R1 = os.path.join(fastq_folder, "processed_R1.flag")
    processed_target_R2 = os.path.join(fastq_folder, "processed_R2.flag")
    rmUMI_target = os.path.join(fastq_folder, "readname_repaired.flag")
    rmUMI_dups_target = os.path.join(fastq_folder, "readname_dups_repaired.flag")
    repair_target = os.path.join(fastq_folder, "repaired.flag")
    dups_repair_target = os.path.join(fastq_folder, "dups_repaired.flag")

    # If single-end, must use cutadapt for plotting purposes
    if not args.paired_end:
        if args.adapter != "cutadapt":
            pm.warning("You set adapter arg to '{}' but you must select "
                       "'cutadapt' to plot the adapter insertion distribution "
                       "for single end data.".format(args.adapter))

    if args.paired_end:
        if not args.complexity and int(args.umi_len) > 0:
            if not os.path.exists(processed_target_R1) or args.new_start:
                unmap_fq1, unmap_fq1_dups = _process_fastq(
                    args, tools, res, False,
                    untrimmed_fastq1, outfolder=param.outfolder)
            cmd = ("touch " + processed_target_R1)
            pm.run(cmd, processed_target_R1)
        else:
            if not os.path.exists(processed_target_R1) or args.new_start:
                unmap_fq1 = _process_fastq(
                    args, tools, res, False,
                    untrimmed_fastq1, outfolder=param.outfolder)
            cmd = ("touch " + processed_target_R1)
            pm.run(cmd, processed_target_R1)

        if not os.path.exists(processed_target_R2) or args.new_start:
            unmap_fq2, unmap_fq2_dups = _process_fastq(
                args, tools, res, True,
                untrimmed_fastq2, outfolder=param.outfolder)
        cmd = ("touch " + processed_target_R2)
        pm.run(cmd, processed_target_R2)

        pm.debug("\n\nunmap_fq1: {}\nunmap_fq2: {}\n\n".format(unmap_fq1, unmap_fq2))

        # Gut check
        # Processing fastq should have trimmed the reads.
        tmp = pm.pipestat.retrieve(None,"Trimmed_reads_R1")
        if tmp:
            tr = float(tmp)
        else:
            tr = 0
        if (tr < 1):
            pm.fail_pipeline(RuntimeError("No reads left after trimming. Check trimmer settings"))

        # Re-pair fastq files
        r1_repair = os.path.join(
            fastq_folder, args.sample_name + "_R1_processed.fastq.paired.fq")
        r2_repair = os.path.join(
            fastq_folder, args.sample_name + "_R2_trimmed.fastq.paired.fq")

        r1_repair_single = os.path.join(
            fastq_folder, args.sample_name + "_R1_processed.fastq.single.fq")
        r2_repair_single = os.path.join(
            fastq_folder, args.sample_name + "_R2_trimmed.fastq.single.fq")

        tmp = float(pm.pipestat.retrieve(None,"Raw_reads"))
        if tmp:
            rr = float(tmp)
        else:
            rr = 0
        if (rr < 1):
            pm.fail_pipeline(RuntimeError("Raw_reads were not reported. Check output ({})".format(param.outfolder)))

        if args.adapter == "fastp" and int(args.umi_len) > 0:
            noUMI_fq1 = os.path.join(fastq_folder,
                args.sample_name + "_R1_processed_noUMI.fastq")
            noUMI_fq2 = os.path.join(fastq_folder,
                args.sample_name + "_R2_trimmed_noUMI.fastq")
            cmd1 = ("sed -e 's|\\:[^:]*\\([[:space:]].*\\)|\\1 |g'" +
                    " " + unmap_fq1 + " > " + noUMI_fq1)
            cmd2 = ("sed -e 's|\\:[^:]*\\([[:space:]].*\\)|\\1 |g'" +
                    " " + unmap_fq2 + " > " + noUMI_fq2)
            pm.run([cmd1, cmd2], rmUMI_target, shell=True)

            cmd1 = ("mv " + noUMI_fq1 + " " + unmap_fq1)
            cmd2 = ("mv " + noUMI_fq2 + " " + unmap_fq2)
            cmd3 = ("touch " + rmUMI_target)
            pm.run([cmd1, cmd2, cmd3], rmUMI_target)

        cmd1 = (tools.fastqpair + " -t " + str(int(0.9*rr)) + " " + 
                unmap_fq1 + " " + unmap_fq2)
        cmd2 = ("mv " + r1_repair + " " + unmap_fq1)
        cmd3 = ("mv " + r2_repair + " " + unmap_fq2)
        cmd4 = ("touch " + repair_target)
        pm.run([cmd1, cmd2, cmd3, cmd4], repair_target)
        pm.clean_add(r1_repair_single)
        pm.clean_add(r2_repair_single)

        # Re-pair the duplicates (but only if we could identify duplicates)
        if int(args.umi_len) > 0:
            r1_dups_repair = os.path.join(
                fastq_folder, args.sample_name + "_R1_trimmed.fastq.paired.fq")
            r2_dups_repair = os.path.join(
                fastq_folder, args.sample_name + "_R2_trimmed_dups.fastq.paired.fq")

            r1_dups_repair_single = os.path.join(
                fastq_folder, args.sample_name + "_R1_trimmed.fastq.single.fq")
            r2_dups_repair_single = os.path.join(
                fastq_folder, args.sample_name + "_R2_trimmed_dups.fastq.single.fq")

            if args.adapter == "fastp" and int(args.umi_len) > 0:
                noUMI_fq1_dups = os.path.join(fastq_folder,
                    args.sample_name + "_R1_trimmed_dups_noUMI.fastq")
                noUMI_fq2_dups = os.path.join(fastq_folder,
                    args.sample_name + "_R2_trimmed_dups_noUMI.fastq")
                cmd1 = ("sed -e 's|\\:[^:]*\\([[:space:]].*\\)|\\1 |g'" +
                        " " + unmap_fq1_dups + " > " + noUMI_fq1_dups)
                cmd2 = ("sed -e 's|\\:[^:]*\\([[:space:]].*\\)|\\1 |g'" +
                        " " + unmap_fq2_dups + " > " + noUMI_fq2_dups)
                pm.run([cmd1, cmd2], [noUMI_fq1_dups, noUMI_fq2_dups], shell=True)
                cmd1 = ("mv " + noUMI_fq1_dups + " " + unmap_fq1_dups)
                cmd2 = ("mv " + noUMI_fq2_dups + " " + unmap_fq2_dups)
                cmd3 = ("touch " + rmUMI_dups_target)
                pm.run([cmd1, cmd2, cmd3], rmUMI_dups_target)

            cmd1 = (tools.fastqpair + " -t " + str(int(0.9*rr)) + " " +
                    unmap_fq1_dups + " " + unmap_fq2_dups)
            cmd2 = ("mv " + r1_dups_repair + " " + unmap_fq1_dups)
            cmd3 = ("mv " + r2_dups_repair + " " + unmap_fq2_dups)
            cmd4 = ("touch " + dups_repair_target)
            pm.run([cmd1, cmd2, cmd3, cmd4], dups_repair_target)
            pm.clean_add(r1_dups_repair_single)
            pm.clean_add(r2_dups_repair_single)
    else:
        if not args.complexity and int(args.umi_len) > 0:
            if not os.path.exists(processed_target_R1) or args.new_start:
                unmap_fq1, unmap_fq1_dups = _process_fastq(
                    args, tools, res, False,
                    untrimmed_fastq1, outfolder=param.outfolder)
                unmap_fq2 = ""
                unmap_fq2_dups = ""
            cmd = ("touch " + processed_target_R1)
            pm.run(cmd, processed_target_R1)
        else:
            if not os.path.exists(processed_target_R1) or args.new_start:
                unmap_fq1 = _process_fastq(
                    args, tools, res, False,
                    untrimmed_fastq1, outfolder=param.outfolder)
                unmap_fq2 = ""
            cmd = ("touch " + processed_target_R1)
            pm.run(cmd, processed_target_R1)

    # NOTE: maintain this functionality for single-end data
    #       for paired-end it has already been generated at this point
    if not args.paired_end:
        pm.timestamp("### Plot adapter insertion distribution")

        if not args.adapter == "cutadapt":
            pm.info("Skipping adapter insertion distribution plotting...")
            pm.info("For SE data, this requires using 'cutadapt' for adapter removal.")
        elif not os.path.exists(cutadapt_report):
            pm.info("Skipping adapter insertion distribution plotting...")
            pm.info("Could not find {}.`".format(cutadapt_report))
        else:
            degradation_pdf = os.path.join(cutadapt_folder,
                args.sample_name + "_R1_adapter_insertion_distribution.pdf")
            degradation_png = os.path.join(cutadapt_folder,
                args.sample_name + "_R1_adapter_insertion_distribution.png")
            cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") + 
                   " cutadapt -i " + cutadapt_report + " -o " + cutadapt_folder)
            if int(args.umi_len) > 0:
                cmd += (" -u " + str(args.umi_len))
                umi_len = args.umi_len
            else:
                umi_len = 0

            pm.run(cmd, degradation_pdf, nofail=True)
            # pm.report_object("Adapter insertion distribution", degradation_pdf,
            #                  anchor_image=degradation_png)
            pm.pipestat.report(values={"Adapter_insertion_distribution": {"path": degradation_pdf, "thumbnail_path": degradation_png, "title": "Adapter insertion distribution"}})

            if not safely_retrieve(pm, "Peak_adapter_insertion_size") or args.new_start:
                # Determine the peak insertion size
                cmd = ("awk '/count/,0' " + cutadapt_report +
                       " | awk 'NR>2 {print prev} {prev=$0}'" +
                       " | awk '{if ($3/$2 < 0.01) print $1, $2}'" +
                       " | awk 'BEGIN{max=   0; max_len=0; len=0}" +
                       "{if ($2>0+max) {max=$2; len=$1}; max_len=$1} " +
                       "END{print max_len-len}'")
                adapter_peak = pm.checkprint(cmd)
                if adapter_peak:
                    ap = int(adapter_peak)
                    # pm.report_result("Peak_adapter_insertion_size", ap)
                    pm.pipestat.report(values={"Peak_adapter_insertion_size": ap})

            # Calculate the degradation ratio
            if not safely_retrieve(pm, "Degradation_ratio") or args.new_start:
                pm.timestamp("###  Calculating degradation ratio")

                cmd = ("awk 'NR>2 {print prev} {prev=$0}' " + cutadapt_report +
                       " | awk '{ if ($1 == 10) {status = 1}} END " + 
                       "{if (status) {print status} else {print 0}}'")
                degraded_lower = pm.checkprint(cmd)
                cmd = ("awk 'NR>2 {print prev} {prev=$0}' " + cutadapt_report +
                       " | awk '{ if ($1 == 20) {status = 1}} END " + 
                       "{if (status) {print status} else {print 0}}'")
                degraded_upper = pm.checkprint(cmd)
                cmd = ("awk 'NR>2 {print prev} {prev=$0}' " + cutadapt_report +
                       " | awk '{ if ($1 == 30) {status = 1}} END " + 
                       "{if (status) {print status} else {print 0}}'")
                intact_lower = pm.checkprint(cmd)
                cmd = ("awk 'NR>2 {print prev} {prev=$0}' " + cutadapt_report +
                       " | awk '{ if ($1 == 40) {status = 1}} END " + 
                       "{if (status) {print status} else {print 0}}'")
                intact_upper = pm.checkprint(cmd)

                if degraded_lower:
                    dl = int(degraded_lower)
                if dl == 1:
                    dl = 10
                else:
                    cmd = ("awk 'NR>2 {print prev} {prev=$0}' " +
                           cutadapt_report + " | awk 'NR==1 {print $1}'")
                    degraded_lower = pm.checkprint(cmd)
                    dl = int(degraded_lower) if degraded_lower else 1

                if degraded_upper:
                    du = int(degraded_upper)
                if du == 1:
                    du = 20
                else:
                    du = int(degraded_lower) + 9

                if intact_upper:
                    iu = int(intact_upper)
                if iu == 1:
                    iu = 40
                else:
                    cmd = ("awk 'NR>2 {print prev} {prev=$0}' " +
                           cutadapt_report + " | awk 'END {print $1}'")
                    intact_upper = pm.checkprint(cmd)
                    dl = int(intact_upper) if intact_upper else 40

                if intact_lower:
                    il = int(intact_lower)
                if il == 1:
                    il = 30
                else:
                    il = int(intact_upper) - 10

                cmd = ("awk '/count/,0' " + cutadapt_report +
                       " | awk 'NR>2 {print prev} {prev=$0}'" +
                       " | awk '{if ($3/$2 < 0.01) print $1, $2}'" +
                       " | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}" +
                       "{if ($1 > max_len) {max_len=$1}} " +
                       "END{ for (i in a) print 1+max_len-a[i], b[i]}'" +
                       " | sort -nk1 | awk '($1 <= " + str(du) +
                       " && $1 >= " + str(dl) + "){degradedSum += $2}; " +
                       "($1 >= " + str(il) + " && $1 <= " + str(iu) +
                       "){intactSum += $2} END {if (intactSum < 1) " +
                       "{intactSum = 1} print degradedSum/intactSum}'")
                degradation_ratio = pm.checkprint(cmd)
                if degradation_ratio:
                    dr = float(degradation_ratio)
                    pm.pipestat.report(values={"Degradation_ratio": round(dr, 4)})
                    # pm.report_result("Degradation_ratio", round(dr, 4))

    pm.clean_add(fastq_folder, conditional=True)

    ############################################################################
    #                  Map to any requested prealignments                      #
    ############################################################################
    # We recommend mapping to human_rDNA first for PRO-seq data
    pm.timestamp("### Prealignments")

    to_compress = []
    #if not pm.pipestat.retrieve(None,"Aligned_reads") or args.new_start:
    if len(args.prealignments) == 0:
        print("You may use `--prealignments` to align to references before "
              "the genome alignment step. See docs.")
    else:
        print("Prealignment assemblies: " + str(args.prealignments))
        # Loop through any prealignment references and map to them sequentially
        for reference in args.prealignments:
            bt2_index = rgc.seek(reference, BT2_IDX_KEY)
            if not bt2_index.endswith(reference):
                bt2_index = os.path.join(
                    rgc.seek(reference, BT2_IDX_KEY), reference)
            if not args.complexity and int(args.umi_len) > 0:
                if args.no_fifo:
                    unmap_fq1, unmap_fq2 = _align_with_bt2(
                        args, tools, args.paired_end, False, unmap_fq1,
                        unmap_fq2, reference,
                        assembly_bt2=bt2_index,
                        outfolder=param.outfolder,
                        aligndir="prealignments",
                        bt2_opts_txt=param.bowtie2_pre.params)

                    unmap_fq1_dups, unmap_fq2_dups = _align_with_bt2(
                        args, tools, args.paired_end, False, unmap_fq1_dups,
                        unmap_fq2_dups, reference,
                        assembly_bt2=bt2_index,
                        outfolder=param.outfolder,
                        aligndir="prealignments",
                        dups=True,
                        bt2_opts_txt=param.bowtie2_pre.params)
                else:
                    unmap_fq1, unmap_fq2 = _align_with_bt2(
                        args, tools, args.paired_end, True, unmap_fq1,
                        unmap_fq2, reference,
                        assembly_bt2=bt2_index,
                        outfolder=param.outfolder,
                        aligndir="prealignments",
                        bt2_opts_txt=param.bowtie2_pre.params)

                    unmap_fq1_dups, unmap_fq2_dups = _align_with_bt2(
                        args, tools, args.paired_end, True, unmap_fq1_dups,
                        unmap_fq2_dups, reference,
                        assembly_bt2=bt2_index,
                        outfolder=param.outfolder,
                        aligndir="prealignments",
                        dups=True,
                        bt2_opts_txt=param.bowtie2_pre.params)

                if args.paired_end:
                    to_compress.append(unmap_fq1_dups)
                    to_compress.append(unmap_fq2_dups)
                    to_compress.append(unmap_fq1)
                    to_compress.append(unmap_fq2)
                else:
                    to_compress.append(unmap_fq1_dups)
                    to_compress.append(unmap_fq1)
            else:
                if args.no_fifo:
                    unmap_fq1, unmap_fq2 = _align_with_bt2(
                        args, tools, args.paired_end, False,
                        unmap_fq1, unmap_fq2, reference,
                        assembly_bt2=bt2_index,
                        outfolder=param.outfolder,
                        aligndir="prealignments",
                        bt2_opts_txt=param.bowtie2_pre.params)
                else:
                    unmap_fq1, unmap_fq2 = _align_with_bt2(
                        args, tools, args.paired_end, True,
                        unmap_fq1, unmap_fq2, reference,
                        assembly_bt2=bt2_index,
                        outfolder=param.outfolder,
                        aligndir="prealignments",
                        bt2_opts_txt=param.bowtie2_pre.params)
                if args.paired_end:
                    to_compress.append(unmap_fq1)
                    to_compress.append(unmap_fq2)
                else:
                    to_compress.append(unmap_fq1)

    ############################################################################
    #                           Map to primary genome                          #
    ############################################################################
    pm.timestamp("### Map to genome")

    # Set up named files and options
    map_genome_folder = os.path.join(
        param.outfolder, "aligned_" + args.genome_assembly)
    ngstk.make_dir(map_genome_folder)

    # Deduplicated alignment files
    mapping_genome_bam = os.path.join(
        map_genome_folder, args.sample_name + "_sort.bam")
    mapping_genome_bam_temp = os.path.join(
        map_genome_folder, args.sample_name + "_temp.bam")
    failQC_genome_bam = os.path.join(
        map_genome_folder, args.sample_name + "_fail_qc.bam")
    unmap_genome_bam = os.path.join(
        map_genome_folder, args.sample_name + "_unmap.bam")

    # Alignment files with duplicates (for library complexity)
    mapping_genome_bam_dups = os.path.join(
        map_genome_folder, args.sample_name + "_sort_dups.bam")
    mapping_genome_bam_temp_dups = os.path.join(
        map_genome_folder, args.sample_name + "_temp_dups.bam")
    failQC_genome_bam_dups = os.path.join(
        map_genome_folder, args.sample_name + "_fail_qc_dups.bam")
    unmap_genome_bam_dups = os.path.join(
        map_genome_folder, args.sample_name + "_unmap_dups.bam")

    # For quality control output
    QC_folder = os.path.join(param.outfolder, "QC_" + args.genome_assembly)
    ngstk.make_dir(QC_folder)

    # For library complexity (for --recover)
    preseq_pdf = os.path.join(
        QC_folder, args.sample_name + "_preseq_plot.pdf")

    temp_mapping_index = os.path.join(mapping_genome_bam_temp + ".bai")
    temp_mapping_index_dups = os.path.join(mapping_genome_bam_temp_dups + ".bai")

    mito_name = ["chrM", "chrMT", "M", "MT", "rCRSd", "rCRSd_3k"]

    if not param.bowtie2.params:
        bt2_options = " --very-sensitive"
        if args.paired_end:
            bt2_options += " -X 2000"
    else:
        bt2_options = param.bowtie2.params

    # samtools sort needs a temporary directory
    tempdir = tempfile.mkdtemp(dir=map_genome_folder)
    os.chmod(tempdir, 0o771)
    pm.clean_add(tempdir)

    # check input for zipped or not
    unmap_fq1_gz = unmap_fq1 + ".gz"
    unmap_fq2_gz = unmap_fq2 + ".gz"

    bt2_index = rgc.seek(args.genome_assembly, BT2_IDX_KEY)
    if not bt2_index.endswith(args.genome_assembly):
        bt2_index = os.path.join(
            rgc.seek(args.genome_assembly, BT2_IDX_KEY),
                     args.genome_assembly)

    if _itsa_file(unmap_fq1_gz) and not _itsa_file(unmap_fq1):
        cmd = (ngstk.ziptool + " -d " + unmap_fq1_gz)
        pm.run(cmd, mapping_genome_bam)
        to_compress.append(unmap_fq1)
    if args.paired_end:
        if _itsa_file(unmap_fq2_gz) and not _itsa_file(unmap_fq2):
            cmd = (ngstk.ziptool + " -d " + unmap_fq2_gz)
            pm.run(cmd, mapping_genome_bam)
        to_compress.append(unmap_fq2)

    cmd = tools.bowtie2 + " -p " + str(pm.cores)
    cmd += bt2_options
    cmd += " --rg-id " + args.sample_name
    cmd += " -x " + bt2_index
    if args.paired_end:
        cmd += " --rf -1 " + unmap_fq1 + " -2 " + unmap_fq2
    else:
        cmd += " -U " + unmap_fq1
    cmd += " | " + tools.samtools + " view -bS - -@ 1 "
    cmd += " | " + tools.samtools + " sort - -@ 1"
    cmd += " -T " + tempdir
    cmd += " -o " + mapping_genome_bam_temp

    if not args.complexity and int(args.umi_len) > 0:
        # check input for zipped or not
        if pypiper.is_gzipped_fastq(unmap_fq1_dups):
            cmd = (ngstk.ziptool + " -d " + (unmap_fq1_dups + ".gz"))
            pm.run(cmd, mapping_genome_bam)
        if args.paired_end:
            if pypiper.is_gzipped_fastq(unmap_fq2_dups):
                cmd = (ngstk.ziptool + " -d " + (unmap_fq2_dups + ".gz"))
                pm.run(cmd, mapping_genome_bam)

        cmd_dups = tools.bowtie2 + " -p " + str(pm.cores)
        cmd_dups += bt2_options
        cmd_dups += " --rg-id " + args.sample_name
        cmd_dups += " -x " + bt2_index
        if args.paired_end:
            cmd_dups += " --rf -1 " + unmap_fq1_dups + " -2 " + unmap_fq2_dups
        else:
            cmd_dups += " -U " + unmap_fq1_dups
        cmd_dups += " | " + tools.samtools + " view -bS - -@ 1 "
        cmd_dups += " | " + tools.samtools + " sort - -@ 1"
        cmd_dups += " -T " + tempdir
        cmd_dups += " -o " + mapping_genome_bam_temp_dups

    # Split genome mapping result bamfile into two: high-quality aligned
    # reads (keepers) and unmapped reads (in case we want to analyze the
    # altogether unmapped reads)
    # -q 10: skip alignments with MAPQ less than 10
    cmd2 = (tools.samtools + " view -q 10 -b -@ " + str(pm.cores) +
            " -U " + failQC_genome_bam + " ")
    cmd2 += mapping_genome_bam_temp + " > " + mapping_genome_bam

    if not args.complexity and int(args.umi_len) > 0:
        cmd2_dups = (tools.samtools + " view -q 10 -b -@ " + str(pm.cores) +
                     " -U " + failQC_genome_bam_dups + " ")
        cmd2_dups += mapping_genome_bam_temp_dups + " > " + mapping_genome_bam_dups
        pm.clean_add(failQC_genome_bam_dups)

    def check_alignment_genome(temp_bam, bam):
        mr = ngstk.count_mapped_reads(temp_bam, args.paired_end)
        ar = ngstk.count_mapped_reads(bam, args.paired_end)

        if float(ar) < 1:
            err_msg = "No aligned reads. Check alignment settings."
            pm.fail_pipeline(RuntimeError(err_msg))
        if args.paired_end:
            ar = float(ar)/2


        tmp = pm.pipestat.retrieve(None,"Raw_reads")
        if tmp:
            rr = float(tmp)
        else:
            rr = 0

        tmp = pm.pipestat.retrieve(None,"Trimmed_reads_R1")
        if tmp:
            tr = float(tmp)
        else:
            tr = 0

        if os.path.exists(res.refgene_pre_mRNA):
            cmd = (tools.samtools + " depth -b " +
                   res.refgene_pre_mRNA + " " + bam +
                   " | awk '{counter++;sum+=$3}END{print sum/counter}'")
            rd = pm.checkprint(cmd)
        else:
            cmd = (tools.samtools + " depth " + bam +
                   " | awk '{counter++;sum+=$3}END{print sum/counter}'")
            rd = pm.checkprint(cmd)

        # pm.report_result("Mapped_reads", mr)
        # pm.report_result("QC_filtered_reads", round(float(mr)) - round(float(ar)))
        # pm.report_result("Aligned_reads", ar)
        # pm.report_result("Alignment_rate", round(float(ar) * 100 / float(tr), 2))
        # pm.report_result("Total_efficiency", round(float(ar) * 100 / float(rr), 2))

        pm.pipestat.report(values={
            "Mapped_reads": mr,
            "QC_filtered_reads": round(float(mr)) - round(float(ar)),
            "Aligned_reads": ar,
            "Alignment_rate": round(float(ar) * 100 / float(tr), 2),
            "Total_efficiency": round(float(ar) * 100 / float(rr), 2)
        })

        if rd and rd.strip():
            pm.pipestat.report(values={"Read_depth": round(float(rd), 2)})
            # pm.report_result("Read_depth", round(float(rd), 2))

    pm.run([cmd, cmd2], mapping_genome_bam,
           follow=lambda: check_alignment_genome(mapping_genome_bam_temp,
                                                 mapping_genome_bam))

    if not args.complexity and int(args.umi_len) > 0:
        if not _itsa_file(preseq_pdf) or args.new_start:
            pm.run([cmd_dups, cmd2_dups], mapping_genome_bam_dups)

    pm.timestamp("### Compress all unmapped read files")
    for unmapped_fq in list(set(to_compress)):
        # Compress unmapped fastq reads
        if not pypiper.is_gzipped_fastq(unmapped_fq) and not unmapped_fq == '':
            if 'unmap_dups' in unmapped_fq:
                pm.clean_add(unmapped_fq)
            else:
                cmd = (ngstk.ziptool + " " + unmapped_fq)
                unmapped_fq = unmapped_fq + ".gz"
                pm.run(cmd, unmapped_fq)

    if not args.prealignments and os.path.exists(mapping_genome_bam_temp):
        # Index the temporary bam file
        cmd = tools.samtools + " index " + mapping_genome_bam_temp
        pm.run(cmd, temp_mapping_index)
        pm.clean_add(temp_mapping_index)

        if not args.complexity and int(args.umi_len) > 0:
            cmd_dups = tools.samtools + " index " + mapping_genome_bam_temp_dups
            pm.run(cmd_dups, temp_mapping_index_dups)
            pm.clean_add(temp_mapping_index_dups)
            pm.clean_add(mapping_genome_bam_temp_dups)

    if not pm.pipestat.retrieve(None,"Mitochondrial_reads") or args.new_start:
        # Determine mitochondrial read counts
        if os.path.exists(mapping_genome_bam_temp):
            if not os.path.exists(temp_mapping_index):
                cmd = tools.samtools + " index " + mapping_genome_bam_temp
                pm.run(cmd, temp_mapping_index)
                pm.clean_add(temp_mapping_index)

            cmd = (tools.samtools + " idxstats " +
                   mapping_genome_bam_temp + " | grep")
            for name in mito_name:
                cmd += " -we '" + name + "'"
            cmd += "| cut -f 3"
            mr = pm.checkprint(cmd)
        
            # If there are mitochondrial reads, report and remove them
            if mr and mr.strip():
                pm.pipestat.report(values={"Mitochondrial_reads": round(float(mr))})
                # pm.report_result("Mitochondrial_reads", round(float(mr)))
                # Index the sort'ed BAM file first
                mapping_genome_index = os.path.join(mapping_genome_bam + ".bai")
                noMT_mapping_genome_bam = os.path.join(
                    map_genome_folder, args.sample_name + "_noMT.bam")
                chr_bed = os.path.join(map_genome_folder, "chr_sizes.bed")

                cmd1 = tools.samtools + " index " + mapping_genome_bam
                cmd2 = (tools.samtools + " idxstats " + mapping_genome_bam +
                        " | cut -f 1-2 | awk '{print $1, 0, $2}' | grep")
                for name in mito_name:
                    cmd2 += " -vwe '" + name + "'"
                cmd2 += (" > " + chr_bed)
                cmd3 = (tools.samtools + " view -L " + chr_bed + " -b -@ " +
                        str(pm.cores) + " " + mapping_genome_bam + " > " +
                        noMT_mapping_genome_bam)
                cmd4 = ("mv " + noMT_mapping_genome_bam +
                        " " + mapping_genome_bam)
                cmd5 = tools.samtools + " index " + mapping_genome_bam
                pm.run([cmd1, cmd2, cmd3, cmd4, cmd5], noMT_mapping_genome_bam)
                pm.clean_add(mapping_genome_index)
                pm.clean_add(chr_bed)

    # Remove PE2 reads
    if args.paired_end:
        pm.timestamp("### Split BAM file")
        mapping_pe1_bam = os.path.join(
            map_genome_folder, args.sample_name + "_PE1.bam")
        mapping_pe2_bam = os.path.join(
            map_genome_folder, args.sample_name + "_PE2.bam")
        cmd1 = (tools.samtools + " view -b -f 64 " + mapping_genome_bam +
                " | " + tools.samtools + " sort - -@ " + str(pm.cores) +
                " > " + mapping_pe1_bam)
        cmd2 = (tools.samtools + " view -b -f 128 " + mapping_genome_bam +
                " | " + tools.samtools + " sort - -@ " + str(pm.cores) +
                " > " + mapping_pe2_bam)
        pm.run([cmd1, cmd2], [mapping_pe1_bam, mapping_pe2_bam])
        mapping_genome_bam = mapping_pe1_bam

    ############################################################################
    #       Determine maximum read length and add seqOutBias resource          #
    ############################################################################

    if not pm.pipestat.retrieve(None,"Maximum_read_length") or args.new_start:
        if int(args.max_len) > 0:
            max_len = int(args.max_len)
        elif _itsa_file(mapping_genome_bam):
            cmd = (tools.samtools + " stats " + mapping_genome_bam +
                   " | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-")
            max_len = int(pm.checkprint(cmd))
        else:
            max_len = int(DEFAULT_MAX_LEN)
        # pm.report_result("Maximum_read_length", max_len)
        pm.pipestat.report(values={"Maximum_read_length": max_len})
    else:
        max_len = int(pm.pipestat.retrieve(None,"Maximum_read_length"))

    # At this point we can check for seqOutBias required indicies.
    # Can't do it earlier because we haven't determined the read_length of 
    # interest for mappability purposes.
    if args.sob:
        pm.debug("max_len: {}".format(max_len))  # DEBUG
        if not args.search_file:
            if max_len == DEFAULT_MAX_LEN:
                search_asset = [{"asset_name":"tallymer_index",
                                 "seek_key":"search_file",
                                 "tag_name":"default",
                                 "arg":"search_file",
                                 "user_arg":"search-file",
                                 "required":True}]
            else:
                search_asset = [{"asset_name":"tallymer_index",
                                 "seek_key":"search_file",
                                 "tag_name":max_len,
                                 "arg":"search_file",
                                 "user_arg":"search-file",
                                 "required":True}]
        elif ((args.search_file) and os.path.isfile(args.search_file) and
                os.stat(args.search_file).st_size > 0):
            res.search_file = args.search_file
        res, rgc = _add_resources(args, res, search_asset)

    # Calculate size of genome
    if not pm.pipestat.retrieve(None,"Genome_size") or args.new_start:
        genome_size = int(pm.checkprint(
            ("awk '{sum+=$2} END {printf \"%.0f\", sum}' " +
             res.chrom_sizes)))
        # pm.report_result("Genome_size", genome_size)
        pm.pipestat.report(values={"Genome_size": genome_size})
    else:
        genome_size = int(pm.pipestat.retrieve(None,"Genome_size"))

    ############################################################################
    #                     Calculate library complexity                         #
    ############################################################################
    preseq_output = os.path.join(
        QC_folder, args.sample_name + "_preseq_out.txt")
    preseq_yield = os.path.join(
        QC_folder, args.sample_name + "_preseq_yield.txt")
    preseq_counts = os.path.join(
        QC_folder, args.sample_name + "_preseq_counts.txt")
    preseq_plot = os.path.join(
        QC_folder, args.sample_name + "_preseq_plot")
    preseq_png = os.path.join(
        QC_folder, args.sample_name + "_preseq_plot.png")

    if not _itsa_file(preseq_pdf) or args.new_start:
        if not args.complexity and int(args.umi_len) > 0:
            if os.path.exists(mapping_genome_bam_temp_dups):
                if not os.path.exists(temp_mapping_index_dups):
                    cmd = tools.samtools + " index " + mapping_genome_bam_temp_dups
                    pm.run(cmd, temp_mapping_index_dups)
                    pm.clean_add(temp_mapping_index_dups)

                cmd_dups = (tools.samtools + " idxstats " +
                            mapping_genome_bam_temp_dups + " | grep")
                for name in mito_name:
                    cmd_dups += " -we '" + name + "'"
                cmd_dups += "| cut -f 3"
                mr_dups = pm.checkprint(cmd_dups)

                if mr_dups and mr_dups.strip():
                    # Index the sort'ed BAM file first
                    mapping_genome_index_dups = os.path.join(
                        mapping_genome_bam_dups + ".bai")
                    noMT_mapping_genome_bam_dups = os.path.join(
                        map_genome_folder, args.sample_name + "_noMT_dups.bam")

                    cmd1 = tools.samtools + " index " + mapping_genome_bam_dups
                    cmd2 = (tools.samtools + " idxstats " +
                            mapping_genome_bam_dups + " | cut -f 1 | grep")
                    for name in mito_name:
                        cmd2 += " -vwe '" + name + "'"
                    cmd2 += ("| xargs " + tools.samtools + " view -b -@ " +
                             str(pm.cores) + " " + mapping_genome_bam_dups +
                             " > " + noMT_mapping_genome_bam_dups)
                    cmd3 = ("mv " + noMT_mapping_genome_bam_dups + " " +
                            mapping_genome_bam_dups)
                    cmd4 = tools.samtools + " index " + mapping_genome_bam_dups
                    pm.run([cmd1, cmd2, cmd3, cmd4], mapping_genome_bam_dups)
                    pm.clean_add(mapping_genome_index_dups)

            # Remove PE2 reads
            if args.paired_end:
                dups_pe1_bam = os.path.join(
                    map_genome_folder, args.sample_name + "_dups_PE1.bam")
                dups_pe2_bam = os.path.join(
                    map_genome_folder, args.sample_name + "_dups_PE2.bam")
                cmd1 = (tools.samtools + " view -b -f 64 " +
                    mapping_genome_bam_dups + " | " + tools.samtools +
                    " sort - -@ " + str(pm.cores) + " > " + dups_pe1_bam)
                cmd2 = (tools.samtools + " view -b -f 128 " +
                    mapping_genome_bam_dups + " | " + tools.samtools +
                    " sort - -@ " + str(pm.cores) + " > " + dups_pe2_bam)
                pm.run([cmd1, cmd2], [dups_pe1_bam, dups_pe2_bam])
                mapping_genome_bam_dups = dups_pe1_bam

            pm.timestamp("### Calculate library complexity")

            cmd1 = (tools.preseq + " c_curve -v -o " + preseq_output +
                    " -B " + mapping_genome_bam_dups)
            pm.run(cmd1, preseq_output)

            cmd2 = (tools.preseq + " lc_extrap -v -o " + preseq_yield +
                    " -B " + mapping_genome_bam_dups)
            pm.run(cmd2, preseq_yield, nofail=True)

            if os.path.exists(preseq_yield):
                # cmd3 = ("bam2mr " + mapping_genome_bam_dups +
                #         " > " + preseq_mr)
                # cmd4 = (tools.preseq + " gc_extrap -v -o " + preseq_cov +
                #         " " + preseq_mr)
                cmd5 = ("echo '" + preseq_yield +
                        " '$(" + tools.samtools + " view -c -F 4 " + 
                        mapping_genome_bam_dups + ")" + "' '" +
                        "$(" + tools.samtools + " view -c -F 4 " +
                        mapping_genome_bam + ") > " + preseq_counts)

                # pm.run([cmd3, cmd4, cmd5],
                #        [preseq_mr, preseq_cov, preseq_counts])
                pm.run(cmd5, preseq_counts)
                #pm.clean_add(preseq_mr)
                pm.clean_add(mapping_genome_bam_dups)
                pm.clean_add(mapping_genome_bam_temp_dups)

                cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") +
                       " preseq " + "-i " + preseq_yield)
                if args.coverage:
                    cmd += (" -c " + str(genome_size) + " -l " + max_len)
                cmd += (" -r " + preseq_counts + " -o " + preseq_plot)

                pm.run(cmd, [preseq_pdf, preseq_png])

                # pm.report_object("Library complexity", preseq_pdf,
                #                  anchor_image=preseq_png)
                pm.pipestat.report(values={"Library_complexity": {"path": preseq_pdf, "thumbnail_path": preseq_png, "title": "Library complexity"}})
                if not pm.pipestat.retrieve(None,'Frac_exp_unique_at_10M') or args.new_start:
                    # Report the expected unique at 10M reads
                    cmd = ("grep -w '10000000' " + preseq_yield +
                           " | awk '{print $2}'")
                    expected_unique = pm.checkprint(cmd)
                    if expected_unique:
                        fraction_unique = float(expected_unique)/float(10000000)
                        # pm.report_result("Frac_exp_unique_at_10M", round(fraction_unique, 4))
                        pm.pipestat.report(values={"Frac_exp_unique_at_10M": round(fraction_unique, 4)})
            else:
                print("Unable to calculate library complexity.")

    ############################################################################
    #         Calculate quality control metrics for the alignment file         #
    ############################################################################
    pm.timestamp("### Calculate NRF, PBC1, and PBC2")

    # Need index for mapping_genome_bam before calculating bamQC metrics
    mapping_genome_index = os.path.join(mapping_genome_bam + ".bai")
    cmd = tools.samtools + " index " + mapping_genome_bam
    pm.run(cmd, mapping_genome_index)

    bamQC = os.path.join(QC_folder, args.sample_name + "_bamQC.tsv")
    cmd = tool_path("bamQC.py")
    cmd += " --silent"
    cmd += " -i " + mapping_genome_bam
    cmd += " -c " + str(pm.cores)
    cmd += " -o " + bamQC

    def report_bam_qc(bamqc_log):
        # Reported BAM QC metrics via the bamQC metrics file
        if os.path.isfile(bamqc_log):
            cmd1 = ("awk '{ for (i=1; i<=NF; ++i) {" +
                    " if ($i ~ \"NRF\") c=i } getline; print $c }' " +
                    bamqc_log)
            cmd2 = ("awk '{ for (i=1; i<=NF; ++i) {" +
                    " if ($i ~ \"PBC1\") c=i } getline; print $c }' " +
                    bamqc_log)
            cmd3 = ("awk '{ for (i=1; i<=NF; ++i) {" +
                    " if ($i ~ \"PBC2\") c=i } getline; print $c }' " +
                    bamqc_log)
            nrf = pm.checkprint(cmd1)
            pbc1 = pm.checkprint(cmd2)
            pbc2 = pm.checkprint(cmd3)
        else:
            # there were no successful chromosomes yielding results
            nrf = 0
            pbc1 = 0
            pbc2 = 0

        # pm.report_result("NRF", round(float(nrf),2))
        # pm.report_result("PBC1", round(float(pbc1),2))
        # pm.report_result("PBC2", round(float(pbc2), 2))

        pm.pipestat.report(values={
            "NRF": round(float(nrf),2),
            "PBC1": round(float(pbc1),2),
            "PBC2": round(float(pbc2), 2),
        })

    pm.run(cmd, bamQC, follow=lambda: report_bam_qc(bamQC))

    ############################################################################
    #                     Produce unmapped reads file                          #
    ############################################################################
    def count_unmapped_reads():
        # Report total number of unmapped reads (-f 4)
        cmd = (tools.samtools + " view -c -f 4 -@ " + str(pm.cores) +
               " " + mapping_genome_bam_temp)
        ur = pm.checkprint(cmd)
        pm.report_result("Unmapped_reads", round(float(ur)))

    unmap_cmd = tools.samtools + " view -b -@ " + str(pm.cores)
    if args.paired_end:
        # require both read and mate unmapped
        unmap_cmd += " -f 12 "
    else:
        # require only read unmapped
        unmap_cmd += " -f 4 "

    unmap_cmd += " " + mapping_genome_bam_temp + " > " + unmap_genome_bam
    pm.run(unmap_cmd, unmap_genome_bam, follow=count_unmapped_reads)

    # Remove temporary bam file from unmapped file production
    pm.clean_add(mapping_genome_bam_temp)

    ############################################################################
    #                          Separate BAM by strand                          #
    ############################################################################
    pm.timestamp("### Split BAM by strand")
    plus_bam = os.path.join(
        map_genome_folder, args.sample_name + "_plus.bam")
    minus_bam = os.path.join(
        map_genome_folder, args.sample_name + "_minus.bam")
    
    cmd1 = build_command([
        tools.samtools,
        "view",
        "-bh",
        ("-F", 20),
        mapping_genome_bam,
        (">", plus_bam)
    ])

    cmd2 = build_command([
        tools.samtools,
        "view",
        "-bh",
        ("-f", 16),
        mapping_genome_bam,
        (">", minus_bam)
    ])
    
    pm.run([cmd1, cmd2], [plus_bam, minus_bam])

    ############################################################################
    #                             TSS enrichment                               #
    ############################################################################
    if not os.path.exists(res.refgene_tss):
        print("Skipping TSS -- TSS enrichment requires TSS annotation file: {}"
              .format(res.refgene_tss))
    else:
        pm.timestamp("### Calculate TSS enrichment")

        Tss_plus = os.path.join(QC_folder, args.sample_name +
                                "_plus_TssEnrichment.txt")
        Tss_minus = os.path.join(QC_folder, args.sample_name +
                                 "_minus_TssEnrichment.txt")
        if not pm.pipestat.retrieve(None,"TSS_non-coding_score") or args.new_start:
            # Split TSS file
            plus_TSS  = os.path.join(QC_folder, "plus_TSS.tsv")
            minus_TSS = os.path.join(QC_folder, "minus_TSS.tsv")
            cmd = ("sed -n -e '/[[:space:]]+/w " +
                   plus_TSS + "' -e '/[[:space:]]-/w " +
                   minus_TSS + "' " + res.refgene_tss)
            pm.run(cmd, [plus_TSS, minus_TSS])

            # pyTssEnrichment requires indexed bam
            if not os.path.exists(mapping_genome_index):
                cmd = build_command([
                    tools.samtools,
                    "index",
                    mapping_genome_bam
                ])
                pm.run(cmd, mapping_genome_index)

            # Plus TSS enrichment
            cmd = tool_path("pyTssEnrichment.py")
            cmd += " -a " + mapping_genome_bam + " -b " + plus_TSS + " -p ends"
            cmd += " -c " + str(pm.cores)
            cmd += " -z -v -s 6 -o " + Tss_plus
            pm.run(cmd, Tss_plus, nofail=True)
            pm.clean_add(plus_TSS)
            pm.clean_add(Tss_plus)

            with open(Tss_plus) as f:
                floats = list(map(float, f))
            try:
                # Catch if the TSS score is trying to divide by 0           
                list_len = 0.05*float(len(floats))
                normTSS = [x / (sum(floats[1:int(list_len)]) /
                           len(floats[1:int(list_len)])) for x in floats]
                max_index = normTSS.index(max(normTSS))

                if (((normTSS[max_index]/normTSS[max_index-1]) > 1.5) and
                    ((normTSS[max_index]/normTSS[max_index+1]) > 1.5)):
                    tmpTSS = list(normTSS)
                    del tmpTSS[max_index]
                    max_index = tmpTSS.index(max(tmpTSS)) + 1

                Tss_score = round(
                    (sum(normTSS[int(max_index-50):int(max_index+50)])) /
                    (len(normTSS[int(max_index-50):int(max_index+50)])), 1)

                # pm.report_result("TSS_coding_score", round(Tss_score, 1))
                pm.pipestat.report(values={"TSS_coding_score": round(Tss_score, 1)})
            except ZeroDivisionError:
                # pm.report_result("TSS_coding_score", 0)
                pm.pipestat.report(values={"TSS_coding_score": 0})
                pass

            # Minus TSS enrichment
            cmd = tool_path("pyTssEnrichment.py")
            cmd += " -a " + mapping_genome_bam + " -b " + minus_TSS + " -p ends"
            cmd += " -c " + str(pm.cores)
            cmd += " -z -v -s 6 -o " + Tss_minus
            pm.run(cmd, Tss_minus, nofail=True)
            pm.clean_add(minus_TSS)
            pm.clean_add(Tss_minus)

            with open(Tss_minus) as f:
                floats = list(map(float, f))
            try:
                # If the TSS enrichment is 0, don't report
                list_len = 0.05*float(len(floats))
                normTSS = [x / (sum(floats[1:int(list_len)]) /
                           len(floats[1:int(list_len)])) for x in floats]
                max_index = normTSS.index(max(normTSS))

                if (((normTSS[max_index]/normTSS[max_index-1]) > 1.5) and
                    ((normTSS[max_index]/normTSS[max_index+1]) > 1.5)):
                    tmpTSS = list(normTSS)
                    del tmpTSS[max_index]
                    max_index = tmpTSS.index(max(tmpTSS)) + 1

                Tss_score = round(
                    (sum(normTSS[int(max_index-50):int(max_index+50)])) /
                    (len(normTSS[int(max_index-50):int(max_index+50)])), 1)

                # pm.report_result("TSS_non-coding_score", round(Tss_score, 1))
                pm.pipestat.report(values={"TSS_non-coding_score": round(Tss_score, 1)})
            except ZeroDivisionError:
                # pm.report_result("TSS_non-coding_score", 0)
                pm.pipestat.report(values={"TSS_non-coding_score": 0})
                pass

        # Call Rscript to plot TSS Enrichment
        TSS_pdf = os.path.join(QC_folder,  args.sample_name +
                               "_TSSenrichment.pdf")
        cmd = (tools.Rscript + " " + tool_path("PEPPRO.R"))
        cmd += " tss -i " + Tss_plus + " " + Tss_minus
        pm.run(cmd, TSS_pdf, nofail=True)

        TSS_png = os.path.join(QC_folder,  args.sample_name +
                               "_TSSenrichment.png")
        # pm.report_object("TSS enrichment", TSS_pdf, anchor_image=TSS_png)
        pm.pipestat.report(values={"TSS_enrichment": {"path": TSS_pdf, "thumbnail_path": TSS_png, "title": "TSS enrichment"}})

    ############################################################################
    #                Pause index calculation and plotting                      #
    ############################################################################
    chr_order = os.path.join(QC_folder, "chr_order.txt")
    chr_keep = os.path.join(QC_folder, "chr_keep.txt")

    if not os.path.exists(chr_order):
        cmd = (tools.samtools + " view -H " + mapping_genome_bam +
               " | grep 'SN:' | awk -F':' '{print $2,$3}' | " +
               "awk -F' ' -v OFS='\t' '{print $1,$3}' > " + chr_order)
        pm.run(cmd, chr_order)
        pm.clean_add(chr_order)

    if not os.path.exists(chr_keep):
        if os.path.exists(chr_order):
            cmd = ("cut -f 1 " + chr_order + " > " + chr_keep)
            pm.run(cmd, chr_keep)
            pm.clean_add(chr_keep)
        else:
            cmd1 = (tools.samtools + " view -H " + mapping_genome_bam +
                    " | grep 'SN:' | awk -F':' '{print $2,$3}' | " +
                    "awk -F' ' -v OFS='\t' '{print $1,$3}' > " + chr_order)
            cmd2 = ("cut -f 1 " + chr_order + " > " + chr_keep)
            pm.run([cmd1, cmd2], [chr_order, chr_keep])
            pm.clean_add(chr_order)
            pm.clean_add(chr_keep)

    if not os.path.exists(res.ensembl_tss):
        if not os.path.exists(res.ensembl_gene_body):
            print("Skipping PI -- Pause index requires 'TSS' and 'gene body' annotation files: {} and {}"
                  .format(res.ensembl_tss, res.ensembl_gene_body))
        else:
            print("Skipping PI -- Pause index requires 'TSS' annotation file: {}"
                  .format(res.ensembl_tss))
    elif not os.path.exists(res.ensembl_gene_body):
        print("Skipping PI -- Pause index requires 'gene body' annotation file: {}"
              .format(res.ensembl_gene_body))
    else:
        pm.timestamp("### Calculate Pause Index (PI)")
        temp = tempfile.NamedTemporaryFile(dir=QC_folder, delete=False)
        pause_index = os.path.join(QC_folder, args.sample_name +
                                   "_pause_index.bed")
        pause_index_gz = os.path.join(QC_folder, args.sample_name +
                                      "_pause_index.bed.gz")
        if not pm.pipestat.retrieve(None,"Pause_index") or args.new_start:
            # Remove missing chr from PI annotations
            tss_local = os.path.join(QC_folder,
                args.genome_assembly + "_ensembl_tss.bed")
            body_local = os.path.join(QC_folder,
                args.genome_assembly + "_ensembl_gene_body.bed")
            cmd1 = ("grep -wf " + chr_keep + " " + res.ensembl_tss + " | " +
                    tools.bedtools + " sort -i stdin -faidx " + chr_order + 
                    " > " + tss_local)
            cmd2 = ("grep -wf " + chr_keep + " " + res.ensembl_gene_body +
                    " | " + tools.bedtools + " sort -i stdin -faidx " +
                    chr_order + " > " + body_local)
            pm.run([cmd1,cmd2], [tss_local, body_local], nofail=True)
            pm.clean_add(tss_local)
            pm.clean_add(body_local)

            # Determine coverage of highest scoring TSS
            TSS_density = os.path.join(QC_folder, args.sample_name +
                                       "_TSS_density.bed")
            cmd = (tools.bedtools + " coverage -sorted -counts -s -a " +
                   tss_local + " -b " + mapping_genome_bam +
                   " -g " + chr_order + " | awk '$7>0' | " + 
                   "sort -k4,4 -k7,7nr | " +
                   "sort -k4,4 -u > " + TSS_density)
            pm.run(cmd, TSS_density, nofail=True)
            pm.clean_add(TSS_density)

            # Determine coverage of gene body
            body_density = os.path.join(QC_folder, args.sample_name +
                                        "_gene_body_density.bed")
            cmd = (tools.bedtools + " coverage -sorted -counts -s -a " +
                   body_local + " -b " + mapping_genome_bam +
                   " -g " + chr_order + " | awk '$7>0' | " + 
                   "sort -k4 > " + body_density)
            pm.run(cmd, body_density, nofail=True)
            pm.clean_add(body_density)

            # Calculate expression and pause indicies
            cmd = ("join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7" +
                   " 2.2 2.3 2.7 " + TSS_density + " " + body_density +
                   " | awk -v OFS='\t' '{ if ($5 == \"+\"){print $1, $2, $8," +
                   " $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), " +
                   "($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} " +
                   "else {print $1, $2, $8, $4, " +
                   "sqrt((($6+$9)/sqrt(($3-$7)^2))^2)," +
                   "($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' " +
                   "| env LC_COLLATE=C sort -k1,1 -k2,2n > " + temp.name)
            pm.run(cmd, pause_index, nofail=True)
            temp.close()

            # Calculate median expression
            cmd = ("awk '{print $5}' " + temp.name + " | sort -n | " +
                   "awk 'BEGIN{i=0} {s[i]=$1; i++;} " + 
                   "END{print s[int(NR*0.5-0.5)]}'")
            median_expr = pm.checkprint(cmd)
            try:
                cmd = ("awk -v OFS='\t' '{ if ($5 > " + str(median_expr) +
                       ") " + "{print $1, $2, $3, $4, $6, $7}}' " + 
                       temp.name + " > " + pause_index)
                pm.run(cmd, pause_index, nofail=True)
            except ZeroDivisionError:
                # Fall back to using all data to determine pause index
                cmd = ("join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 " +
                       "2.2 2.3 2.7 " + TSS_density + " " + body_density +
                       " | awk -v OFS='\t' '{print $1, $2, $3, $4, " +
                       "($6/($3-$2))/($9/($8-$7)), $5}' " +
                       "| env LC_COLLATE=C sort -k1,1 -k2,2n > " + pause_index)
                pm.run(cmd, pause_index, nofail=True)
                pass

            # Median pause index
            cmd = ("sort -k5,5n " + pause_index +
                   " | awk ' { a[i++]=$5; } END " + 
                   "{ x=int((i+1)/2); if (x < (i+1)/2) " +
                   "print (a[x-1]+a[x])/2; else print a[x-1]; }'")
            val = pm.checkprint(cmd)
            if val and val.strip():
                pi = float(val)
                pm.pipestat.report(values={"Pause_index": round(pi, 2)})
                # pm.report_result("Pause_index", round(pi, 2))

        # Plot pause index distribution
        pi_pdf = os.path.join(QC_folder, args.sample_name +
                              "_pause_index.pdf")
        pi_png = os.path.join(QC_folder, args.sample_name +
                              "_pause_index.png")

        if _itsa_file(pause_index_gz) and not _itsa_file(pi_pdf):
            cmd = (ngstk.ziptool + " -d " + pause_index_gz)
            pm.run(cmd, pause_index)

        cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") + 
               " pi --annotate -i " + pause_index)
        pm.run(cmd, pi_pdf, nofail=True)
        # pm.report_object("Pause index", pi_pdf, anchor_image=pi_png)
        # TODO: had to changed key since it was masked by prev stat
        pm.pipestat.report(values={"Pause_index_plot": {"path": pi_pdf, "thumbnail_path": pi_png, "title": "TSS Pause_index_plot"}})

        if not is_gzipped(pause_index):
            cmd = (ngstk.ziptool + " -f " + pause_index)
            pause_index = pause_index + ".gz"
            pm.run(cmd, pause_index)

    ############################################################################
    #           Calculate Fraction of Reads in Pre-mRNA (FRiP)                 #
    ############################################################################
    signal_folder = os.path.join(
        param.outfolder, "signal_" + args.genome_assembly)
    ngstk.make_dir(signal_folder)

    if not os.path.exists(res.refgene_pre_mRNA):
        print("Skipping FRiP and gene coverage calculation which require the "
              "pre-mRNA annotation file: {}"
              .format(res.refgene_pre_mRNA))
    else:
        pm.timestamp("### Calculate Fraction of Reads in pre-mature mRNA")
        if not pm.pipestat.retrieve(None,'Plus_FRiP') or args.new_start:
            # Plus
            plus_frip = calc_frip(plus_bam, res.refgene_pre_mRNA,
                                  frip_func=ngstk.simple_frip,
                                  pipeline_manager=pm)
            pm.pipestat.report(values={"Plus_FRiP": round(plus_frip, 2)})
            # pm.report_result("Plus_FRiP", round(plus_frip, 2))

        if not pm.pipestat.retrieve(None,'Minus_FRiP') or args.new_start:
            # Minus
            minus_frip = calc_frip(minus_bam, res.refgene_pre_mRNA,
                                   frip_func=ngstk.simple_frip,
                                   pipeline_manager=pm)
            # pm.report_result("Minus_FRiP", round(minus_frip, 2))
            pm.pipestat.report(values={"Minus_FRiP": round(minus_frip, 2)})

        # Calculate gene coverage
        gene_cov = os.path.join(signal_folder,
                                args.sample_name + "_gene_coverage.bed")
        gene_sort = os.path.join(QC_folder, args.genome_assembly +
                                 "_gene_sort.bed")
        cmd1 = ("grep -wf " + chr_keep + " " + res.refgene_pre_mRNA +
                " | " + tools.bedtools + " sort -i stdin -faidx " +
                chr_order + " > " + gene_sort)
        cmd2 = (tools.bedtools + " coverage -sorted -counts -s -a " +
                gene_sort + " -b " + mapping_genome_bam +
                " -g " + chr_order + " > " + gene_cov)
        pm.run([cmd1, cmd2], gene_cov)
        pm.clean_add(gene_sort)

    ############################################################################
    #                        Extract genomic features                          #
    ############################################################################
    # Generate local unzipped annotation file
    anno_local = os.path.join(raw_folder,
                              args.genome_assembly + "_annotations.bed")
    anno_zip = os.path.join(raw_folder,
                            args.genome_assembly + "_annotations.bed.gz")

    if (not os.path.exists(anno_local) and
        not os.path.exists(anno_zip) and
        os.path.exists(res.feat_annotation) or
        args.new_start):

        if res.feat_annotation.endswith(".gz"):
            cmd1 = ("ln -sf " + res.feat_annotation + " " + anno_zip)
            cmd2 = (ngstk.ziptool + " -d -c " + anno_zip +
                    " > " + anno_local)
            pm.run([cmd1, cmd2], anno_local)
            pm.clean_add(anno_local)
        elif res.feat_annotation.endswith(".bed"):
            cmd = ("ln -sf " + res.feat_annotation + " " + anno_local)
            pm.run(cmd, anno_local)
            pm.clean_add(anno_local)
        else:
            print("Skipping read and peak annotation...")
            print("This requires a {} annotation file."
                  .format(args.genome_assembly))
            print("Could not find {}.`"
                  .format(str(os.path.dirname(res.feat_annotation))))
    elif os.path.exists(anno_zip) or args.new_start:
        cmd = (ngstk.ziptool + " -d -c " + anno_zip +
               " > " + anno_local)
        pm.run(cmd, anno_local)
        pm.clean_add(anno_local)

    ############################################################################ 
    #                  Determine genomic feature coverage                      #
    ############################################################################
    pm.timestamp("### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF)")

    # Cummulative Fraction of Reads in Features (cFRiF)
    cFRiF_PDF = os.path.join(QC_folder, args.sample_name + "_cFRiF.pdf")
    cFRiF_PNG = os.path.join(QC_folder, args.sample_name + "_cFRiF.png")

    # Fraction of Reads in Feature (FRiF)
    FRiF_PDF = os.path.join(QC_folder, args.sample_name + "_FRiF.pdf")
    FRiF_PNG = os.path.join(QC_folder, args.sample_name + "_FRiF.png")

    if not os.path.exists(cFRiF_PDF) or args.new_start:
        anno_files = list()
        anno_list_plus = list()
        anno_list_minus = list()

        if os.path.isfile(anno_local):
            # Get list of features
            if args.prioritize:
                cmd1 = ("cut -f 4 " + anno_local + " | uniq")
            else:
                cmd1 = ("cut -f 4 " + anno_local + " | sort -u")
            ft_list = pm.checkprint(cmd1, shell=True)
            ft_list = ft_list.splitlines()
            if param.precedence.params:
                p_list = param.precedence.params.split(",")
                p_list = [feature.strip() for feature in p_list]
                if all(feature in ft_list for feature in p_list):
                    ft_list = p_list
                else:
                    pm.warning("The provided precedence list ({}) of features "
                               "are not all present in your annotation file "
                               "({})".format(str(p_list), anno_local))
                    pm.warning("Defaulting to the order of features in "
                               "{}".format(anno_local))

            # Split annotation file on features
            cmd2 = ("awk -F'\t' '{print>\"" + QC_folder + "/\"$4}' " +
                    anno_local)

            if args.prioritize:
                if len(ft_list) >= 1:
                    for pos, anno in enumerate(ft_list):
                        # working files
                        anno_file = os.path.join(QC_folder, str(anno))
                        valid_name = str(re.sub('[^\w_.)( -]', '', anno).strip().replace(' ', '_'))
                        file_name = os.path.join(QC_folder, valid_name)
                        anno_sort = os.path.join(QC_folder,
                                                 valid_name + "_sort.bed")
                        anno_cov_plus = os.path.join(QC_folder,
                                                     args.sample_name + "_" +
                                                     valid_name +
                                                     "_plus_coverage.bed")
                        anno_cov_minus = os.path.join(QC_folder,
                                                      args.sample_name + "_" +
                                                      valid_name +
                                                      "_minus_coverage.bed")

                        # Extract feature files
                        pm.run(cmd2, anno_file)

                        # Rename files to valid file_names
                        # Avoid 'mv' "are the same file" error
                        if not os.path.exists(file_name):
                            cmd = 'mv "{old}" "{new}"'.format(old=anno_file,
                                                              new=file_name)
                            pm.run(cmd, file_name)

                        # Sort files (ensure only aligned chromosomes are kept)
                        # Need to cut -f 1-6 if you want strand information
                        # Not all features are stranded
                        # TODO: check for strandedness (*only works on some features)
                        cmd3 = ("cut -f 1 " + chr_order + " | grep -wf - " +
                                file_name + " | cut -f 1-3 | " +
                                "bedtools sort -i stdin -faidx " +
                                chr_order + " | bedtools merge -i stdin > " +                           
                                anno_sort)
                        # for future stranded possibilities include for merge
                        # "-c 4,5,6 -o collapse,collapse,collapse > " +
                        pm.run(cmd3, anno_sort)
                        
                        anno_files.append(anno_sort)
                        anno_list_plus.append(anno_cov_plus)
                        anno_list_minus.append(anno_cov_minus)

                        pm.clean_add(file_name)
                        pm.clean_add(anno_sort)
                        pm.clean_add(anno_cov_plus)
                        pm.clean_add(anno_cov_minus)

                    # Iteratively prioritize annotations by order presented
                    anno_files.reverse()
                    if len(anno_files) >= 1:
                        idx = list(range(0,len(anno_files)))
                        #idx.reverse()
                        file_count = 1
                        for annotation in anno_files:
                            del idx[0]
                            if file_count < len(anno_files):
                                file_count += 1
                                for i in idx:
                                    if annotation is not anno_files[i]:
                                        os.path.join(QC_folder)
                                        temp = tempfile.NamedTemporaryFile(dir=QC_folder, delete=False)
                                        #os.chmod(temp.name, 0o771)
                                        cmd1 = ("bedtools subtract -a " +
                                                annotation + " -b " +
                                                anno_files[i] + " > " + 
                                                temp.name)
                                        cmd2 = ("mv " + temp.name +
                                                " " + annotation)
                                        pm.run([cmd1, cmd2], cFRiF_PDF)
                                        temp.close()

                    anno_list_plus.reverse()
                    anno_list_minus.reverse()
                    if len(anno_files) >= 1:
                        for idx, annotation in enumerate(anno_files):
                            # Identifies unstranded coverage
                            # Would need to use '-s' flag to be stranded
                            if _itsa_file(annotation):
                                cmd4 = (tools.bedtools +
                                        " coverage -sorted -a " +
                                        annotation + " -b " + plus_bam +
                                        " -g " + chr_order + " > " +
                                        anno_list_plus[idx])
                                pm.run(cmd4, cFRiF_PDF)
                            if _itsa_file(annotation):
                                cmd5 = (tools.bedtools +
                                        " coverage -sorted -a " +
                                        annotation + " -b " + minus_bam +
                                        " -g " + chr_order + " > " +
                                        anno_list_minus[idx])
                                pm.run(cmd5, cFRiF_PDF)
            else:
                if len(ft_list) >= 1:
                    for pos, anno in enumerate(ft_list):
                        # working files
                        anno_file = os.path.join(QC_folder, str(anno))
                        valid_name = str(re.sub('[^\w_.)( -]', '', anno).strip().replace(' ', '_'))
                        file_name = os.path.join(QC_folder, valid_name)
                        anno_sort = os.path.join(QC_folder,
                                                 valid_name + "_sort.bed")
                        anno_cov_plus = os.path.join(QC_folder,
                                                     args.sample_name + "_" +
                                                     valid_name +
                                                     "_plus_coverage.bed")
                        anno_cov_minus = os.path.join(QC_folder,
                                                      args.sample_name + "_" +
                                                      valid_name +
                                                      "_minus_coverage.bed")

                        # Extract feature files
                        pm.run(cmd2, anno_file)

                        # Rename files to valid file_names
                        # Avoid 'mv' "are the same file" error
                        if not os.path.exists(file_name):
                            cmd = 'mv "{old}" "{new}"'.format(old=anno_file,
                                                              new=file_name)
                            pm.run(cmd, file_name)

                        # Sort files (ensure only aligned chromosomes are kept)
                        # Need to cut -f 1-6 if you want strand information
                        # Not all features are stranded
                        # TODO: check for strandedness
                        cmd3 = ("cut -f 1 " + chr_order + " | grep -wf - " +
                                file_name + " | cut -f 1-3 | " +
                                "bedtools sort -i stdin -faidx " +
                                chr_order + " > " + anno_sort)
                        pm.run(cmd3, anno_sort)
                        
                        anno_list_plus.append(anno_cov_plus)
                        anno_list_minus.append(anno_cov_minus)
                        # Identifies unstranded coverage
                        # Would need to use '-s' flag to be stranded
                        cmd4 = (tools.bedtools + " coverage -sorted " +
                                " -a " + anno_sort + " -b " + plus_bam +
                                " -g " + chr_order + " > " +
                                anno_cov_plus)
                        cmd5 = (tools.bedtools + " coverage -sorted " +
                                " -a " + anno_sort + " -b " + minus_bam +
                                " -g " + chr_order + " > " +
                                anno_cov_minus)
                        pm.run(cmd4, anno_cov_plus)
                        pm.run(cmd5, anno_cov_minus)

                        pm.clean_add(file_name)
                        pm.clean_add(anno_sort)
                        pm.clean_add(anno_cov_plus)
                        pm.clean_add(anno_cov_minus)

    ############################################################################
    #                            Plot cFRiF/FRiF                               #
    ############################################################################
    pm.timestamp("### Plot cFRiF/FRiF")

    # Plus
    if not os.path.exists(cFRiF_PDF) or args.new_start:
        if args.prioritize:
            # Count bases, not reads
            # return to original priority ranked order
            anno_list_plus.reverse()
            anno_list_minus.reverse()
            count_cmd = (tools.bedtools + " genomecov -ibam " + plus_bam +
                         " -bg | awk '{sum+=($3-$2)}END{print sum}'")
        else:
            # Count reads
            count_cmd = (tools.samtools + " view -@ " + str(pm.cores) + " " +
                         param.samtools.params + " -c -F4 " + plus_bam)

        plus_read_count = pm.checkprint(count_cmd)
        plus_read_count = str(plus_read_count).rstrip()

        cFRiF_cmd = [tools.Rscript, tool_path("PEPPRO.R"), "frif",
                     "-s", args.sample_name, "-z", str(genome_size).rstrip(),
                     "-n", plus_read_count, "-y", "cfrif"]

        FRiF_cmd = [tools.Rscript, tool_path("PEPPRO.R"), "frif",
                    "-s", args.sample_name, "-z", str(genome_size).rstrip(),
                    "-n", plus_read_count, "-y", "frif"]

        if not args.prioritize:
            # Use reads for calculation
            cFRiF_cmd.append("--reads")
            FRiF_cmd.append("--reads")

        cFRiF_cmd.append("-o")
        cFRiF_cmd.append(cFRiF_PDF)
        cFRiF_cmd.append("--bed")

        FRiF_cmd.append("-o")
        FRiF_cmd.append(FRiF_PDF)
        FRiF_cmd.append("--bed")

        if anno_list_plus:
            for cov in anno_list_plus:
                if _itsa_file(cov):
                    cFRiF_cmd.append(cov)
                    FRiF_cmd.append(cov)
            cmd = build_command(cFRiF_cmd)
            pm.run(cmd, cFRiF_PDF, nofail=False)
            pm.pipestat.report(values={"cFRiF": {"path": cFRiF_PDF, "thumbnail_path": cFRiF_PNG, "title": "cFRiF"}})
            # pm.report_object("cFRiF", cFRiF_PDF, anchor_image=cFRiF_PNG)

            cmd = build_command(FRiF_cmd)
            pm.run(cmd, FRiF_PDF, nofail=False)
            # pm.report_object("FRiF", FRiF_PDF, anchor_image=FRiF_PNG)
            pm.pipestat.report(values={"FRiF": {"path": FRiF_PDF, "thumbnail_path": FRiF_PNG, "title": "FRiF"}})

    ############################################################################
    #                         Report mRNA contamination                        #
    ############################################################################
    if (os.path.exists(res.refgene_exon) and
        os.path.exists(res.refgene_intron)):

        pm.timestamp("### Calculate mRNA contamination")
        intron_exon = os.path.join(QC_folder, args.sample_name +
                                   "_exon_intron_ratios.bed")
        intron_exon_gz = os.path.join(QC_folder, args.sample_name +
                                      "_exon_intron_ratios.bed.gz")

        if not pm.pipestat.retrieve(None,"mRNA_contamination") or args.new_start:
            # Sort exons and introns
            exons_sort = os.path.join(QC_folder, args.genome_assembly +
                                      "_exons_sort.bed")
            introns_sort = os.path.join(QC_folder, args.genome_assembly +
                                        "_introns_sort.bed")
            cmd1 = ("grep -wf " + chr_keep + " " + res.refgene_exon +
                    " | " + tools.bedtools + " sort -i stdin -faidx " +
                    chr_order + " > " + exons_sort)
            # a single sort fails to sort a 1 bp different start position intron
            cmd2 = ("grep -wf " + chr_keep + " " + res.refgene_intron +
                    " | " + tools.bedtools + " sort -i stdin -faidx " +
                    chr_order + " | " + tools.bedtools +
                    " sort -i stdin -faidx " + chr_order + " > " + introns_sort)
            pm.run([cmd1, cmd2], [exons_sort, introns_sort], nofail=True)
            pm.clean_add(exons_sort)
            pm.clean_add(introns_sort)
            
            # Determine coverage of exons/introns
            exons_cov = os.path.join(QC_folder, args.sample_name +
                                     "_exons_coverage.bed")
            introns_cov = os.path.join(QC_folder, args.sample_name +
                                       "_introns_coverage.bed")
            cmd1 = (tools.bedtools + " coverage -sorted -counts -s -a " +
                    exons_sort + " -b " + mapping_genome_bam +
                    " -g " + chr_order + " > " + exons_cov)
            cmd2 = (tools.bedtools + " coverage -sorted -counts -s -a " +
                    introns_sort + " -b " + mapping_genome_bam +
                    " -g " + chr_order + " > " + introns_cov)
            pm.run([cmd1, cmd2], [exons_cov, introns_cov], nofail=True)
            pm.clean_add(exons_cov)
            pm.clean_add(introns_cov)

            # need Total Reads divided by 1M
            ar = float(pm.pipestat.retrieve(None,"Aligned_reads"))
            scaling_factor = float(ar/1000000)

            exons_rpkm = os.path.join(QC_folder, args.sample_name +
                                      "_exons_rpkm.bed")
            introns_rpkm = os.path.join(QC_folder, args.sample_name +
                                        "_introns_rpkm.bed")

            # determine exonic RPKM for individual genes
            if os.path.exists(exons_cov):
                cmd = ("awk -v OFS='\t' '{chrom[$4] = $1; " +
                       "if($4!=prev4) {chromStart[$4] = $2} " +
                       "strand[$4] = $6; " +
                       "readCount[$4] += $7; " + 
                       "exonCount[$4] += 1; " +
                       "geneSizeKB[$4] += " +
                       "(sqrt(($3-$2+0.00000001)^2)/1000); " +
                       "gene[$4] = $4; " +
                       "chromEnd[$4]=$3; " +
                       "prev4=$4} END " +
                       "{ for (a in readCount) " +
                       "{ print chrom[a], chromStart[a], " +
                       "chromEnd[a], gene[a], " +
                       "(readCount[a]/" + str(scaling_factor) +
                       ")/geneSizeKB[a], strand[a]}}' " +
                        exons_cov + " | awk '$5>0' | sort -k4 > " +
                        exons_rpkm)
                pm.run(cmd, exons_rpkm, nofail=True)
                pm.clean_add(exons_rpkm)

            # determine intronic RPKM for individual genes
            if os.path.exists(introns_cov):
                cmd = ("awk -v OFS='\t' '{chrom[$4] = $1; " +
                       "if($4!=prev4) {chromStart[$4] = $2} " +
                       "strand[$4] = $6; " +
                       "readCount[$4] += $7; " + 
                       "exonCount[$4] += 1; " +
                       "geneSizeKB[$4] += " +
                       "(sqrt(($3-$2+0.00000001)^2)/1000); " +
                       "gene[$4] = $4; " +
                       "chromEnd[$4]=$3; " +
                       "prev4=$4} END " +
                       "{ for (a in readCount) " +
                       "{ print chrom[a], chromStart[a], " +
                       "chromEnd[a], gene[a], " +
                       "(readCount[a]/" + str(scaling_factor) +
                       ")/geneSizeKB[a], strand[a]}}' " +
                        introns_cov + " | awk '$5>0' | sort -k4 > " +
                        introns_rpkm)
                pm.run(cmd, introns_rpkm, nofail=True)
                pm.clean_add(introns_rpkm)

            # join intron, exon RPKM on gene name and calculate ratio
            if os.path.exists(exons_rpkm) and os.path.exists(introns_rpkm):
                cmd = ("join --nocheck-order -a1 -a2 -j4 " +
                       introns_rpkm + " " + exons_rpkm + " | " +
                       "awk -v OFS='\t' " +
                       "'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}'" +
                       " | sort -k1,1 -k2,2n > " + intron_exon)
                pm.run(cmd, intron_exon, nofail=True)

            # report median ratio
            if os.path.exists(intron_exon):
                cmd = ("awk '{print $5}' " + intron_exon +
                       " | sort -n | awk ' { a[i++]=$1; }" +
                       " END { x=int((i+1)/2);" +
                       " if (x < (i+1)/2) print (a[x-1]+a[x])/2;" +
                       " else print a[x-1]; }'")
                mrna_con = float(pm.checkprint(cmd))
                # pm.report_result("mRNA_contamination", round(mrna_con, 2))
                pm.pipestat.report(values={"mRNA_contamination": round(mrna_con, 2)})

        # plot mRNA contamination distribution
        mRNApdf = os.path.join(QC_folder,
            args.sample_name + "_mRNA_contamination.pdf")
        mRNApng = os.path.join(QC_folder,
            args.sample_name + "_mRNA_contamination.png")

        if _itsa_file(intron_exon_gz) and not _itsa_file(mRNApdf):
            cmd = (ngstk.ziptool + " -d " + intron_exon_gz)
            pm.run(cmd, intron_exon)

        mRNAplot = [tools.Rscript, tool_path("PEPPRO.R"), "mrna",
                    "-i", intron_exon, "--annotate"]
        cmd = build_command(mRNAplot)
        pm.run(cmd, mRNApdf, nofail=False)
        # pm.report_object("mRNA contamination", mRNApdf, anchor_image=mRNApng)
        pm.pipestat.report(values={"mRNA_contamination_plot": {"path": mRNApdf, "thumbnail_path": mRNApng, "title": "mRNA contamination"}})

        if not is_gzipped(intron_exon):
            cmd = (ngstk.ziptool + " -f " + intron_exon)
            intron_exon = intron_exon + ".gz"
            pm.run(cmd, intron_exon)

    ############################################################################
    #                             Produce BigWigs                              #
    ############################################################################
    genome_fq = rgc.seek(args.genome_assembly,
                              asset_name="fasta",
                              seek_key="fasta")
    plus_exact_bw = os.path.join(
        signal_folder, args.sample_name + "_plus_exact_body_0-mer.bw")
    plus_smooth_bw = os.path.join(
        signal_folder, args.sample_name + "_plus_smooth_body_0-mer.bw")
    minus_exact_bw = os.path.join(
        signal_folder, args.sample_name + "_minus_exact_body_0-mer.bw")
    minus_smooth_bw = os.path.join(
        signal_folder, args.sample_name + "_minus_smooth_body_0-mer.bw")
    
    if not args.sob:
        # If not scaling we don't need to use seqOutBias to generate the
        # separate strand bigWigs; just convert the BAM's directly with 
        # bamSitesToWig.py which uses UCSC wigToBigWig
        pm.timestamp("### Produce bigWig files")
        
        # need Total Reads divided by 1M
        ar = float(pm.pipestat.retrieve(None,"Aligned_reads"))
        scaling_factor = float(ar/1000000)

        wig_cmd_callable = ngstk.check_command("wigToBigWig")

        if wig_cmd_callable:
            cmd1 = tools.samtools + " index " + plus_bam
            cmd2 = tool_path("bamSitesToWig.py")
            cmd2 += " -i " + plus_bam
            cmd2 += " -c " + res.chrom_sizes
            cmd2 += " -o " + plus_exact_bw  # DEBUG formerly smoothed " -w " + plus_bw
            cmd2 += " -w " + plus_smooth_bw  
            cmd2 += " -p " + str(int(max(1, int(pm.cores) * 2/3)))
            cmd2 += " --variable-step"
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                cmd2 += " --tail-edge"
            if args.scale:
                cmd2 += " --scale " + str(ar)
            pm.run([cmd1, cmd2], [plus_exact_bw, plus_smooth_bw])

            cmd3 = tools.samtools + " index " + minus_bam
            cmd4 = tool_path("bamSitesToWig.py")
            cmd4 += " -i " + minus_bam
            cmd4 += " -c " + res.chrom_sizes
            cmd4 += " -o " + minus_exact_bw # DEBUG formerly smoothed " -w " + minus_bw
            cmd4 += " -w " + minus_smooth_bw  
            cmd4 += " -p " + str(int(max(1, int(pm.cores) * 2/3)))
            cmd4 += " --variable-step"
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                cmd4 += " --tail-edge"
            if args.scale:
                cmd4 += " --scale " + str(ar)
            pm.run([cmd3, cmd4], [minus_exact_bw, minus_smooth_bw])
        else:
            print("Skipping signal track production -- Could not call \'wigToBigWig\'.")
            print("Check that you have the required UCSC tools in your PATH.")
    else:
        pm.debug("The max read length is {}".format(str(max_len)))

        # seqOutBias needs a working directory, we'll make that temporary
        tempdir = tempfile.mkdtemp(dir=signal_folder)
        os.chmod(tempdir, 0o771)
        pm.clean_add(tempdir)

        pm.timestamp("### Use seqOutBias to produce bigWig files")

        seqtable = os.path.join(signal_folder, (args.genome_assembly + ".tbl"))
        seqtable_flag = os.path.join(signal_folder, "seqtable_complete.flag")

        if not os.path.exists(seqtable_flag) and os.path.exists(seqtable):
            os.remove(seqtable)
            
        seqtable_cmd = build_command([
            (tools.seqoutbias, "seqtable"),
            res.fasta,
            str("--tallymer=" + res.search_file),
            str("--gt-workdir=" + tempdir),
            str("--read-size=" + str(max_len)),
            str("--out=" + seqtable)
        ])
        
        complete_cmd = "touch " + seqtable_flag

        pm.run([seqtable_cmd, complete_cmd], seqtable_flag)
        pm.clean_add(seqtable)

        if args.scale:
            scale_plus_chunks = [
                (tools.seqoutbias, "scale"),
                seqtable,
                plus_bam,
                "--skip-bed",
                str("--bw=" + plus_exact_bw)
            ]
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                scale_plus_chunks.extend([("--tail-edge")])
            scale_plus_cmd = build_command(scale_plus_chunks)

            scale_minus_chunks = [
                (tools.seqoutbias, "scale"),
                seqtable,
                minus_bam,
                "--skip-bed",
                str("--bw=" + minus_exact_bw),
            ]
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                scale_minus_chunks.extend([("--tail-edge")])
            scale_minus_cmd = build_command(scale_minus_chunks)
        else:
            scale_plus_chunks = [
                (tools.seqoutbias, "scale"),
                seqtable,
                plus_bam,
                "--no-scale",
                "--skip-bed",
                str("--bw=" + plus_exact_bw)
            ]
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                scale_plus_chunks.extend([("--tail-edge")])
            scale_plus_cmd = build_command(scale_plus_chunks)

            scale_minus_chunks = [
                (tools.seqoutbias, "scale"),
                seqtable,
                minus_bam,
                "--no-scale",
                "--skip-bed",
                str("--bw=" + minus_exact_bw),
            ]
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                scale_minus_chunks.extend([("--tail-edge")])
            scale_minus_cmd = build_command(scale_minus_chunks)

        pm.run([scale_plus_cmd, scale_minus_cmd], minus_exact_bw)

    # Remove potentially empty folders
    if os.path.exists(raw_folder) and os.path.isdir(raw_folder):
        if not os.listdir(raw_folder):
            pm.clean_add(raw_folder)
    if os.path.exists(fastqc_folder) and os.path.isdir(fastqc_folder):
        if not os.listdir(fastqc_folder):
            pm.clean_add(fastqc_folder)

    ############################################################################
    #                            PIPELINE COMPLETE!                            #
    ############################################################################
    pm.stop_pipeline()


if __name__ == '__main__':
    pm = None
    ngstk = None
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
