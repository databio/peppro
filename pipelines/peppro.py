#!/usr/bin/env python
"""
PEPPRO - Run-on sequencing pipeline
"""

__author__ = ["Jason Smith", "Nathan Sheffield", "Mike Guertin"]
__email__ = "jasonsmith@virginia.edu"
__version__ = "0.8.1"


from argparse import ArgumentParser
import os
import sys
import re
import tempfile
import tarfile
import pypiper
from pypiper import build_command
from refgenconf import RefGenConf as RGC, select_genome_config

TOOLS_FOLDER = "tools"
RUNON_SOURCE_PRO = ["PRO", "pro", "PRO-SEQ", "PRO-seq", "proseq", "PROSEQ"]
RUNON_SOURCE_GRO = ["GRO", "gro", "groseq", "GROSEQ", "GRO-SEQ", "GRO-seq"]
RUNON_SOURCE = RUNON_SOURCE_PRO + RUNON_SOURCE_GRO

ADAPTER_REMOVERS = ["fastp", "cutadapt"]
DEDUPLICATORS = ["seqkit", "fqdedup"]
TRIMMERS = ["seqtk", "fastx"]

DEFAULT_REMOVER = "fastp"
DEFAULT_DEDUPLICATOR = "seqkit"
DEFAULT_TRIMMER = "seqtk"

BT2_IDX_KEY = "bowtie2_index"
DEFAULT_UMI_LEN = 0
DEFAULT_MAX_LEN = 30

def parse_arguments():
    """
    Parse command-line arguments passed to the pipeline.
    """
    # Argument Parsing from yaml file
    ###########################################################################
    parser = ArgumentParser(description='PEPPRO version ' + __version__)
    parser = pypiper.add_pypiper_args(parser, groups=
        ['pypiper', 'looper', 'ngs'],
        required=["input", "genome", "sample-name", "output-parent"])

    # Pipeline-specific arguments
    parser.add_argument("--protocol", dest="protocol",
                        default="pro", choices=RUNON_SOURCE,
                        help="Run on sequencing type.")

    parser.add_argument("--adapter", dest="adapter",
                        default=DEFAULT_REMOVER, choices=ADAPTER_REMOVERS,
                        help="Name of adapter removal program."
                        "Default: {}".format(DEFAULT_REMOVER))

    parser.add_argument("--dedup", dest="dedup",
                        default=DEFAULT_DEDUPLICATOR, choices=DEDUPLICATORS,
                        help="Program to use to duplicate reads. "
                            "Default: {}".format(DEFAULT_DEDUPLICATOR))

    parser.add_argument("--trimmer", dest="trimmer",
                        default=DEFAULT_TRIMMER, choices=TRIMMERS,
                        help="Name of read trimming program. "
                            "Default: {}".format(DEFAULT_TRIMMER))

    parser.add_argument("--umi-fastp", action='store_true', default=False,
                        help="Use fastp to remove UMIs. Default: Use "
                        "tool selected for 'trimmer' option.")
    
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
                        help="Scale output with seqOutBias when producing"
                             " signal tracks.")

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

    parser.add_argument("--complexity", action='store_false', default=True,
                        dest="complexity",
                        help="Disable library complexity calculation (faster).")

    parser.add_argument("-V", "--version", action="version",
                        version="%(prog)s {v}".format(v=__version__))

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        raise SystemExit

    return args


def _remove_adapters(args, tools, read2, fq_file, outfolder):
    """
    A helper function to build a command for adapter removal.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param bool read2: if True, do not deduplicate and do not retain
        intermediate files
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to remove adapters
    """

    sname = args.sample_name  # for concise code

    cutadapt_folder = os.path.join(outfolder, "cutadapt")
    cutadapt_report = os.path.join(cutadapt_folder, sname + "_cutadapt.txt")

    fastp_folder = os.path.join(outfolder, "fastp")
    fastq_folder = os.path.join(outfolder, "fastq")
    noadap_fastq = os.path.join(fastq_folder, sname + "_R1_noadap.fastq")

    fastp_pfx = os.path.join(fastp_folder, sname + "_R1_fastp_adapter")
    fastp_pfx_R2 = os.path.join(fastp_folder, sname + "_R2_fastp_adapter")

    fastp_report_txt = fastp_pfx + ".txt"
    fastp_report_html = fastp_pfx + ".html"
    fastp_report_json = fastp_pfx + ".json"
    fastp_report_html_R2 = fastp_pfx_R2 + ".html"
    fastp_report_json_R2 = fastp_pfx_R2 + ".json"

    # Setup report output folders
    if args.adapter == "cutadapt":
        ngstk.make_dir(cutadapt_folder)
        adapter_report = cutadapt_report
    elif args.adapter == "fastp":
        ngstk.make_dir(fastp_folder)
        adapter_report = fastp_report_txt

    # Check quality encoding for use with FastX_Tools
    if args.trimmer == "fastx":
        encoding = _guess_encoding(fq_file)

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
                ("--adapter_sequence", "GATCGTCGGACTGTAGAACTCTGAAC"),
                ("--length_required", (18 + int(float(args.umi_len)))),
                ("--html", fastp_report_html_R2),
                ("--json", fastp_report_json_R2),
                ("--report_title", ("'" + sname + "'"))
            ])
        else:
            adapter_cmd_chunks.extend([
                ("--adapter_sequence", "TGGAATTCTCGGGTGCCAAGG"),
                ("--length_required", (18 + int(float(args.umi_len)))),
                ("--html", fastp_report_html),
                ("--json", fastp_report_json),
                ("--report_title", ("'" + sname + "'"))
            ])
        # If calculating library complexity and this is read 1 or single-end,
        # must produce an intermediate file.
        if args.complexity and args.umi_len > 0 and not read2 :
            adapter_cmd_chunks.extend([("-o", noadap_fastq)])
        else:
            adapter_cmd_chunks.extend([("--stdout")])

        adapter_cmd_chunks.extend([
            (") 2>", fastp_report_txt)
        ])

        adapter_cmd = build_command(adapter_cmd_chunks)

    elif args.adapter == "cutadapt":
        cut_version = float(pm.checkprint("cutadapt --version"))
        if read2:
            adapter_cmd_chunks = [tools.cutadapt]
            # old versions of cutadapt can not use multiple cores
            if cut_version >= 1.15:
                adapter_cmd_chunks.extend([("-j", str(pm.cores))])
            adapter_cmd_chunks.extend([
                    ("-m", (18 + int(float(args.umi_len)))),
                    ("-a", "GATCGTCGGACTGTAGAACTCTGAAC"),
                    fq_file
            ])
        else:
            if args.complexity and args.umi_len > 0:
                adapter_cmd_chunks = ["(" + tools.cutadapt]
                # old versions of cutadapt can not use multiple cores
                if cut_version >= 1.15:
                    adapter_cmd_chunks.extend([("-j", str(pm.cores))])
                adapter_cmd_chunks.extend([
                    ("-m", (18 + int(float(args.umi_len)))),
                    ("-a", "TGGAATTCTCGGGTGCCAAGG"),
                    fq_file,
                    ("-o", noadap_fastq + ")"),
                    (">", cutadapt_report)
                ])
            else:
                adapter_cmd_chunks = [tools.cutadapt]
                # old versions of cutadapt can not use multiple cores
                if cut_version >= 1.15:
                    adapter_cmd_chunks.extend([("-j", str(pm.cores))])
                adapter_cmd_chunks.extend([
                    ("-m", (18 + int(float(args.umi_len)))),
                    ("-a", "TGGAATTCTCGGGTGCCAAGG"),
                    fq_file
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
                ("--adapter_sequence", "GATCGTCGGACTGTAGAACTCTGAAC"),
                ("--length_required", (18 + int(float(args.umi_len)))),
                ("--html", fastp_report_html_R2),
                ("--json", fastp_report_json_R2),
                ("--report_title", ("'" + sname + "'"))
            ])
        else:
            adapter_cmd_chunks.extend([
                ("--adapter_sequence", "TGGAATTCTCGGGTGCCAAGG"),
                ("--length_required", (18 + int(float(args.umi_len)))),
                ("--html", fastp_report_html),
                ("--json", fastp_report_json),
                ("--report_title", ("'" + sname + "'"))
            ])
        # If calculating library complexity and this is read 1 or single-end,
        # must produce a physical output file.
        if args.complexity and not read2:
            adapter_cmd_chunks.extend([("-o", noadap_fastq)])
        else:
            adapter_cmd_chunks.extend([("--stdout")])

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
    if args.complexity and args.umi_len > 0:
        if args.dedup == "seqkit":
            dedup_cmd_chunks = [
                (tools.seqkit, "rmdup"),
                ("--threads", str(pm.cores)),
                "--by-seq",
                "--ignore-case",
                "-o"           
            ]
            if args.complexity:
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

    # Only call this when args.complexity and args.umi_len > 0
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    dedup_fastq = os.path.join(fastq_folder, sname + "_R1_dedup.fastq")
    #processed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")
    #trimmed_fastq = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")
    processed_fastq = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")

    fastp_folder = os.path.join(outfolder, "fastp")    
    umi_report = os.path.join(fastp_folder, sname + "_R1_rmUmi.html")
    umi_json = os.path.join(fastp_folder, sname + "_R1_rmUmi.json")

    # Check quality encoding for use with FastX_Tools
    if args.trimmer == "fastx":
        encoding = _guess_encoding(fq_file)

    if args.umi_fastp and args.adapter == "fastp":
        trim_cmd_chunks = [
            tools.fastp,
            ("--thread", str(pm.cores)),
            ("-i", dedup_fastq),
            "--stdout",
            "--umi",
            ("--umi_loc", "read1"),
            ("--umi_len", args.umi_len),
            ("--html", umi_report),
            ("--json", umi_json),
            "|",
            (tools.seqtk, "trimfq")
        ]

        if args.max_len != -1:
            trim_cmd_chunks.extend([("-L", args.max_len)])

        trim_cmd_chunks.extend(["-"])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([(">", processed_fastq)])
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                ("-r", "-"),
                (">", processed_fastq)
            ])
    else:
        if args.umi_fastp and args.adapter != "fastp":
            print("To remove UMI intelligently, you must process adapters using 'fastp'")
            print("Defaulting to removing the first {} "
                  "bp instead via trimming".format(str(args.umi_len)))
        if args.trimmer == "seqtk":
            trim_cmd_chunks = [
                tools.seqtk,
                "trimfq",
                ("-b", str(args.umi_len))
            ]

            if args.max_len != -1:
                trim_cmd_chunks.extend([("-L", str(args.max_len))])

            trim_cmd_chunks.extend([dedup_fastq])

            # Do not reverse complement for GRO-seq
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([(">", processed_fastq)])                           
            else:
                trim_cmd_chunks.extend([
                    "|",
                    tools.seqtk,
                    ("seq", "-r"),
                    ("-", ">"),
                    processed_fastq
                ])
        elif args.trimmer == "fastx":
            trim_tool = tools.fastx + "_trimmer"
            rc_tool = tools.fastx + "_reverse_complement"
            trim_cmd_chunks = [trim_tool]

            if encoding == "Illumina-1.8":
                trim_cmd_chunks.extend([
                    ("-Q", str(33))
                ])

            trim_cmd_chunks.extend([("-f", str(int(float(args.umi_len)) + 1))])

            if args.max_len != -1:
                trim_cmd_chunks.extend([
                    ("-l", (str(int(float(args.max_len)) + int(float(args.umi_len)))))
                ])

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

            if args.max_len != -1:
                trim_cmd_chunks.extend([("-L", str(args.max_len))])

            trim_cmd_chunks.extend([dedup_fastq])

            # Do not reverse complement for GRO-seq
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([(">", processed_fastq)])                           
            else:
                trim_cmd_chunks.extend([
                    "|",
                    tools.seqtk,
                    ("seq", "-r"),
                    ("-", ">"),
                    processed_fastq
                ])

    trim_cmd = build_command(trim_cmd_chunks)
    pm.debug("trim_deduplicated_cmd: {}".format(build_command(trim_cmd_chunks)))

    return trim_cmd


def _trim_adapter_files(args, tools, fq_file, outfolder):
    """
    A helper function to build a command for read trimming using fastq files
    without deduplication.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return str: command to trim adapter trimmed files
    """
    # Need undeduplicated results for complexity calculation
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    noadap_fastq = os.path.join(fastq_folder, sname + "_R1_noadap.fastq")
    trimmed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")
    #trimmed_fastq = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")
    #processed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")

    fastp_folder = os.path.join(outfolder, "fastp") 
    umi_report = os.path.join(fastp_folder, sname + "_R1_rmUmi.html")
    umi_json = os.path.join(fastp_folder, sname + "_R1_rmUmi.json")

    # Check quality encoding for use with FastX_Tools
    if args.trimmer == "fastx":
        encoding = _guess_encoding(fq_file)
    
    if args.umi_fastp and args.adapter == "fastp":
        # Use FastP
        # Still requires seqtk for reverse complementation
        trim_cmd_chunks = [
            tools.fastp,
            ("--thread", str(pm.cores)),
            ("-i", noadap_fastq),
            "--stdout",
            "--umi",
            ("--umi_loc", "read1"),
            ("--umi_len", args.umi_len),
            ("--html", umi_report),
            ("--json", umi_json),
            "|",
            (tools.seqtk, "trimfq")
        ]

        if args.max_len != -1:
            trim_cmd_chunks.extend([("-L", args.max_len)])

        trim_cmd_chunks.extend(["-"])

        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([
                (">", trimmed_fastq)
            ])
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                ("-r", "-"),
                (">", trimmed_fastq)
            ])
    else:
        if args.umi_fastp and args.adapter != "fastp":
            print("To remove UMI intelligently, you must process adapters using 'fastp'")
            print("Defaulting to removing the first {} "
                  "bp instead via trimming".format(str(args.umi_len)))
        if args.trimmer == "seqtk":
            trim_cmd_chunks = [
                tools.seqtk,
                "trimfq",
                ("-b", str(args.umi_len))
            ]
            if args.max_len != -1:
                trim_cmd_chunks.extend([
                    ("-L", str(args.max_len))
                ])
            
            trim_cmd_chunks.extend([noadap_fastq])
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([
                    (">", trimmed_fastq)
                ])
            else:
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "seq"),
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
            trim_cmd_chunks.extend([
                ("-f", str(int(float(args.umi_len)) + 1))
            ])
            if args.max_len != -1:
                trim_cmd_chunks.extend([
                    ("-l", (str(int(float(args.max_len)) + int(float(args.umi_len)))))
                ])
            # Need undeduplicated results for complexity calculation
            
            trim_cmd_chunks.extend([
                ("-i", noadap_fastq)
            ])
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([
                    ("-o", trimmed_fastq)
                ])
            else:
                if encoding == "Illumina-1.8":
                    trim_cmd_chunks.extend([
                        ("-Q", str(33))
                    ])
                trim_cmd_chunks.extend([
                    ("-o", trimmed_fastq)
                ])
        # Default to seqtk
        else:
            trim_cmd_chunks = [
                tools.seqtk,
                "trimfq",
                ("-b", str(args.umi_len))
            ]
            if args.max_len != -1:
                trim_cmd_chunks.extend([
                    ("-L", str(args.max_len))
                ])
            
            trim_cmd_chunks.extend([noadap_fastq])
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([
                    (">", trimmed_fastq)
                ])
            else:
                trim_cmd_chunks.extend([
                    "|",
                    (tools.seqtk, "seq"),
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

    # Only call this when NOT args.complexity or NOT args.umi_len > 0
    sname = args.sample_name  # for concise code

    fastq_folder = os.path.join(outfolder, "fastq")
    if read2:
        processed_fastq = os.path.join(fastq_folder, sname + "_R2_trimmed.fastq")
    else:
        processed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")

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

    if args.umi_fastp and args.adapter == "fastp":
        # There are no intermediate files, just pipes
        trim_cmd_chunks = [
            tools.fastp,
            ("--thread", str(pm.cores)),
            ("--stdin", "--stdout"),
            "--umi"
        ]

        if read2:
            trim_cmd_chunks.extend([("--umi_loc", "paired_end")])
        else:
            trim_cmd_chunks.extend([("--umi_loc", "read1")])

        trim_cmd_chunks.extend([
            ("--umi_len", args.umi_len),
            ("--html", umi_report),
            ("--json", umi_json),
            "|",
            (tools.seqtk, "trimfq")
        ])

        if args.max_len != -1:
            trim_cmd_chunks.extend([("-L", args.max_len)])

        trim_cmd_chunks.extend(["-"])

        # Do not reverse complement for GRO-seq
        if args.protocol.lower() in RUNON_SOURCE_GRO:
            trim_cmd_chunks.extend([(">", processed_fastq)])
        else:
            trim_cmd_chunks.extend([
                "|",
                (tools.seqtk, "seq"),
                ("-r", "-"),
                (">", processed_fastq)
            ])
    # If args.complexity and args.umi_len > 0 retain intermediate files
    else:
        if args.umi_fastp and args.adapter != "fastp":
            print("To remove UMI intelligently, you must process adapters using 'fastp'")
            print("Defaulting to removing the first {} "
                  "bp instead via trimming".format(str(args.umi_len)))
        if args.trimmer == "seqtk":
            trim_cmd_chunks = [
                tools.seqtk,
                "trimfq"
            ]

            if read2:
                trim_cmd_chunks.extend([("-e", str(args.umi_len))])
            else:
                trim_cmd_chunks.extend([("-b", str(args.umi_len))])
                if args.max_len != -1:
                    trim_cmd_chunks.extend([("-L", str(args.max_len))])

            trim_cmd_chunks.extend(["-"])

            # Do not reverse complement for GRO-seq
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([(">", processed_fastq)])                           
            else:
                trim_cmd_chunks.extend([
                    "|",
                    tools.seqtk,
                    ("seq", "-r"),
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
                if args.max_len != -1:
                    trim_cmd_chunks.extend([
                        ("-l", (str(int(float(args.max_len)) + int(float(args.umi_len)))))
                    ])

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
                if args.max_len != -1:
                    trim_cmd_chunks.extend([("-L", str(args.max_len))])

            trim_cmd_chunks.extend(["-"])

            # Do not reverse complement for GRO-seq
            if args.protocol.lower() in RUNON_SOURCE_GRO:
                trim_cmd_chunks.extend([(">", processed_fastq)])                           
            else:
                trim_cmd_chunks.extend([
                    "|",
                    tools.seqtk,
                    ("seq", "-r"),
                    ("-", ">"),
                    processed_fastq
                ])

    trim_cmd = build_command(trim_cmd_chunks)
    pm.debug("trim_pipes_cmd: {}".format(build_command(trim_cmd_chunks)))
    pm.debug("trim_pipes_cmd read2 status: {}".format(read2))

    return trim_cmd


def _process_fastq(args, tools, read2, fq_file, outfolder):
    """
    A helper function to prepare read files for downstream processing.

    :param argparse.Namespace args: binding between option name and argument,
        e.g. from parsing command-line options
    :param looper.models.AttributeDict tools: binding between tool name and
        value, e.g. for tools/resources used by the pipeline
    :param bool read2: if True, do not deduplicate and do not retain
        intermediate files
    :param str fq_file: path to FASTQ file
    :param str outfolder: path to output directory for the pipeline
    :return (str, str): pair (R1, R2) of paths to FASTQ files
    """
    # Create names for processed FASTQ files.
    fastq_folder = os.path.join(outfolder, "fastq")
    fastqc_folder = os.path.join(outfolder, "fastqc")
    fastp_folder = os.path.join(outfolder, "fastp")

    sname = args.sample_name  # for concise code

    noadap_fastq = os.path.join(fastq_folder, sname + "_R1_noadap.fastq")
    dedup_fastq = os.path.join(fastq_folder, sname + "_R1_dedup.fastq")
    trimmed_fastq = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")
    trimmed_fastq_R2 = os.path.join(fastq_folder, sname + "_R2_trimmed.fastq")
    trimmed_dups_fastq_R2 = os.path.join(fastq_folder, sname + "_R2_trimmed_dups.fastq")
    processed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")

    if args.adapter == "cutadapt":
        cutadapt_folder = os.path.join(outfolder, "cutadapt")
        cutadapt_report = os.path.join(cutadapt_folder, sname + "_cutadapt.txt")
        adapter_report = cutadapt_report
    else:
        adapter_report = os.path.join(fastqc_folder, sname + "_R1_rmAdapter.txt")

    fastp_pfx = os.path.join(fastp_folder, sname + "_R1_fastp_adapter")
    fastp_report_txt = fastp_pfx + ".txt"
    fastp_report_html = fastp_pfx + ".html"

    adapter_command = _remove_adapters(args, tools, read2, fq_file, outfolder)
    pm.debug("Adapter command: {}".format(adapter_command))
    pm.debug("Read2 status: {}".format(read2))

    if read2:
        trim_command = _trim_pipes(args, tools, True, fq_file, outfolder)
    elif args.complexity and args.umi_len > 0:
        deduplicate_command = _deduplicate(args, tools, fq_file, outfolder)
        pm.debug("Dedup command: {}".format(deduplicate_command))
        trim_command = _trim_adapter_files(args, tools, fq_file, outfolder)
        trim_command2 = _trim_deduplicated_files(args, tools, fq_file, outfolder)
    else:
        trim_command = _trim_pipes(args, tools, False, fq_file, outfolder)

    def report_fastq():
        """
        Report QC metrics on intermediate steps of fastq file preparation
        """
        if args.adapter == "cutadapt":
            adapter_term = "Reads with adapters:"
            too_short_term = "Reads that were too short:"
            total_bases_term = "Total basepairs processed:"

            ac_cmd = ("grep '" + adapter_term + "' " +
                      cutadapt_report + " | awk '{print $(NF-1)}'")
            ts_cmd = ("grep '" + too_short_term + "' " +
                      cutadapt_report + " | awk '{print $(NF-1)}'")
            bases = ("grep '" + total_bases_term + "' " +
                     cutadapt_report + " | awk '{print $(NF-1)}'")
            adapter_bases = ("awk '{sum+=$1*$2} END {printf \"%.0f\", sum}' " +
                             cutadapt_report)

        else:  # default to fastp
            adapter_term = "reads with adapter trimmed:"
            too_short_term = "reads failed due to too short:"
            total_bases_term = "total bases:"

            ac_cmd = ("grep '" + adapter_term + "' " +
                       fastp_report_txt + " | head -n 1 | awk '{print $NF}'")
            ts_cmd = ("grep '" + too_short_term + "' " +
                       fastp_report_txt + " | head -n 1 | awk '{print $NF}'")
            bases = ("grep '" + total_bases_term + "' " +
                     fastp_report_txt + " | head -n 1 | awk '{print $NF}'")
            adapter_bases = ("grep 'bases trimmed due to adapters:' " +
                             fastp_report_txt + " | awk '{print $NF}'")

            pm.report_object("FastP_report", fastp_report_html)

        ac = float(pm.checkprint(ac_cmd).replace(',',''))
        pm.report_result("Reads_with_adapter", ac)
        total_bases = float(pm.checkprint(bases).replace(',',''))
        total_adapter = float(pm.checkprint(adapter_bases).replace(',',''))
        pm.report_result("Pct_adapter_contamination",
                         round(float(total_adapter/total_bases), 2))

        ts = float(pm.checkprint(ts_cmd).replace(',',''))
        pm.report_result("Reads_too_short", ts)

        tr = int(ngstk.count_lines(noadap_fastq).strip())
        dr = int(ngstk.count_lines(dedup_fastq).strip())
        dups = max(0, (float(tr)/4 - float(dr)/4))
        pm.report_result("Duplicate_reads", dups)

    # Put it all together
    if read2:
        # cutadapt directs its report to stderr if the command lacks
        # a -o and the actual reads are directed to stdout.
        process_fastq_cmd2 = build_command([
            "(", adapter_command, "|", trim_command, ") 2> ", adapter_report])
        pm.debug("process_fastq_cmd2: {}".format(process_fastq_cmd2))
        pm.run(process_fastq_cmd2, trimmed_fastq_R2)
        cp_cmd = ("cp " + trimmed_fastq_R2 + " " + trimmed_dups_fastq_R2)
        pm.run(cp_cmd, trimmed_dups_fastq_R2)
        return trimmed_fastq_R2, trimmed_dups_fastq_R2
    else:
        if args.complexity and args.umi_len > 0:
            # This trim command DOES need the adapter file...
            pm.debug("\ntrim_command1: {} + {}\n".format(adapter_command, trim_command))
            pm.run([adapter_command, trim_command], processed_fastq,
                   follow=ngstk.check_trim(processed_fastq, False, None))
            # This needs to produce the trimmed_fastq file
            pm.debug("\ntrim_command2: {} + {} + {}\n".format(adapter_command, deduplicate_command, trim_command2))
            pm.run([adapter_command, deduplicate_command, trim_command2],
                   trimmed_fastq, follow=report_fastq)
            pm.clean_add(noadap_fastq)
            pm.clean_add(dedup_fastq)
            pm.clean_add(trimmed_fastq)
            return processed_fastq, trimmed_fastq
        else:
            pm.debug("\nELSE: trim_command: {} + {}\n".format(adapter_command, trim_command))
            process_fastq_cmd = build_command([
                "(", adapter_command, "|", trim_command, ") 2> ", adapter_report])
            pm.run(process_fastq_cmd, processed_fastq,
               follow=ngstk.check_trim(processed_fastq, False, None))
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
            pm.clean_add(cmd, out_fastq_tmp)

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
        else:
            cmd += " > /dev/null"
        cmd += ") 2>" + summary_file

        if paired:
            if args.keep or not useFIFO:
                pm.run([cmd, filter_pair], mapped_bam)
            else:
                pm.wait = False
                pm.run(filter_pair, [summary_file, out_fastq_r2_gz],
                       container=pm.container)
                pm.wait = True
                pm.run(cmd, [summary_file, out_fastq_r2_gz],
                       container=pm.container)
        else:
            if args.keep:
                pm.run(cmd, mapped_bam)
            else:
                # TODO: switch to this once filter_paired_fq works with SE
                #pm.run(cmd2, summary_file)
                #pm.run(cmd1, out_fastq_r1)
                pm.run(cmd, out_fastq_tmp_gz)

        pm.clean_add(out_fastq_tmp)

        if not dups:
            cmd = ("grep 'aligned exactly 1 time' " + summary_file +
                   " | awk '{print $1}'")
            align_exact = pm.checkprint(cmd)
            if align_exact:
                ar = float(align_exact)*2
            else:
                ar = 0

            # report aligned reads
            pm.report_result("Aligned_reads_" + assembly_identifier, ar)
            try:
                # wrapped in try block in case Trimmed_reads is not reported 
                # in this pipeline.
                tr = float(pm.get_stat("Trimmed_reads"))
            except:
                print("Trimmed reads is not reported.")
            else:
                res_key = "Alignment_rate_" + assembly_identifier
                pm.report_result(res_key, round(float(ar) * 100 / float(tr), 2))
        
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
                res[asset] = rgc.get_asset(reference, asset)
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
                    res[seek_key] = rgc.get_asset(args.genome_assembly,
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


###############################################################################
def main():
    """
    Main pipeline process.
    """

    args = parse_arguments()

    args.paired_end = args.single_or_paired == "paired"

    # Initialize, creating global PipelineManager and NGSTk instance for
    # access in ancillary functions outside of main().
    outfolder = os.path.abspath(
        os.path.join(args.output_parent, args.sample_name))
    global pm
    pm = pypiper.PipelineManager(
        name="PEPPRO", outfolder=outfolder, args=args, version=__version__)
    global ngstk
    ngstk = pypiper.NGSTk(pm=pm)

    # Convenience alias
    tools = pm.config.tools
    param = pm.config.parameters
    res   = pm.config.resources
    #sstructure = pm.sample_structure  # maybe possible in the future?

    # Check that the required tools are callable by the pipeline
    tool_list = [v for k,v in tools.items()]    # extract tool list
    tool_list = [t.replace('fastx', 'fastx_trimmer') for t in tool_list]
    tool_list = [t.replace('seqoutbias', 'seqOutBias') for t in tool_list]
    opt_tools = ["fqdedup", "fastx_trimmer", "seqOutBias"]
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
         "required":True},
        {"asset_name":"refgene_anno", "seek_key":"refgene_tss",
         "tag_name":"default", "arg":"TSS_name", "user_arg":"TSS-name",
         "required":False},
        {"asset_name":"ensembl_gtf", "seek_key":"ensembl_tss",
         "tag_name":"default", "arg":"ensembl_tss", "user_arg":"pi-tss",
         "required":False},
        {"asset_name":"ensembl_gtf", "seek_key":"ensembl_gene_body",
         "tag_name":"default", "arg":"ensembl_gene_body", "user_arg":"pi-body",
         "required":False},
        {"asset_name":"refgene_anno", "seek_key":"refgene_pre_mRNA",
         "tag_name":"default", "arg":"pre_name", "user_arg":"pre-name",
         "required":False},
        {"asset_name":"feat_annotation", "seek_key":"feat_annotation",
         "tag_name":"default", "arg":"anno_name", "user_arg":"anno-name",
         "required":False},
        {"asset_name":"refgene_anno", "seek_key":"refgene_exon",
         "tag_name":"default", "arg":"exon_name", "user_arg":"exon-name",
         "required":False},
        {"asset_name":"refgene_anno", "seek_key":"refgene_intron",
         "tag_name":"default", "arg":"intron_name", "user_arg":"intron-name",
         "required":False}
    ]
    res, rgc = _add_resources(args, res, check_list)

    # Adapter file can be set in the config; if left null, we use a default.
    # TODO: use this option or just specify directly the adapter sequence as I do now
    res.adapters = res.adapters or tool_path("PRO-seq_adapter.fa")

    param.outfolder = outfolder

    # Check that the input file(s) exist before continuing
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

    container = None

    ###########################################################################

    pm.report_result(
        "File_mb",
        round(ngstk.get_file_size(
            [x for x in [args.input, args.input2] if x is not None]), 2))
    pm.report_result("Read_type", args.single_or_paired)
    pm.report_result("Genome", args.genome_assembly)

    # PRO-seq pipeline
    # Each (major) step should have its own subfolder
    if args.protocol.lower() in RUNON_SOURCE_GRO:
        pm.info("Detected GRO input")
    elif args.protocol.lower() in RUNON_SOURCE_PRO:
        pm.info("Detected PRO input")
    else:
        pm.fail_pipeline(RuntimeError("Input protocol must be GRO or PRO."))

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
    if not pm.get_stat("Raw_reads") or args.new_start:
        pm.run(cmd, unaligned_fastq,
               follow=ngstk.check_fastq(
                   local_input_files, unaligned_fastq, args.paired_end),
               container=pm.container)
        pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)
    pm.info(local_input_files)
    untrimmed_fastq1 = out_fastq_pre + "_R1.fastq"
    untrimmed_fastq2 = out_fastq_pre + "_R2.fastq" if args.paired_end else None

    ############################################################################
    #                          Process read files                              #
    ############################################################################
    pm.timestamp("### FASTQ processing: ")
    cutadapt_folder = os.path.join(outfolder, "cutadapt")
    cutadapt_report = os.path.join(cutadapt_folder, args.sample_name + "_cutadapt.txt")
    repair_target = os.path.join(fastq_folder, "repaired.flag")
    dups_repair_target = os.path.join(fastq_folder, "dups_repaired.flag")

    # If single-end, must use cutadapt for plotting purposes
    if not args.paired_end:
        if args.adapter != "cutadapt":
            pm.warning("You set adapter arg to '{}' but you must select 'cutadapt'" 
                " for single end data. Overriding.".format(args.adapter))
        args.adapter = "cutadapt"

    # If we've already aligned to the primary genome, skip these steps unless
    # it's a --new-start
    if not pm.get_stat("Aligned_reads") or args.new_start:
        if args.paired_end:
            if args.complexity and args.umi_len > 0:
                unmap_fq1, unmap_fq1_dups = _process_fastq(args, tools, False,
                                                           untrimmed_fastq1,
                                                           outfolder=param.outfolder)
            else:
                unmap_fq1 = _process_fastq(args, tools, False,
                                           untrimmed_fastq1,
                                           outfolder=param.outfolder)
            unmap_fq2, unmap_fq2_dups = _process_fastq(args, tools, True,
                                                       untrimmed_fastq2,
                                                       outfolder=param.outfolder)

            # Gut check
            # Processing fastq should have trimmed the reads.
            tr = float(pm.get_stat("Trimmed_reads"))
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

            rr = float(pm.get_stat("Raw_reads"))
            cmd = (tools.fastqpair + " -t " + str(int(0.9*rr)) + " " + unmap_fq1 +
                   " " + unmap_fq2)
            pm.run(cmd, [r1_repair, r2_repair])
            pm.clean_add(r1_repair_single)
            pm.clean_add(r2_repair_single)
            cmd1 = ("mv " + r1_repair + " " + unmap_fq1)
            cmd2 = ("mv " + r2_repair + " " + unmap_fq2)
            cmd3 = ("touch repaired.flag")
            pm.run([cmd1, cmd2, cmd3], repair_target)

            # Re-pair the duplicates (but only if we could identify duplicates)
            if args.umi_len > 0:
                r1_dups_repair = os.path.join(
                    fastq_folder, args.sample_name + "_R1_trimmed.fastq.paired.fq")
                r2_dups_repair = os.path.join(
                    fastq_folder, args.sample_name + "_R2_trimmed_dups.fastq.paired.fq")

                r1_dups_repair_single = os.path.join(
                    fastq_folder, args.sample_name + "_R1_trimmed.fastq.single.fq")
                r2_dups_repair_single = os.path.join(
                    fastq_folder, args.sample_name + "_R2_trimmed_dups.fastq.single.fq")

                cmd = (tools.fastqpair + " -t " + str(int(0.9*rr)) + " " +
                       unmap_fq1_dups + " " + unmap_fq2_dups)
                pm.run(cmd, [r1_dups_repair, r2_dups_repair])
                pm.clean_add(r1_dups_repair_single)
                pm.clean_add(r2_dups_repair_single)
                cmd1 = ("mv " + r1_dups_repair + " " + unmap_fq1_dups)
                cmd2 = ("mv " + r2_dups_repair + " " + unmap_fq2_dups)
                cmd3 = ("touch dups_repaired.flag")
                pm.run([cmd1, cmd2, cmd3], dups_repair_target)
        else:
            if args.complexity and args.umi_len > 0:
                unmap_fq1, unmap_fq1_dups = _process_fastq(args, tools, False,
                                                           untrimmed_fastq1,
                                                           outfolder=param.outfolder)
                unmap_fq2 = ""
                unmap_fq2_dups = ""
            else:
                unmap_fq1 = _process_fastq(args, tools, False,
                                           untrimmed_fastq1,
                                           outfolder=param.outfolder)
                unmap_fq2 = ""

    pm.timestamp("### Plot adapter insertion distribution")
    if not args.adapter == "cutadapt":
        pm.info("Skipping sample degradation plotting...")
        pm.info("This requires using 'cutadapt' for adapter clipping.")
    elif not os.path.exists(cutadapt_report):
        pm.info("Skipping sample degradation plotting...")
        pm.info("Could not find {}.`".format(cutadapt_report))
    else:
        degradation_pdf = os.path.join(cutadapt_folder,
            args.sample_name + "_adapter_insertion_distribution.pdf")
        degradation_png = os.path.join(cutadapt_folder,
            args.sample_name + "_adapter_insertion_distribution.png")
        cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") + 
               " cutadapt -i " + cutadapt_report + " -o " + cutadapt_folder)
        pm.run(cmd, degradation_pdf, nofail=True)
        pm.report_object("Adapter insertion distribution", degradation_pdf,
                         anchor_image=degradation_png)

    pm.clean_add(os.path.join(fastq_folder, "*.fq"), conditional=True)
    pm.clean_add(os.path.join(fastq_folder, "*.fastq"), conditional=True)
    pm.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)

    ############################################################################
    #                  Map to any requested prealignments                      #
    ############################################################################
    # We recommend mapping to human_rDNA first for PRO-seq data
    pm.timestamp("### Prealignments")

    to_compress = []
    if not pm.get_stat("Aligned_reads") or args.new_start:
        if len(args.prealignments) == 0:
            print("You may use `--prealignments` to align to references before "
                  "the genome alignment step. See docs.")
        else:
            print("Prealignment assemblies: " + str(args.prealignments))
            # Loop through any prealignment references and map to them sequentially
            for reference in args.prealignments:
                if args.complexity and args.umi_len > 0:
                    if args.no_fifo:
                        unmap_fq1, unmap_fq2 = _align_with_bt2(
                            args, tools, args.paired_end, False, unmap_fq1,
                            unmap_fq2, reference,
                            assembly_bt2=os.path.join(
                                rgc.get_asset(reference, BT2_IDX_KEY), reference),
                            outfolder=param.outfolder,
                            aligndir="prealignments")

                        unmap_fq1_dups, unmap_fq2_dups = _align_with_bt2(
                            args, tools, args.paired_end, False, unmap_fq1_dups,
                            unmap_fq2_dups, reference,
                            assembly_bt2=os.path.join(
                                rgc.get_asset(reference, BT2_IDX_KEY), reference),
                            outfolder=param.outfolder,
                            aligndir="prealignments",
                            dups=True)
                        
                    else:
                        unmap_fq1, unmap_fq2 = _align_with_bt2(
                            args, tools, args.paired_end, True, unmap_fq1,
                            unmap_fq2, reference,
                            assembly_bt2=os.path.join(
                                rgc.get_asset(reference, BT2_IDX_KEY), reference),
                            outfolder=param.outfolder,
                            aligndir="prealignments")

                        unmap_fq1_dups, unmap_fq2_dups = _align_with_bt2(
                            args, tools, args.paired_end, True, unmap_fq1_dups,
                            unmap_fq2_dups, reference,
                            assembly_bt2=os.path.join(
                                rgc.get_asset(reference, BT2_IDX_KEY), reference),
                            outfolder=param.outfolder,
                            aligndir="prealignments",
                            dups=True)
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
                            assembly_bt2=os.path.join(
                                rgc.get_asset(reference, BT2_IDX_KEY), reference),
                            outfolder=param.outfolder,
                            aligndir="prealignments")
                    else:
                        unmap_fq1, unmap_fq2 = _align_with_bt2(
                            args, tools, args.paired_end, True,
                            unmap_fq1, unmap_fq2, reference,
                            assembly_bt2=os.path.join(
                                rgc.get_asset(reference, BT2_IDX_KEY), reference),
                            outfolder=param.outfolder,
                            aligndir="prealignments")
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

    mapping_genome_bam = os.path.join(
        map_genome_folder, args.sample_name + "_sort.bam")
    mapping_genome_bam_temp = os.path.join(
        map_genome_folder, args.sample_name + "_temp.bam")
    failQC_genome_bam = os.path.join(
        map_genome_folder, args.sample_name + "_fail_qc.bam")
    unmap_genome_bam = os.path.join(
        map_genome_folder, args.sample_name + "_unmap.bam")

    mapping_genome_bam_dups = os.path.join(
        map_genome_folder, args.sample_name + "_sort_dups.bam")
    mapping_genome_bam_temp_dups = os.path.join(
        map_genome_folder, args.sample_name + "_temp_dups.bam")
    failQC_genome_bam_dups = os.path.join(
        map_genome_folder, args.sample_name + "_fail_qc_dups.bam")
    unmap_genome_bam_dups = os.path.join(
        map_genome_folder, args.sample_name + "_unmap_dups.bam")

    temp_mapping_index = os.path.join(mapping_genome_bam_temp + ".bai")
    temp_mapping_index_dups = os.path.join(mapping_genome_bam_temp_dups + ".bai")

    mito_name = ["chrM", "chrMT", "M", "MT", "rCRSd", "rCRSd_3k"]

    bt2_options = " --very-sensitive"
    bt2_options += " -X 2000"

    # samtools sort needs a temporary directory
    tempdir = tempfile.mkdtemp(dir=map_genome_folder)
    os.chmod(tempdir, 0o771)
    pm.clean_add(tempdir)

    # Skip if this is a recovery and the following has already occurred
    if not pm.get_stat("Aligned_reads") or args.new_start:
        # check input for zipped or not
        if pypiper.is_gzipped_fastq(unmap_fq1):
            cmd = (ngstk.ziptool + " -d " + (unmap_fq1 + ".gz"))
            pm.run(cmd, mapping_genome_bam)
        if args.paired_end:
            if pypiper.is_gzipped_fastq(unmap_fq2):
                cmd = (ngstk.ziptool + " -d " + (unmap_fq2 + ".gz"))
                pm.run(cmd, mapping_genome_bam)

        cmd = tools.bowtie2 + " -p " + str(pm.cores)
        cmd += bt2_options
        cmd += " --rg-id " + args.sample_name
        cmd += " -x " + os.path.join(
            rgc.get_asset(args.genome_assembly, BT2_IDX_KEY),
            args.genome_assembly)
        if args.paired_end:
            cmd += " --rf -1 " + unmap_fq1 + " -2 " + unmap_fq2
        else:
            cmd += " -U " + unmap_fq1
        cmd += " | " + tools.samtools + " view -bS - -@ 1 "
        cmd += " | " + tools.samtools + " sort - -@ 1"
        cmd += " -T " + tempdir
        cmd += " -o " + mapping_genome_bam_temp

        if args.complexity and args.umi_len > 0:
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
            cmd_dups += " -x " + os.path.join(
                rgc.get_asset(args.genome_assembly, BT2_IDX_KEY),
                args.genome_assembly)
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

        if args.complexity and args.umi_len > 0:
            cmd2_dups = (tools.samtools + " view -q 10 -b -@ " + str(pm.cores) +
                " -U " + failQC_genome_bam_dups + " ")
            cmd2_dups += mapping_genome_bam_temp_dups + " > " + mapping_genome_bam_dups
            pm.clean_add(failQC_genome_bam_dups)

        def check_alignment_genome(temp_bam, bam):
            mr = ngstk.count_mapped_reads(temp_bam, args.paired_end)
            ar = ngstk.count_mapped_reads(bam, args.paired_end)
            if args.paired_end:
                ar = float(ar)/2
            rr = float(pm.get_stat("Raw_reads"))
            tr = float(pm.get_stat("Trimmed_reads"))
            if os.path.exists(res.refgene_pre_mRNA):
                cmd = (tools.samtools + " depth -b " +
                       res.refgene_pre_mRNA + " " + bam +
                       " | awk '{counter++;sum+=$3}END{print sum/counter}'")
                rd = pm.checkprint(cmd)
            else:
                cmd = (tools.samtools + " depth " + bam +
                       " | awk '{counter++;sum+=$3}END{print sum/counter}'")
                rd = pm.checkprint(cmd)
            pm.report_result("Mapped_reads", mr)
            pm.report_result("QC_filtered_reads",
                             round(float(mr)) - round(float(ar)))
            pm.report_result("Aligned_reads", ar)
            pm.report_result("Alignment_rate", round(float(ar) * 100 /
                             float(tr), 2))
            pm.report_result("Total_efficiency", round(float(ar) * 100 /
                             float(rr), 2))
            if rd and rd.strip():
                pm.report_result("Read_depth", round(float(rd), 2))

        pm.run([cmd, cmd2], mapping_genome_bam,
               follow=lambda: check_alignment_genome(mapping_genome_bam_temp,
                                                     mapping_genome_bam),
               container=pm.container)

        if args.complexity and args.umi_len > 0:
            pm.run([cmd_dups, cmd2_dups], mapping_genome_bam_dups,
                   container=pm.container)

        pm.timestamp("### Compress all unmapped read files")
        for unmapped_fq in to_compress:
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

            if args.complexity and args.umi_len > 0:
                cmd_dups = tools.samtools + " index " + mapping_genome_bam_temp_dups
                pm.run(cmd_dups, temp_mapping_index_dups)
                pm.clean_add(temp_mapping_index_dups)
                pm.clean_add(mapping_genome_bam_temp_dups)

        # Determine mitochondrial read counts
        if os.path.exists(mapping_genome_bam_temp):
            if not os.path.exists(temp_mapping_index):
                cmd = tools.samtools + " index " + mapping_genome_bam_temp
                pm.run(cmd, temp_mapping_index)
                pm.clean_add(temp_mapping_index)

            cmd = (tools.samtools + " idxstats " + mapping_genome_bam_temp +
                   " | grep")
            for name in mito_name:
                cmd += " -we '" + name + "'"
            cmd += "| cut -f 3"
            mr = pm.checkprint(cmd)
        
            # If there are mitochondrial reads, report and remove them
            if mr and mr.strip():
                pm.report_result("Mitochondrial_reads", round(float(mr)))
                # Index the sort'ed BAM file first
                mapping_genome_index = os.path.join(mapping_genome_bam + ".bai")
                noMT_mapping_genome_bam = os.path.join(
                    map_genome_folder, args.sample_name + "_noMT.bam")

                cmd1 = tools.samtools + " index " + mapping_genome_bam
                cmd2 = (tools.samtools + " idxstats " + mapping_genome_bam +
                        " | cut -f 1 | grep")
                for name in mito_name:
                    cmd2 += " -vwe '" + name + "'"
                cmd2 += ("| xargs " + tools.samtools + " view -b -@ " +
                         str(pm.cores) + " " + mapping_genome_bam + " > " +
                         noMT_mapping_genome_bam)
                cmd3 = ("mv " + noMT_mapping_genome_bam + " " + mapping_genome_bam)
                cmd4 = tools.samtools + " index " + mapping_genome_bam
                pm.run([cmd1, cmd2, cmd3, cmd4], noMT_mapping_genome_bam)
                pm.clean_add(mapping_genome_index)

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

    if int(args.max_len) != -1:
        max_len = args.max_len
    elif _itsa_file(mapping_genome_bam):
        cmd = (tools.samtools + " stats " + mapping_genome_bam +
               " | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-")
        max_len = int(pm.checkprint(cmd))
    else:
        max_len = DEFAULT_MAX_LEN

    # At this point we can check for seqOutBias required indicies.
    # Can't do it earlier because we haven't determined the read_length of 
    # interest for mappability purposes.
    if args.sob:
        pm.debug("max_len: {}".format(max_len))  # DEBUG
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
        res, rgc = _add_resources(args, res, search_asset)

    ############################################################################
    #                     Calculate library complexity                         #
    ############################################################################
    QC_folder = os.path.join(param.outfolder, "QC_" + args.genome_assembly)
    ngstk.make_dir(QC_folder)
    preseq_output = os.path.join(
        QC_folder, args.sample_name + "_preseq_out.txt")
    preseq_yield = os.path.join(
        QC_folder, args.sample_name + "_preseq_yield.txt")
    preseq_counts = os.path.join(
        QC_folder, args.sample_name + "_preseq_counts.txt")
    preseq_plot = os.path.join(
        QC_folder, args.sample_name + "_preseq_plot")
    preseq_pdf = os.path.join(
        QC_folder, args.sample_name + "_preseq_plot.pdf")
    preseq_png = os.path.join(
        QC_folder, args.sample_name + "_preseq_plot.png")

    if not _itsa_file(preseq_plot) or args.new_start:
        if args.complexity and args.umi_len > 0:
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
                    pm.run([cmd1, cmd2, cmd3, cmd4], noMT_mapping_genome_bam_dups)
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

                cmd = ("awk '{sum+=$2} END {printf \"%.0f\", sum}' " + res.chrom_sizes)
                genome_size = int(pm.checkprint(cmd))

                cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") +
                       " preseq " + "-i " + preseq_yield)
                if args.coverage:
                    cmd += (" -c " + str(genome_size) + " -l " + max_len)
                cmd += (" -r " + preseq_counts + " -o " + preseq_plot)

                pm.run(cmd, [preseq_pdf, preseq_png])

                pm.report_object("Library complexity", preseq_pdf,
                                 anchor_image=preseq_png)
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
    cmd += " --silent"

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

        pm.report_result("NRF", round(float(nrf),2))
        pm.report_result("PBC1", round(float(pbc1),2))
        pm.report_result("PBC2", round(float(pbc2), 2))

    pm.run(cmd, bamQC, follow=lambda: report_bam_qc(bamQC),
           container=pm.container)

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
    pm.run(unmap_cmd, unmap_genome_bam, follow=count_unmapped_reads,
           container=pm.container)

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

        if not pm.get_stat("TSS_Minus_Score") or args.new_start:
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
                # If the TSS enrichment is 0, don't report            
                Tss_score = (
                    (sum(floats[int(floats.index(max(floats))-49):
                                int(floats.index(max(floats))+51)]) / 100) /
                    (sum(floats[1:int(len(floats)*0.05)]) / int(len(floats)*0.05)))
                pm.report_result("TSS_Plus_Score", round(Tss_score, 1))
            except ZeroDivisionError:
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
                Tss_score = (
                    (sum(floats[int(floats.index(max(floats))-49):
                                int(floats.index(max(floats))+51)]) / 100) /
                    (sum(floats[1:int(len(floats)*0.05)]) / int(len(floats)*0.05)))
                pm.report_result("TSS_Minus_Score", round(Tss_score, 1))
            except ZeroDivisionError:
                pass

        # Call Rscript to plot TSS Enrichment
        TSS_pdf = os.path.join(QC_folder,  args.sample_name +
                               "_TSSenrichment.pdf")
        cmd = (tools.Rscript + " " + tool_path("PEPPRO.R"))
        cmd += " tss -i " + Tss_plus + " " + Tss_minus
        pm.run(cmd, TSS_pdf, nofail=True)

        TSS_png = os.path.join(QC_folder,  args.sample_name +
                               "_TSSenrichment.png")
        pm.report_object("TSS enrichment", TSS_pdf, anchor_image=TSS_png)

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
        pause_index = os.path.join(QC_folder, args.sample_name +
                                   "_pause_index.bed")
        if not pm.get_stat("Pause_index") or args.new_start:
            # Remove missing chr from PI annotations
            tss_local = os.path.join(QC_folder,
                                     args.genome_assembly + "_ensembl_tss.bed")
            body_local = os.path.join(QC_folder,
                                      args.genome_assembly + "_ensembl_gene_body.bed")
            cmd1 = ("grep -wf " + chr_keep + " " + res.ensembl_tss + " | " +
                    tools.bedtools + " sort -i stdin -faidx " + chr_order + 
                    " > " + tss_local)
            cmd2 = ("grep -wf " + chr_keep + " " + res.ensembl_gene_body + " | " +
                    tools.bedtools + " sort -i stdin -faidx " + chr_order + 
                    " > " + body_local)
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

            # Determine pause index
            cmd = ("join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 " +
                   TSS_density + " " + body_density +
                   " | awk -v OFS='\t' '{print $1, $2, $3, $4, ($6/($3-$2))" + 
                   "/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > " +
                   pause_index)
            pm.run(cmd, pause_index, nofail=True)

            # Median pause index
            cmd = ("sort -k5,5n " + pause_index +
                   " | awk ' { a[i++]=$5; } END " + 
                   "{ x=int((i+1)/2); if (x < (i+1)/2) " +
                   "print (a[x-1]+a[x])/2; else print a[x-1]; }'")
            val = pm.checkprint(cmd)
            if val and val.strip():
                pi = float(val)
                pm.report_result("Pause_index", round(pi, 2))

        # Plot pause index distribution
        pi_pdf = os.path.join(QC_folder, args.sample_name +
                              "_pause_index.pdf")
        pi_png = os.path.join(QC_folder, args.sample_name +
                              "_pause_index.png")
        cmd = (tools.Rscript + " " + tool_path("PEPPRO.R") + 
               " pi -i " + pause_index)
        pm.run(cmd, pi_pdf, nofail=True)
        pm.report_object("Pause index", pi_pdf, anchor_image=pi_png)

    ############################################################################
    #           Calculate Fraction of Reads in Pre-mRNA (FRiP)                 #
    ############################################################################
    if not os.path.exists(res.refgene_pre_mRNA):
        print("Skipping FRiP -- Fraction of reads in pre-mRNA requires "
              "pre-mRNA annotation file: {}"
              .format(res.refgene_pre_mRNA))
    else:
        pm.timestamp("### Calculate FRiP")
        # Plus
        plus_frip = calc_frip(plus_bam, res.refgene_pre_mRNA,
                              frip_func=ngstk.simple_frip,
                              pipeline_manager=pm)
        pm.report_result("Plus FRiP", round(plus_frip, 2))
        # Minus
        minus_frip = calc_frip(minus_bam, res.refgene_pre_mRNA,
                               frip_func=ngstk.simple_frip,
                               pipeline_manager=pm)
        pm.report_result("Minus FRiP", round(minus_frip, 2))

    ############################################################################
    #             Plot fragment distribution (for SE data)                     #
    ############################################################################
    if args.paired_end:
        pm.timestamp("### Plot fragment distribution")
        frag_len = os.path.join(QC_folder, args.sample_name + "_fragLen.txt")
        frag_dist_tool = tool_path("fragment_length_dist.pl")
        cmd = build_command([tools.perl, frag_dist_tool,
                             mapping_genome_bam, frag_len])

        fragL_counts_file = args.sample_name + "_fragCount.txt"
        fragL_count = os.path.join(QC_folder, fragL_counts_file)
        cmd1 = "sort -n  " + frag_len + " | uniq -c  > " + fragL_count

        fragL_dis1 = os.path.join(QC_folder, args.sample_name +
                                  "_fragLenDistribution.pdf")
        fragL_dis2 = os.path.join(QC_folder, args.sample_name +
                                  "_fragLenDistribution.txt")

        cmd2 = (tools.Rscript + " " + tool_path("PEPPRO.R"))
        cmd2 += " frag -l " + frag_len + " -c " + fragL_count
        cmd2 += " -p " + fragL_dis1 + " -t " + fragL_dis2

        pm.run([cmd, cmd1, cmd2], fragL_dis1, nofail=True)
        pm.clean_add(frag_len)
        pm.clean_add(fragL_count)

        fragL_png = os.path.join(QC_folder, args.sample_name +
                                 "_fragLenDistribution.png")
        pm.report_object("Fragment distribution", fragL_dis1,
                         anchor_image=fragL_png)
    else:
        pass
        # Used to plot adapter distribution here, but moved to cutadapt.

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

    ############################################################################ 
    #                  Determine genomic feature coverage                      #
    ############################################################################
    pm.timestamp("### Calculate fraction of reads in features (FRiF)")

    frif_plus_PDF = os.path.join(QC_folder, args.sample_name + "_plus_frif.pdf")
    frif_plus_PNG = os.path.join(QC_folder, args.sample_name + "_plus_frif.png")
    frif_minus_PDF = os.path.join(QC_folder,
                                  args.sample_name + "_minus_frif.pdf")
    frif_minus_PNG = os.path.join(QC_folder,
                                  args.sample_name + "_minus_frif.png")

    if not os.path.exists(frif_plus_PDF) or args.new_start:
        anno_list_plus = list()
        anno_list_minus = list()

        if os.path.isfile(anno_local):
            # Get list of features
            cmd1 = ("cut -f 4 " + anno_local + " | sort -u")
            ft_list = pm.checkprint(cmd1, shell=True)
            ft_list = ft_list.splitlines()

            # Split annotation file on features
            cmd2 = ("awk -F'\t' '{print>\"" + QC_folder + "/\"$4}' " +
                    anno_local)
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
                                                 valid_name + "_plus_coverage.bed")
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
                    cmd3 = ("cut -f 1 " + chr_order + " | grep -wf - " +
                            file_name + " | cut -f 1-3 | " +
                            "bedtools sort -i stdin -faidx " +
                            chr_order + " > " + anno_sort)
                    pm.run(cmd3, anno_sort)
                    
                    anno_list_plus.append(anno_cov_plus)
                    anno_list_minus.append(anno_cov_minus)
                    cmd4 = (tools.bedtools + " coverage -sorted -counts -a " +
                            anno_sort + " -b " + plus_bam +
                            " -g " + chr_order + " > " +
                            anno_cov_plus)
                    cmd5 = (tools.bedtools + " coverage -sorted -counts -a " +
                            anno_sort + " -b " + minus_bam +
                            " -g " + chr_order + " > " +
                            anno_cov_minus)
                    pm.run(cmd4, anno_cov_plus)
                    pm.run(cmd5, anno_cov_minus)

                    pm.clean_add(file_name)
                    pm.clean_add(anno_sort)
                    pm.clean_add(anno_cov_plus)
                    pm.clean_add(anno_cov_minus)

    ############################################################################
    #                                 Plot FRiF                                #
    ############################################################################
    pm.timestamp("### Plot FRiF")
    # Plus
    if not os.path.exists(frif_plus_PDF) or args.new_start:
        count_cmd = (tools.samtools + " view -@ " + str(pm.cores) + " " +
                     param.samtools.params + " -c -F4 " + plus_bam)
        plus_read_count = pm.checkprint(count_cmd)
        plus_read_count = str(plus_read_count).rstrip()

        frif_cmd = [tools.Rscript, tool_path("PEPPRO.R"), "frif",
                   "-n", args.sample_name, "-r", plus_read_count,
                   "-o", frif_plus_PDF, "--bed"]
        if anno_list_plus:
            for cov in anno_list_plus:
                frif_cmd.append(cov)
            cmd = build_command(frif_cmd)
            pm.run(cmd, frif_plus_PDF, nofail=False)
            pm.report_object("Plus FRiF", frif_plus_PDF,
                             anchor_image=frif_plus_PNG)

    # Minus
    if not os.path.exists(frif_minus_PDF) or args.new_start:
        count_cmd = (tools.samtools + " view -@ " + str(pm.cores) + " " +
                     param.samtools.params + " -c -F4 " + minus_bam)
        minus_read_count = pm.checkprint(count_cmd)
        minus_read_count = str(minus_read_count).rstrip()

        frif_cmd = [tools.Rscript, tool_path("PEPPRO.R"), "frif",
                   "-n", args.sample_name, "-r", minus_read_count,
                   "-o", frif_minus_PDF, "--bed"]
        if anno_list_minus:
            for cov in anno_list_minus:
                frif_cmd.append(cov)
            cmd = build_command(frif_cmd)
            pm.run(cmd, frif_minus_PDF, nofail=False)
            pm.report_object("Minus FRiF", frif_minus_PDF,
                             anchor_image=frif_minus_PNG)

    ############################################################################
    #                         Report mRNA contamination                        #
    ############################################################################
    if (os.path.exists(res.refgene_exon) and
        os.path.exists(res.refgene_intron)):

        pm.timestamp("### Calculate mRNA contamination")
        intron_exon = os.path.join(QC_folder, args.sample_name +
                                   "_exon_intron_ratios.bed")
        if not pm.get_stat("mRNA_contamination") or args.new_start:
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
                    chr_order + " | " + tools.bedtools + " sort -i stdin -faidx " +
                    chr_order + " > " + introns_sort)
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
            ar = float(pm.get_stat("Aligned_reads"))
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
                       "geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); " +
                       "gene[$4] = $4; " +
                       "chromEnd[$4]=$3; " +
                       "prev4=$4} END " +
                       "{ for (a in readCount) " +
                       "{ print chrom[a], chromStart[a], chromEnd[a], gene[a], " +
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
                       "geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); " +
                       "gene[$4] = $4; " +
                       "chromEnd[$4]=$3; " +
                       "prev4=$4} END " +
                       "{ for (a in readCount) " +
                       "{ print chrom[a], chromStart[a], chromEnd[a], gene[a], " +
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
                pm.report_result("mRNA_contamination", round(mrna_con, 2))

        # plot mRNA contamination distribution
        mRNApdf = os.path.join(QC_folder,
            args.sample_name + "_mRNA_contamination.pdf")
        mRNApng = os.path.join(QC_folder,
            args.sample_name + "_mRNA_contamination.png")
        mRNAplot = [tools.Rscript, tool_path("PEPPRO.R"), "mrna",
                    "-i", intron_exon, "--raw"]
        cmd = build_command(mRNAplot)
        pm.run(cmd, mRNApdf, nofail=False)
        pm.report_object("mRNA contamination", mRNApdf, anchor_image=mRNApng)

    ############################################################################
    #                        Shift and produce BigWigs                         #
    ############################################################################
    genome_fq = rgc.get_asset(args.genome_assembly,
                              asset_name="fasta",
                              seek_key="fasta")
    signal_folder = os.path.join(
        param.outfolder, "signal_" + args.genome_assembly)
    ngstk.make_dir(signal_folder)
    plus_bw = os.path.join(
        signal_folder, args.sample_name + "_plus_body_0-mer.bw")
    minus_bw = os.path.join(
        signal_folder, args.sample_name + "_minus_body_0-mer.bw")
    
    if not args.sob:
        # If not scaling we don't need to use seqOutBias to generate the
        # separate strand bigWigs; just convert the BAM's directly with 
        # bamSitesToWig.py which uses UCSC wigToBigWig
        pm.timestamp("### Produce bigWig files")

        wig_cmd_callable = ngstk.check_command("wigToBigWig")

        if wig_cmd_callable:
            cmd1 = tools.samtools + " index " + plus_bam
            cmd2 = tool_path("bamSitesToWig.py")
            cmd2 += " -i " + plus_bam
            cmd2 += " -c " + res.chrom_sizes
            cmd2 += " -o " + plus_bw  # DEBUG formerly smoothed " -w " + plus_bw
            cmd2 += " -p " + str(int(max(1, int(pm.cores) * 2/3)))
            cmd2 += " --variable-step"
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                cmd2 += " --tail-edge"
            pm.run([cmd1, cmd2], plus_bw)

            cmd3 = tools.samtools + " index " + minus_bam
            cmd4 = tool_path("bamSitesToWig.py")
            cmd4 += " -i " + minus_bam
            cmd4 += " -c " + res.chrom_sizes
            cmd4 += " -o " + minus_bw # DEBUG formerly smoothed " -w " + minus_bw
            cmd4 += " -p " + str(int(max(1, int(pm.cores) * 2/3)))
            cmd4 += " --variable-step"
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                cmd4 += " --tail-edge"
            pm.run([cmd3, cmd4], minus_bw)
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

        seqtable_cmd = build_command([
            (tools.seqoutbias, "seqtable"),
            res.fasta,
            str("--tallymer=" + res.search_file),
            str("--gt-workdir=" + tempdir),
            str("--read-size=" + str(max_len)),
            str("--out=" + seqtable)
        ])

        pm.run(seqtable_cmd, seqtable)
        pm.clean_add(seqtable)

        plus_table = os.path.join(
            signal_folder, (args.genome_assembly + "_plus_tbl.txt"))

        table_plus_cmd = build_command([
            (tools.seqoutbias, "table"),
            seqtable,
            plus_bam,
            (">", plus_table)
        ])

        minus_table = os.path.join(
            signal_folder, (args.genome_assembly + "_minus_tbl.txt"))

        table_minus_cmd = build_command([
            (tools.seqoutbias, "table"),
            seqtable,
            minus_bam,
            (">", minus_table)
        ])

        pm.run([table_plus_cmd, table_minus_cmd], minus_table)
        pm.clean_add(plus_table)
        pm.clean_add(minus_table)

        if args.scale:
            scale_plus_chunks = [
                (tools.seqoutbias, "scale"),
                seqtable,
                plus_bam,
                "--skip-bed",
                str("--bw=" + plus_bw)
            ]
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                scale_plus_chunks.extend([("--tail-edge")])
            scale_plus_cmd = build_command(scale_plus_chunks)

            scale_minus_chunks = [
                (tools.seqoutbias, "scale"),
                seqtable,
                minus_bam,
                "--skip-bed",
                str("--bw=" + minus_bw),
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
                str("--bw=" + plus_bw)
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
                str("--bw=" + minus_bw),
            ]
            if args.protocol.lower() in RUNON_SOURCE_PRO:
                scale_minus_chunks.extend([("--tail-edge")])
            scale_minus_cmd = build_command(scale_minus_chunks)

        pm.run([scale_plus_cmd, scale_minus_cmd], minus_bw)

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
