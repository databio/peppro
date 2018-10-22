#!/usr/bin/env python
"""
PEPPRO - PRO-seq pipeline
"""

__author__ = ["Jason Smith", "Nathan Sheffield", "Mike Guertin"]
__email__ = "jasonsmith@virginia.edu"
__version__ = "0.1.0"


from argparse import ArgumentParser
import os
import sys
import re
import tempfile
import tarfile
import pypiper
from pypiper import build_command

TOOLS_FOLDER = "tools"
DEDUPLICATORS = ["fqdedup", "seqkit"]
ADAPTER_REMOVAL = ["cutadapt", "fastp"]
TRIMMERS = ["fastx", "seqtk"]

# TODO: if the input is PE, have to merge the files after the prep process
# fastp can handle both as input, will automatically interleave, then I need to de-interlace...
# actually I think bowtie2 can handle interleaved just fine
# cutadapt can use --interleaved option (Yes it works)
# I can deinterleave using the "deinterleave.sh" script

# TODO: can I skip the fastq dedup and do that at alignment? how much time do you lose
# because of that?

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
    parser.add_argument("--keep", action='store_true',
                        dest="keep",
                        help="Keep prealignment BAM files")
                    
    parser.add_argument("--noFIFO", action='store_true',
                        dest="no_fifo",
                        help="Do NOT use named pipes during prealignments")

    parser.add_argument("--umi", action='store_true',
                        dest="umi",
                        help="Remove umi with fastp")

    parser.add_argument("--adapter", dest="adapter",
                        default="fastp", choices=ADAPTER_REMOVAL,
                        help="Name of adapter removal program")

    parser.add_argument("--trimmer", dest="trimmer",
                        default="seqtk", choices=TRIMMERS,
                        help="Name of read trimming program")

    parser.add_argument("--dedup", dest="dedup",
                        default="seqkit", choices=DEDUPLICATORS,
                        help="Name of program that removes duplicate reads")

    parser.add_argument("--prealignments", default=[], type=str, nargs="+",
                        help="Space-delimited list of reference genomes to "
                             "align to before primary alignment.")

    parser.add_argument("-V", "--version", action="version",
                        version="%(prog)s {v}".format(v=__version__))

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        raise SystemExit

    return args


def _align_with_bt2(args, tools, paired, useFIFO, unmap_fq1, unmap_fq2,
                    assembly_identifier, assembly_bt2, outfolder,
                    aligndir=None, bt2_opts_txt=None):
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
        bamname = "{}_{}.bam".format(args.sample_name, assembly_identifier)
        mapped_bam = os.path.join(sub_outdir, bamname)
        summary_name = "{}_{}_bt_aln_summary.log".format(args.sample_name,
                                                          assembly_identifier)
        summary_file = os.path.join(sub_outdir, summary_name)
        out_fastq_pre = os.path.join(
            sub_outdir, args.sample_name + "_" + assembly_identifier)
        # bowtie2 unmapped filename format
        if paired:
            out_fastq_bt2 = out_fastq_pre + '_unmap_R%.fq.gz'
        else:
            out_fastq_bt2 = out_fastq_pre + '_unmap_R1.fq.gz'

        if not bt2_opts_txt:
            # Default options
            bt2_opts_txt = " -k 1"  # Return only 1 alignment
            bt2_opts_txt += " -D 20 -R 3 -N 1 -L 20 -i S,1,0.50"
            if paired and args.keep:
                bt2_opts_txt += " -X 2000"

        # samtools sort needs a temporary directory
        tempdir = tempfile.mkdtemp(dir=sub_outdir)
        pm.clean_add(tempdir)
   
        # Build bowtie2 command
        if args.keep:
            cmd = "(" + tools.bowtie2 + " -p " + str(pm.cores)
            cmd += bt2_opts_txt
            cmd += " -x " + assembly_bt2
            cmd += " --rg-id " + args.sample_name
            if paired:
                cmd += " -1 " + unmap_fq1 + " -2 " + unmap_fq2
                cmd += " --un-conc-gz " + out_fastq_bt2
            else:
                cmd += " -U " + unmap_fq1
                cmd += " --un-gz " + out_fastq_bt2
            cmd += " | " + tools.samtools + " view -bS - -@ 1"  # convert to bam
            cmd += " | " + tools.samtools + " sort - -@ 1"  # sort output
            cmd += " -T " + tempdir
            cmd += " -o " + mapped_bam
            cmd += ") 2>" + summary_file
            
        # In this samtools sort command we print to stdout and then use > to
        # redirect instead of  `+ " -o " + mapped_bam` because then samtools
        # uses a random temp file, so it won't choke if the job gets
        # interrupted and restarted at this step.            
        else:
            if useFIFO:
                out_fastq_tmp = os.path.join(sub_outdir,
                    assembly_identifier + "_bt2")
                if os.path.isfile(out_fastq_tmp):
                    out_fastq_tmp = os.path.join(sub_outdir,
                        assembly_identifier + "_bt2_2")
                cmd = "mkfifo " + out_fastq_tmp
                pm.run(cmd, out_fastq_tmp, container=pm.container)
            else:
                out_fastq_tmp = out_fastq_pre + '_unmap.fq'

            out_fastq_r1 = out_fastq_pre + '_unmap_R1.fq'
            out_fastq_r2 = out_fastq_pre + '_unmap_R2.fq'

            if paired:
                cmd1 = build_command([tools.perl,
                        tool_path("filter_paired_fq.pl"), out_fastq_tmp,
                        unmap_fq1, unmap_fq2, out_fastq_r1, out_fastq_r2])
            else:
                cmd1 = build_command([tools.perl,
                        tool_path("filter_paired_fq.pl"), out_fastq_tmp,
                        unmap_fq1, out_fastq_r1])

            cmd2 = "(" + tools.bowtie2 + " -p " + str(pm.cores)
            cmd2 += bt2_opts_txt
            cmd2 += " -x " + assembly_bt2
            cmd2 += " --rg-id " + args.sample_name
            cmd2 += " -U " + unmap_fq1
            cmd2 += " --un " + out_fastq_tmp
            cmd2 += " > /dev/null"
            cmd2 += ") 2>" + summary_file

        if args.keep:
            pm.run(cmd, mapped_bam, container=pm.container)
        else:
            if useFIFO and paired:
                pm.wait = False
                pm.run(cmd1, out_fastq_r2, container=pm.container)
                pm.wait = True
                pm.run(cmd2, summary_file, container=pm.container)
            else:
                pm.run(cmd2, summary_file, container=pm.container)
                pm.run(cmd1, out_fastq_r1, container=pm.container)

            pm.clean_add(out_fastq_tmp)
        
        # get concordant aligned read pairs
        if args.keep and paired:
            cmd = ("grep 'aligned concordantly exactly 1 time' " +
                   summary_file + " | awk '{print $1}'")
        else:
            cmd = ("grep 'aligned exactly 1 time' " +
                   summary_file + " | awk '{print $1}'")
        concordant = pm.checkprint(cmd)
        if concordant:
            ar = float(concordant)*2
        else:
            ar = 0

        # report concordant aligned reads
        pm.report_result("Aligned_reads_" + assembly_identifier, ar)
        try:
            # wrapped in try block in case Trimmed_reads is not reported in this
            # pipeline.
            tr = float(pm.get_stat("Trimmed_reads"))
        except:
            print("Trimmed reads is not reported.")
        else:
            res_key = "Alignment_rate_" + assembly_identifier
            pm.report_result(res_key, round(float(ar) * 100 / float(tr), 2))
        
        # filter genome reads not mapped
        if args.keep:
            unmap_fq1 = out_fastq_pre + "_unmap_R1.fq.gz"
            unmap_fq2 = out_fastq_pre + "_unmap_R2.fq.gz"
        else:
            unmap_fq1 = out_fastq_r1
            unmap_fq2 = out_fastq_r2

        return unmap_fq1, unmap_fq2
    else:
        msg = "No {} index found in {}; skipping.".format(
            assembly_identifier, os.path.dirname(assembly_bt2))
        print(msg)
        return unmap_fq1, unmap_fq2


def _count_alignment(assembly_identifier, aligned_bam, paired_end):
    """
    This function counts the aligned reads and alignment rate and reports
    statistics. You must have previously reported a "Trimmed_reads" result to
    get alignment rates. It is useful as a follow function after any alignment
    step to quantify and report the number of reads aligning, and the alignment
    rate to that reference.

    :param str  aligned_bam: Path to the aligned bam file.
    :param str assembly_identifier: String identifying the reference to which
                                    you aligned (can be anything)
    :param bool paired_end: Whether the sequencing employed a paired-end
                            strategy.
    """
	# count concordantly aligned reads ONLY
    ar = ngstk.count_concordant(aligned_bam)
	# Count all aligned reads
    #ar = ngstk.count_mapped_reads(aligned_bam, paired_end)
    pm.report_result("Aligned_reads_" + assembly_identifier, ar)
    try:
        # wrapped in try block in case Trimmed_reads is not reported in this
        # pipeline.
        tr = float(pm.get_stat("Trimmed_reads"))
    except:
        print("Trimmed reads is not reported.")
    else:
        res_key = "Alignment_rate_" + assembly_identifier
        pm.report_result(res_key, round(float(ar) * 100 / float(tr), 2))


def _get_bowtie2_index(genomes_folder, genome_assembly):
    """
    Create path to genome assembly folder with refgenie structure.

    Convenience function that returns the bowtie2 index prefix (to be passed
    to bowtie2) for a genome assembly that follows the folder structure
    produced by the RefGenie reference builder.

    :param str genomes_folder: path to central genomes directory, i.e. the
        root for multiple assembly subdirectories
    :param str genome_assembly: name of the specific assembly of interest,
        e.g. 'mm10'
    :return str: path to bowtie2 index subfolder within central assemblies
        home, for assembly indicated
    """   
    return os.path.join(genomes_folder, genome_assembly,
                        "indexed_bowtie2", genome_assembly)


def _check_bowtie2_index(genomes_folder, genome_assembly):
    """
    Confirm bowtie2 index is present.

    Checks by simple file count whether the bowtie2 index for a genome
    assembly (as produced by the RefGenie reference builder) contains the
    correct number of non-empty files.

    :param str genomes_folder: path to central genomes directory, i.e. the
        root for multiple assembly subdirectories
    :param str genome_assembly: name of the specific assembly of interest,
        e.g. 'mm10'
    """
    bt2_path = os.path.join(genomes_folder, genome_assembly, "indexed_bowtie2")
    
    if os.path.isdir(bt2_path):
        if not os.listdir(bt2_path):
            err_msg = "{} does not contain any files."
            pm.fail_pipeline(IOError(err_msg.format(bt2_path)))
        else:
            path, dirs, files = next(os.walk(bt2_path))
    elif os.path.isfile(os.path.join(genomes_folder, (genome_assembly + ".tar.gz"))):
        print("Did you mean this: {}".format(os.path.join(
            genomes_folder, (genome_assembly + ".tar.gz"))))
        err_msg = "Extract {} before proceeding."
        pm.fail_pipeline(IOError(err_msg.format(genome_assembly + ".tar.gz")))
    else:
        err_msg = "Could not find the {} index located at: {}"
        pm.fail_pipeline(IOError(err_msg.format(genome_assembly, bt2_path)))
    # check for bowtie small index
    if [bt for bt in files if bt.endswith('bt2')]:
        bt = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2',
              '.rev.1.bt2', '.rev.2.bt2']
    # check for bowtie large index
    elif [bt for bt in files if bt.endswith('bt2l')]:
        bt = ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l',
              '.rev.1.bt2l', '.rev.2.bt2l']
    # if neither file type present, fail
    else:
        err_msg = "{} does not contain any bowtie2 index files."
        pm.fail_pipeline(IOError(err_msg.format(bt2_path)))

    bt_expected = [genome_assembly + s for s in bt]
    bt_present  = [bt for bt in files if any(s in bt for s in bt_expected)]
    bt_missing  = list(set(bt_expected) - set(bt_present))
    # if there are any missing files (bowtie2 file naming is constant), fail
    if bt_missing:
        err_msg = "The {} bowtie2 index is missing the following file(s): {}"
        pm.fail_pipeline(IOError(
            err_msg.format(genome_assembly,
                           ', '.join(str(s) for s in bt_missing))))
    else:
        for f in bt_present:
            # If any bowtie2 files are empty, fail
            if os.stat(os.path.join(bt2_path, f)).st_size == 0:
                err_msg = "The bowtie2 index file, {}, is empty."
                pm.fail_pipeline(IOError(err_msg.format(f)))

    genome_file = genome_assembly + ".fa"
    fa_files = [fa for fa in files if genome_file in fa]
    if not fa_files:
        # The fasta file does not exist
        err_msg = "Could not find {}.fa in {}."
        pm.fail_pipeline(IOError(
            err_msg.format(genome_assembly, bt2_path)))
    for f in fa_files:
        if os.stat(os.path.join(bt2_path, f)).st_size == 0:
            pm.fail_pipeline(IOError("{} is an empty file.".format(f)))

def tool_path(tool_name):
    """
    Return the path to a tool used by this pipeline.

    :param str tool_name: name of the tool (e.g., a script filename)
    :return str: real, absolute path to tool (expansion and symlink resolution)
    """

    return os.path.join(os.path.dirname(os.path.dirname(__file__)),
                        TOOLS_FOLDER, tool_name)


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
    res = pm.config.resources

    # Set up reference resource according to genome prefix.
    gfolder = os.path.join(res.genomes, args.genome_assembly)
    res.chrom_sizes = os.path.join(
        gfolder, args.genome_assembly + ".chromSizes")

    # Get bowtie2 indexes
    res.bt2_genome = _get_bowtie2_index(res.genomes, args.genome_assembly)
    _check_bowtie2_index(res.genomes, args.genome_assembly)
    for reference in args.prealignments:
        _check_bowtie2_index(res.genomes, reference)

    # Adapter file can be set in the config; if left null, we use a default.
    res.adapters = res.adapters or tool_path("PRO-seq_adapter.fa")

    param.outfolder = outfolder

    print("Local input file: " + args.input[0])
    if args.input2:
        print("Local input file: " + args.input2[0])

    ###########################################################################

    pm.report_result(
        "File_mb",
        ngstk.get_file_size(
            [x for x in [args.input, args.input2] if x is not None]))
    pm.report_result("Read_type", args.single_or_paired)
    pm.report_result("Genome", args.genome_assembly)

    # PRO-seq pipeline
    # Each (major) step should have its own subfolder

    raw_folder = os.path.join(param.outfolder, "raw")
    fastq_folder = os.path.join(param.outfolder, "fastq")

    pm.timestamp("### Merge/link and fastq conversion: ")
    # This command will merge multiple inputs so you can use multiple
    # sequencing lanes in a single pipeline run.
    local_input_files = ngstk.merge_or_link(
        [args.input, args.input2], raw_folder, args.sample_name)
    cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(
        local_input_files, args.sample_name, args.paired_end, fastq_folder)
    pm.run(cmd, unaligned_fastq,
           follow=ngstk.check_fastq(
               local_input_files, unaligned_fastq, args.paired_end),
           container=pm.container)
    pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)
    print(local_input_files)
    untrimmed_fastq1 = out_fastq_pre + "_R1.fastq"
    untrimmed_fastq2 = out_fastq_pre + "_R2.fastq" if args.paired_end else None

    ########################
    # Begin fastq trimming
    ########################
    pm.timestamp("### FASTQ processing: ")

    # Create names for processed FASTQ files.
    processed_fastq = os.path.join(
        fastq_folder, args.sample_name + "_R1.trim.fastq")
    processed_fastq_R2 = os.path.join(
        fastq_folder, args.sample_name + "_R2.trim.fastq")
    fastqc_folder=os.path.join(param.outfolder, "fastqc")
    adapter_report = os.path.join(
        fastqc_folder, args.sample_name + "_rmAdapter.html")
    umi_report = os.path.join(
        fastqc_folder, args.sample_name + "_rmUmi.html")

    # Create adapter trimming command(s).
    if args.adapter == "fastp":
        if args.paired_end:
            adapter_cmd_chunks = [
                tools.fastp,
                ("--thread", str(pm.cores)),
                ("--in1", untrimmed_fastq1),
                ("--in2", untrimmed_fastq2),
                "--stdout",
                ("--adapter_sequence", "TGGAATTCTCGGGTGCCAAGG"),
                ("--length_required", 26),
                ("--html", adapter_report)
            ]
        else:
            adapter_cmd_chunks = [
                tools.fastp,
                ("--thread", str(pm.cores)),
                ("--in1", untrimmed_fastq1),
                "--stdout",
                ("--adapter_sequence", "TGGAATTCTCGGGTGCCAAGG"),
                ("--length_required", 26),
                ("--html", adapter_report)
            ]
        adapter_cmd = build_command(adapter_cmd_chunks)

    elif args.adapter == "cutadapt":
        if args.paired_end:
            adapter_cmd_chunks = [
                tools.cutadapt,
                ("-m", 26),
                "--interleaved",
                untrimmed_fastq1,
                untrimmed_fastq2
            ]
        else:
            adapter_cmd_chunks = [
                tools.cutadapt,
                ("-m", 26),
                untrimmed_fastq1
            ]
        adapter_cmd = build_command(adapter_cmd_chunks)

    else:
        # Default to fastp
        if args.paired_end:
            adapter_cmd_chunks = [
                tools.fastp,
                ("--thread", str(pm.cores)),
                ("--in1", untrimmed_fastq1),
                ("--in2", untrimmed_fastq2),
                "--stdout",
                ("--adapter_sequence", "TGGAATTCTCGGGTGCCAAGG"),
                ("--length_required", 26),
                ("--html", adapter_report)
            ]
        else:
            adapter_cmd_chunks = [
                tools.fastp,
                ("--thread", str(pm.cores)),
                ("--in1", untrimmed_fastq1),
                "--stdout",
                ("--adapter_sequence", "TGGAATTCTCGGGTGCCAAGG"),
                ("--length_required", 26),
                ("--html", adapter_report)
            ]
        adapter_cmd = build_command(adapter_cmd_chunks)
    
    # Create deduplication command(s).
    if args.dedup == "seqkit":
        dedup_cmd_chunks = [
            (tools.seqkit, "rmdup"),
            ("--threads", str(pm.cores)),
            "--by-seq",
            "--ignore-case",
            ("-o", "-")
        ]
        dedup_cmd = build_command(dedup_cmd_chunks)

    elif args.dedup == "fqdedup":
        dedup_cmd_chunks = [
            tools.fqdedup,
            ("-i", "-"),
            ("-o", "-")
        ]
        dedup_cmd = build_command(dedup_cmd_chunks)

    else:
        # Default to seqkit
        dedup_cmd_chunks = [
            (tools.seqkit, "rmdup"),
            ("--threads", str(pm.cores)),
            "--by-seq",
            "--ignore-case",
            ("-o", "-")
        ]
        dedup_cmd = build_command(dedup_cmd_chunks)

    # Create trimming and reverse complementing command(s).
    # TODO: Can also use seqkit for these steps instead of seqtk...
    if args.umi:
        if args.adapter != "fastp":
            print("To remove UMI intelligently, you must process your reads using 'fastp'")
            print("Defaulting to removing the first 8 bp instead via trimming")
            if args.trimmer == "seqtk":
                trim_cmd_chunks = [
                    (tools.seqtk, "trimfq"),
                    ("-b", 8),
                    ("-l", 30),
                    ("-", "|"),
                    (tools.seqtk, "seq"),
                    ("-r", "-"),
                    (">", processed_fastq)
                ]
                trim_cmd = build_command(trim_cmd_chunks)

            elif args.trimmer == "fastx":
                trim_tool = tools.fastx + "_trimmer"
                rc_tool = tools.fastx + "_reverse_complement"
                trim_cmd_chunks = [
                    trim_tool,
                    ("-Q", 33),
                    ("-f", 9),
                    ("-l", 38),
                    "|",
                    rc_tool,
                    ("-Q", 33),
                    ("-o", processed_fastq)
                ]
                trim_cmd = build_command(trim_cmd_chunks)

            else:
                # Default to seqtk
                trim_cmd_chunks = [
                    (tools.seqtk, "trimfq"),
                    ("-b", 8),
                    ("-l", 30),
                    "-",
                    "|",
                    (tools.seqtk, "seq"),
                    ("-r", "-"),
                    (">", processed_fastq)
                ]
                trim_cmd = build_command(trim_cmd_chunks)

        else:
            if args.paired_end:
                trim_cmd_chunks = [
                    tools.fastp,
                    ("--thread", str(pm.cores)),
                    ("--stdin", "--stdout"),
                    "--interleaved_in",
                    "--umi",
                    ("--umi_loc", "read1"),
                    ("--umi_len", 8),
                    ("--html", umi_report),
                    "|",
                    (tools.seqtk, "trimfq"),
                    ("-L", 30),
                    "-",
                    "|",
                    (tools.seqtk, "seq"),
                    ("-r", "-"),
                    (">", processed_fastq)
                ]

            else:
                trim_cmd_chunks = [
                    tools.fastp,
                    ("--thread", str(pm.cores)),
                    ("--stdin", "--stdout"),
                    "--umi",
                    ("--umi_loc", "read1"),
                    ("--umi_len", 8),
                    ("--html", umi_report),
                    "|",
                    (tools.seqtk, "trimfq"),
                    ("-L", 30),
                    "-",
                    "|",
                    (tools.seqtk, "seq"),
                    ("-r", "-"),
                    (">", processed_fastq)
                ]
            trim_cmd = build_command(trim_cmd_chunks)

    else:
        if args.trimmer == "seqtk":
            trim_cmd_chunks = [
                (tools.seqtk, "trimfq"),
                ("-b", 8),
                ("-l", 30),
                ("-", "|"),
                (tools.seqtk, "seq"),
                ("-r", "-"),
                (">", processed_fastq)
            ]
            trim_cmd = build_command(trim_cmd_chunks)

        elif args.trimmer == "fastx":
            trim_tool = tools.fastx + "_trimmer"
            rc_tool = tools.fastx + "_reverse_complement"
            trim_cmd_chunks = [
                trim_tool,
                ("-Q", 33),
                ("-f", 9),
                ("-l", 38),
                "|",
                rc_tool,
                ("-Q", 33),
                ("-o", processed_fastq)
            ]
            trim_cmd = build_command(trim_cmd_chunks)

        else:
            # Default to seqtk
            trim_cmd_chunks = [
                (tools.seqtk, "trimfq"),
                ("-b", 8),
                ("-l", 30),
                "-",
                "|",
                (tools.seqtk, "seq"),
                ("-r", "-"),
                (">", processed_fastq)
            ]
            trim_cmd = build_command(trim_cmd_chunks)
    
    # Put it all together
    process_fastq_cmd = build_command([
        adapter_cmd, "|", dedup_cmd, "|", trim_cmd])
    
    pm.run(process_fastq_cmd, processed_fastq)

    pm.clean_add(os.path.join(fastq_folder, "*.fq"), conditional=True)
    pm.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)
    #########################
    # End fastq processing
    #########################

    #########################
    # Deinterleave if paired-end
    #########################
    if args.paired_end:
        pm.timestamp("### Deinterleave processed FASTQ files")
        deinterleave_fq1 = os.path.join(
            fastq_folder, args.sample_name + "_R1.processed.fastq.gz")
        deinterleave_fq2 = os.path.join(
            fastq_folder, args.sample_name + "_R1.processed.fastq.gz")
        deinterleave_cmd_chunks = [
            tool_path("deinterleave.sh"),
            ("<", processed_fastq),
            (deinterleave_fq1, deinterleave_fq2),
            "compress"    
        ]
        cmd = build_command(deinterleave_cmd_chunks)
        pm.run(cmd, deinterleave_fq2)

        # Prepare variables for alignment step
        unmap_fq1 = deinterleave_fq1
        unmap_fq2 = deinterleave_fq2
    else:
        unmap_fq1 = processed_fastq

    # Map to any requested prealignments
    # We recommend mapping to chrM first for PRO-seq data
    pm.timestamp("### Prealignments")
    if len(args.prealignments) == 0:
        print("You may use `--prealignments` to align to references before "
              "the genome alignment step. See docs.")
    else:
        print("Prealignment assemblies: " + str(args.prealignments))
        # Loop through any prealignment references and map to them sequentially
        for reference in args.prealignments:
            if args.no_fifo:
                unmap_fq1, unmap_fq2 = _align_with_bt2(
                args, tools, args.paired_end, False, unmap_fq1, unmap_fq2, reference,
                assembly_bt2=_get_bowtie2_index(res.genomes, reference),
                outfolder=param.outfolder, aligndir="prealignments")
            else:
                unmap_fq1, unmap_fq2 = _align_with_bt2(
                args, tools, args.paired_end, True, unmap_fq1, unmap_fq2, reference,
                assembly_bt2=_get_bowtie2_index(res.genomes, reference),
                outfolder=param.outfolder, aligndir="prealignments")            

    pm.timestamp("### Map to genome")
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

    bt2_options = " --very-sensitive"
    bt2_options += " -X 2000"

    # samtools sort needs a temporary directory
    tempdir = tempfile.mkdtemp(dir=map_genome_folder)
    pm.clean_add(tempdir)

    cmd = tools.bowtie2 + " -p " + str(pm.cores)
    cmd += bt2_options
    cmd += " --rg-id " + args.sample_name
    cmd += " -x " + res.bt2_genome
    if args.paired_end:
        cmd += " -1 " + unmap_fq1 + " -2 " + unmap_fq2
    else:
        cmd += " -U " + unmap_fq1
    cmd += " | " + tools.samtools + " view -bS - -@ 1 "
    cmd += " | " + tools.samtools + " sort - -@ 1"
    cmd += " -T " + tempdir
    cmd += " -o " + mapping_genome_bam_temp

    # Split genome mapping result bamfile into two: high-quality aligned
    # reads (keepers) and unmapped reads (in case we want to analyze the
    # altogether unmapped reads)
    # -q 10: skip alignments with MAPQ less than 10
    cmd2 = (tools.samtools + " view -q 10 -b -@ " + str(pm.cores) +
            " -U " + failQC_genome_bam + " ")
    if args.paired_end:
        # add a step to accept only reads mapped in proper pair
        cmd2 += "-f 2 "

    cmd2 += mapping_genome_bam_temp + " > " + mapping_genome_bam

    def check_alignment_genome():
        mr = ngstk.count_mapped_reads(mapping_genome_bam_temp, args.paired_end)
        ar = ngstk.count_mapped_reads(mapping_genome_bam, args.paired_end)
        rr = float(pm.get_stat("Raw_reads"))
        tr = float(pm.get_stat("Trimmed_reads"))
        pm.report_result("Mapped_reads", mr)
        pm.report_result("QC_filtered_reads",
                         round(float(mr)) - round(float(ar)))
        pm.report_result("Aligned_reads", ar)
        pm.report_result("Alignment_rate", round(float(ar) * 100 /
                         float(tr), 2))
        pm.report_result("Total_efficiency", round(float(ar) * 100 /
                         float(rr), 2))

    pm.run([cmd, cmd2], mapping_genome_bam,
           follow=check_alignment_genome, container=pm.container)

    if not args.prealignments:
        # Index the temporary bam file
        temp_mapping_index = os.path.join(mapping_genome_bam_temp + ".bai")
        cmd = tools.samtools + " index " + mapping_genome_bam_temp
        pm.run(cmd, temp_mapping_index, container=pm.container)
        pm.clean_add(temp_mapping_index)

        # Determine mitochondrial read counts
        mito_name = ["chrM", "chrMT", "M", "MT"]
        cmd1 = (tools.samtools + " idxstats " + mapping_genome_bam_temp +
                " | grep")
        for name in mito_name:
            cmd1 += " -we '" + name + "'"
        cmd1 += "| cut -f 3"
        mr = pm.checkprint(cmd1)

        # If there are mitochondrial reads, report and remove them
        if mr and mr.strip():
            pm.report_result("Mitochondrial_reads", round(float(mr)))
            noMT_mapping_genome_bam = os.path.join(
                map_genome_folder, args.sample_name + "_noMT.bam")
            cmd2 = (tools.samtools + " idxstats " + mapping_genome_bam +
                    " | cut -f 1 | grep")
            for name in mito_name:
                cmd2 += " -vwe '" + name + "'"
            cmd2 += ("| xargs " + tools.samtools + " view -b -@ " +
                     str(pm.cores) + " " + mapping_genome_bam + " > " +
                     noMT_mapping_genome_bam)
            cmd3 = ("mv " + noMT_mapping_genome_bam + " " + mapping_genome_bam)
            cmd4 = tools.samtools + " index " + mapping_genome_bam
            pm.run([cmd2, cmd3, cmd4], noMT_mapping_genome_bam,
                   container=pm.container)

    # Calculate quality control metrics for the alignment file
    pm.timestamp("### Calculate NRF, PBC1, and PBC2")
    QC_folder = os.path.join(param.outfolder, "QC_" + args.genome_assembly)
    ngstk.make_dir(QC_folder)

    # Need index for mapping_genome_bam before calculating bamQC metrics
    mapping_genome_index = os.path.join(mapping_genome_bam + ".bai")
    cmd = tools.samtools + " index " + mapping_genome_bam
    pm.run(cmd, mapping_genome_index, container=pm.container)

    bamQC = os.path.join(QC_folder, args.sample_name + "_bamQC.tsv")
    cmd = tool_path("bamQC.py")
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

        pm.report_result("NRF", round(float(nrf),2))
        pm.report_result("PBC1", round(float(pbc1),2))
        pm.report_result("PBC2", round(float(pbc2), 2))

    pm.run(cmd, bamQC, follow=lambda: report_bam_qc(bamQC),
           container=pm.container)

    # Now produce the unmapped file
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

    # Separate by strand
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
        ("-f", 0x10),
        mapping_genome_bam,
        (">", minus_bam)
    ])
    
    pm.run([cmd1,cmd2], minus_bam)

    # Shift and produce BigWig's
    pm.timestamp("### Shift and generate BigWig files")
    genome_fq = os.path.join(
        genomes_folder, genome_assembly, (genome_assembly + ".fa"))
    signal_folder = os.path.join(
        param.outfolder, "signal_" + args.genome_assembly)
    plus_bw = os.path.join(
        signal_folder, args.sample_name + "_plus_body_0-mer.bw")
    minus_bw = os.path.join(
        signal_folder, args.sample_name + "_minus_body_0-mer.bw")
    
    cmd1 = build_command([
        tools.seqoutbias,
        genome_fq,
        plus_bam,
        "--no-scale",
        "--skip-bed",
        ("--bw=", plus_bw),
        "--tail-edge",
        ("--read-size=", 30),
    ])
    cmd2 = build_command([
        tools.seqoutbias,
        genome_fq,
        minus_bam,
        "--no-scale",
        "--skip-bed",
        ("--bw=", minus_bw),
        "--tail-edge",
        ("--read-size=", 30),
    ])
    
    pm.run([cmd1,cmd2], minus_bam)

    # COMPLETE!
    pm.stop_pipeline()


if __name__ == '__main__':
    pm = None
    ngstk = None
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
