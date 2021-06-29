#!/usr/bin/env python3
"""
PEPPRO Collator - Run-on sequencing project-level pipeline 
"""

__author__ = ["Michal Stolarczyk", "Jason Smith"]
__email__ = "jasonsmith@virginia.edu"
__version__ = "0.0.3"

import os
import sys
import pypiper
from ubiquerg import VersionInHelpParser


def tool_path(tool_name):
    """
    Return the path to a tool used by this pipeline.

    :param str tool_name: name of the tool (e.g., a script file_name)
    :return str: real, absolute path to tool (expansion and symlink resolution)
    """

    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "tools", tool_name)


def parse_arguments():
    """
    Creat parser instance and parse command-line arguments passed to the pipeline

    :return argparse.Namespace: parsed arguments namespace
    """
    parser = VersionInHelpParser(
        prog="PEPPRO ollator", description="PEPPRO collator", version=__version__
    )
    parser = pypiper.add_pypiper_args(
        parser,
        groups=["pypiper", "looper", "pipestat"],
    )
    parser.add_argument("-n", "--name", help="Name of the project to use.", type=str)
    parser.add_argument(
        "-r", "--results", help="Output results sub directory path.", type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    outfolder = os.path.abspath(os.path.join(args.output_parent, "summary"))

    pm = pypiper.PipelineManager(
        name="PEPPRO_collator",
        outfolder=outfolder,
        args=args,
        version=__version__,
        pipestat_schema=args.pipestat_schema,
        pipestat_results_file=args.pipestat_results_file,
        pipestat_record_id=args.pipestat_record_id,
        pipestat_namespace=args.pipestat_namespace,
        pipestat_config=args.pipestat_config,
    )

    pm.pipestat.report(
        values={
            "log": {
                "path": pm.pipeline_log_file,
                "title": "Pipeline log file",
            },
            "profile": {
                "path": pm.pipeline_profile_file,
                "title": "Pipeline profile file",
            },
            "commands": {
                "path": pm.pipeline_commands_file,
                "title": "Pipeline commands file",
            },
            "version": pm.pl_version,
        }
    )

    # pm.info(f"{args=})

    cmd = "Rscript {R_file} {config_file} {output_dir} {results_subdir}".format(
        R_file=tool_path("PEPPRO_summarizer.R"),
        config_file=args.config_file,
        output_dir=args.output_parent,
        results_subdir=args.results,
    )
    if args.new_start:
        cmd += " --new-start"

    complexity_path = os.path.join(outfolder, f"{args.name}_libComplexity.pdf")
    complexity_thumbnail_path = os.path.join(
        outfolder, f"{args.name}_libComplexity.png"
    )
    counts_file = os.path.join(outfolder, f"{args.name}_countData.csv")
    pm.run(cmd, [complexity_path, counts_file])
    to_report = {}
    if all([os.path.exists(x) for x in [complexity_thumbnail_path, complexity_path]]):
        to_report.update(
            {
                "library_complexity_file": {
                    "path": complexity_path,
                    "thumbnail_path": complexity_thumbnail_path,
                    "title": "Library complexity file",
                }
            }
        )
    if os.path.exists(counts_file):
        to_report.update(
            {
                "counts_table": {
                    "path": counts_file,
                    "title": "Gene counts table",
                }
            }
        )
    if to_report:
        pm.pipestat.report(values=to_report)
    pm.stop_pipeline()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
