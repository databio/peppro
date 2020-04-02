#!/usr/bin/env python3
"""
PEPPRO Collator - Run-on sequencing project-level pipeline 
"""

__author__ = ["Michal Stolarczyk"]
__email__ = "michal@virginia.edu"
__version__ = "0.0.1"

from argparse import ArgumentParser
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

    return os.path.join(os.path.dirname(os.path.dirname(__file__)),
                        "tools", tool_name)


def parse_arguments():
    """
    Creat parser instance and parse command-line arguments passed to the pipeline

    :return argparse.Namespace: parsed arguments namespace
    """
    parser = VersionInHelpParser(prog="PEPPRO collator", description='PEPPRO collator' , version=__version__)
    parser.add_argument("-c", "--config", help="Path to the project config file.", type=str)
    parser.add_argument("-n", "--name", help="Name of the project to use.", type=str)
    parser.add_argument("-o", "--output", help="Output dir path.", type=str)
    args = parser.parse_args()
    return args

args = parse_arguments()

outfolder = os.path.abspath(os.path.join(args.output, "summary"))


def main():
    pm = pypiper.PipelineManager(name="PEPPRO collator", outfolder=outfolder, version=__version__)

    cmd1 = "Rscript {} {}".format(tool_path("PEPPRO_complexity_curves.R"), args.config)
    cmd2 = "Rscript {} {}".format(tool_path("PEPPRO_counts.R"), args.config)

    pm.run(cmd1, target=os.path.join(outfolder, "{name}_libComplexity.pdf".format(name=args.name)))
    pm.run(cmd2, target=os.path.join(outfolder, "{name}_countData.csv".format(name=args.name)))

    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
