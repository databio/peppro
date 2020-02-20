#!/usr/bin/env python
"""
Package PEPPRO reports
"""

__author__ = ["Jason Smith"]
__email__ = "jasonsmith@virginia.edu"
__version__ = "0.0.1"


from argparse import ArgumentParser
import os
import sys
import tempfile
import tarfile
import shutil
import glob
import peppy

################################################################################

def parse_arguments():
    """
    Parse command-line arguments passed to the pipeline.
    """
    
    # Argument Parsing from yaml file
    parser = ArgumentParser(description='Package PEPPRO reports version ' + __version__)

    # Pipeline-specific arguments
    parser.add_argument("-p", "--project", dest="project",
                        help="PEP configuration file (*.yaml).")

    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {v}".format(v=__version__))

    args = parser.parse_args()

    if not args.project:
        parser.print_help()
        raise SystemExit

    return args

################################################################################

def main():
    """
    Main pipeline process.
    """

    args = parse_arguments()

    project = peppy.Project(args.project)

    # convenience alias
    projectName = project.name
    sampleNames = []
    for sample in project.samples:
        sampleNames.append(sample.name)
    projectDir = project.metadata.output_dir
    outputName = os.path.join(projectDir, projectName + ".tar.gz")

    tempdir = tempfile.mkdtemp(dir=projectDir)
    os.chmod(tempdir, 0o771)
    tempName = os.path.basename(tempdir)
    exclude = set([tempdir, os.path.join(tempName, tempName)])

    print("-- Copy project directory structure --")
    for (path, dirs, files) in os.walk(projectDir, topdown=True):
        if path not in exclude:
            newDir = path.replace(projectDir, (tempdir + "/"))
            if os.path.join(tempName, tempName) not in newDir:
                os.makedirs(newDir, exist_ok=True)

    print("-- Copy project reports --")
    projectReports = ['PEPPRO_stats_summary.tsv', 'PEPPRO_summary.html',
                      'PEPPRO_objs_summary.tsv', 'reports/*.html',
                      'summary/*.p[dn][fg]']
    for item in projectReports:
        dirName = os.path.dirname(item)
        path = os.path.join(tempdir, dirName)
        report = os.path.join(projectDir, item.strip('.'))
        for file in glob.glob(report):
            if os.path.exists(file) and 'tmp' not in str(file):
                shutil.copy2(file, path.strip('.'))

    print("-- Copy sample reports --")
    sampleReports = ['objects.tsv', 'stats.tsv', 'PEPPRO_log.md',
                     'PEPPRO_commands.sh', 'PEPPRO_profile.tsv',
                     'fastqc/*.html', 'QC_hg38/*.p[dn][fg]', 
                     'cutadapt/*.p[dn][fg]', 'fastp/*.p[dn][fg]',
                     'fastp/*.html']
    for item in sampleReports:
        dirName = os.path.dirname(item)
        for sample in sampleNames:
            sampleDir = os.path.join(projectDir, "results_pipeline", sample)
            if os.path.isdir(sampleDir):
                report = os.path.join(sampleDir, item.strip('.'))
                path = os.path.join(tempdir, "results_pipeline", sample, dirName)
                for file in glob.glob(report):
                    if os.path.exists(file) and 'tmp' not in str(file):
                        shutil.copy2(file, path)

    print("-- Archive report --")
    with tarfile.open(outputName, "w:gz") as tar:
        tar.add(tempdir, arcname='.')

    print("-- Clean up workspace --")
    shutil.rmtree(tempdir)

    print("-- Package project report complete! --")
    print("Report pacakge: {}".format(outputName))

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("PEPPRO report packager aborted.")
        sys.exit(1)