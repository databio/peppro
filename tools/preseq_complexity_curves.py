#!/usr/bin/env python
"""
preseq_complexity_curves.py

Takes results preseq and plots complexity curves.
"""


import sys
import os
import yaml
import glob
import subprocess
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from cycler import cycler
import numpy as np

def plot_complexity_curves(ccurves, coverage=0, read_length=0, real_counts_path=False, use_unique=True, output_name='complexity_curves', x_min=0, x_max=500000000):
    """
    This script plots the complexity curves generated for one or several libraries. The script is designed to work using
    the output produced by preseq (http://smithlabresearch.org/software/preseq/). preseq version 1.0.0 is currently
    supported by this script (the script is compatible also with version 0.1.0). Preseq is a tool used to estimate the
    library complexity and/or to predict the library complexity. In the first case "preseq c_curve" should be use. In
    the second case "preseq lc_extrap" should be usued. Please, refer to preseq manual available at
    http://smithlabresearch.org/wp-content/uploads/manual.pdf for examples (pages 12 to 14 are the most informatives ones)
    """

    if x_min < 0 or x_max <= x_min:
        sys.exit("problem with x-min or x-max ({}, {}). x-min must be equal or higher to 0 and less than x-max".format(x_min, x_max))

    # Covert limit counts to coverage
    if coverage > 0:
        if read_length == 0:
            raise RuntimeError ("Error: --coverage specified but not --read-length")
        else:
            coverage = float(coverage) / float(read_length)
        x_max = float(x_max) / coverage
        x_min = float(x_min) / coverage

    # Get the real counts if we have them
    real_counts_total = {}
    real_counts_unique = {}
    if real_counts_path is not False:
        real_counts_path = os.path.realpath(real_counts_path)
        try:
            with open(real_counts_path) as fh:
                for line in fh:
                    line = line.strip()
                    cols = line.split() # Split on any whitespace
                    for fn in [cols[0], "{}.preseq".format(cols[0]), "{}.bam".format(cols[0])]:
                        if fn in ccurves:
                            try:
                                if cols[1].isdigit():
                                    real_counts_total[fn] = int(cols[1])
                                if cols[2].isdigit() and use_unique:
                                    real_counts_unique[fn] = int(cols[2])
                            except IndexError:
                                pass
        except IOError as e:
            print("Error loading real counts file: {}".format(real_counts_path))
            raise IOError(e)
        else:
            print("Found {} matching sets of counts from {}".format(len(real_counts_total), real_counts_path))

    # Covert real counts to coverage
    if coverage > 0:
        for f, c in real_counts_total.items():
            real_counts_total[f] = float(c) / coverage
        for f, c in real_counts_unique.items():
            real_counts_unique[f] = float(c) / coverage


    # Set up plot params
    legend = [[],[]]
    global_x_max_ccurve_limit = 0
    global_y_max_ccurve_limit = 0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    max_label_length = 0

    # Each ccurve will get a different color
    colormap = plt.cm.gist_ncar
    colors = [colormap(i) for i in np.linspace(0, 0.9, len(ccurves))]
    plt.gca().set_prop_cycle(cycler('color', colors))

    # Go through inputs and plot line
    for ccurve in ccurves:
        print("Processing {}".format(ccurve))
        ccurve_table             = pd.io.parsers.read_csv(ccurve, sep='\t', header=0)
        ccurve_TOTAL_READS       = []
        ccurve_EXPECTED_DISTINCT = []
        if "TOTAL_READS" in ccurve_table:
            ccurve_TOTAL_READS       = ccurve_table["TOTAL_READS"].tolist()
            ccurve_EXPECTED_DISTINCT = ccurve_table["EXPECTED_DISTINCT"].tolist()
        elif "total_reads" in ccurve_table:
            ccurve_TOTAL_READS       = ccurve_table["total_reads"].tolist()
            ccurve_EXPECTED_DISTINCT = ccurve_table["distinct_reads"].tolist()
        else:
            sys.exit("Error, table {} is not in the expected format... has been generated with preseq?".format(ccurve))

        # Covert real counts to coverage
        if coverage > 0:
            ccurve_TOTAL_READS = [float(c) / coverage for c in ccurve_TOTAL_READS]
            ccurve_EXPECTED_DISTINCT = [float(c) / coverage for c in ccurve_EXPECTED_DISTINCT]

        # I need to find the interpolation point to print the plots
        x_min_ccurve_limit = computeLimit(x_min, ccurve_TOTAL_READS)
        x_max_ccurve_limit = computeLimit(x_max, ccurve_TOTAL_READS)
        try:
            if x_max_ccurve_limit > global_x_max_ccurve_limit:
                global_x_max_ccurve_limit = x_max_ccurve_limit
            if ccurve_EXPECTED_DISTINCT[x_max_ccurve_limit] > global_y_max_ccurve_limit:
                global_y_max_ccurve_limit = ccurve_EXPECTED_DISTINCT[x_max_ccurve_limit]
            x_max_ccurve_limit += 3 # Add a few points to be sure
        except IndexError:
            x_max_ccurve_limit = len(ccurve_EXPECTED_DISTINCT)

        # Plot the curve
        p, = ax.plot(ccurve_TOTAL_READS[x_min_ccurve_limit:x_max_ccurve_limit], ccurve_EXPECTED_DISTINCT[x_min_ccurve_limit:x_max_ccurve_limit])

        # Get the information for the legend
        sample_name = os.path.splitext(ccurve)[0]
        sample_name_raw = sample_name
        if(sample_name[:32] == 'accepted_hits_sorted_dupRemoved_'):
            sample_name = sample_name[32:]
        if(sample_name[:14] == 'accepted_hits_'):
            sample_name = sample_name[14:]
        if(len(sample_name) > max_label_length):
            max_label_length = len(sample_name)
        legend[0].append(p)
        legend[1].append(sample_name)

        # Plot the real data if we have it
        real_key = False
        if sample_name in real_counts_total:
            real_key = sample_name
        if ccurve in real_counts_total:
            real_key = ccurve
        if sample_name_raw + ".bam" in real_counts_total:
            real_key = sample_name_raw + ".bam"
        if real_key:

            t_reads = float(real_counts_total[real_key])

            if real_key in real_counts_unique:
                u_reads = int(real_counts_unique[real_key])
                ax.plot(t_reads, u_reads, 'o', color=p.get_color())
                print("INFO: Found real counts for {} - Total: {}, Unique: {}".format(sample_name, t_reads, u_reads))
            else:
                xvalues = p.get_xdata()
                yvalues = p.get_ydata()
                if t_reads > max(xvalues):
                    print("WARNING: Total reads for {} ({}) > max preseq value ({}) - skipping this point..".format(sample_name, t_reads, max(xvalues)))
                else:
                    interp = np.interp(t_reads, xvalues, yvalues)
                    ax.plot(t_reads, interp, 'o', color=p.get_color())
                    print("INFO: Found real count for {} - Total: {:.2f} (preseq unique reads: {:.2f})".format(sample_name, t_reads, interp))

    # plot perfect library as dashed line
    ax.plot([0, x_max], [0, x_max], color='black', linestyle='--', linewidth=1)

    # Set the axis limits
    max_total = 0
    if len(real_counts_total) > 0:
        max_total = int(max(d for d in real_counts_total.values()))
    if x_max < max_total:
        print("WARNING: x-max value {} is less than max real data {}".format(x_max, max_total))
    plt.xlim(x_min, x_max)

    max_unique = 0
    if len(real_counts_unique) > 0:
        max_unique = int(max(d for d in real_counts_unique.values()))
        max_unique += max_unique * 0.1
    preseq_ymax = global_y_max_ccurve_limit
    preseq_ymax += global_y_max_ccurve_limit * 0.1

    default_ylim = 100000
    if coverage > 0:
        default_ylim = float(default_ylim) / coverage

    plt.ylim(default_ylim, max(preseq_ymax, max_unique))
    if preseq_ymax < max_unique:
        print("WARNING: y-max value changed from default {} to the max real data {}".format(int(preseq_ymax), max_unique))

    # label the axis
    plt.ylabel('Unique Molecules')
    plt.xlabel('Total Molecules (including duplicates)')
    plt.title("Complexity Curve: preseq")

    # Change labels if we're using coverage
    if coverage > 0:
        plt.ylabel('Unique Coverage')
        plt.xlabel('Total Coverage (including duplicates)')

    if len(real_counts_unique) > 0:
        plt.text(0.5, -0.15, 'Points show read count versus deduplicated read counts (externally calculated)',
            horizontalalignment='center', fontsize=8, transform = ax.transAxes)
    elif len(real_counts_total) > 0:
        plt.text(0.5, -0.15, 'Points show externally calculated read counts on the curves',
            horizontalalignment='center', fontsize=8, transform = ax.transAxes)

    # Sort out some of the nastier plotting defaults
    ax.tick_params(top=False, right=False, direction='out')

    # Move the subplot around to fit in the legend
    box = ax.get_position()
    ax.set_position([0.08, box.y0+0.05, (box.width * 0.78)-0.02, box.height-0.05])
    # Set the font size according to how big the legend text is
    font = {'size': 5}
    if len(legend[1]) <= 20 and max_label_length <= 45:
        font = {'size': 6}
    if len(legend[1]) <= 20 and max_label_length <= 30:
        font = {'size': 8}
    if len(legend[1]) <= 20 and max_label_length <= 10:
        font = {'size': 12}
    ax.legend(legend[0], legend[1],loc='center left', bbox_to_anchor=(1.01, 0.5), prop=font)

    # now save the plot
    png_fn = "{}.png".format(output_name)
    pdf_fn = "{}.pdf".format(output_name)
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)
    plt.close(fig)
    return 0


def computeLimit(value, ccurve_TOTAL_READS):
    """This function returns the index of ccurve_TOTAL_READS containing the closest value to x_max"""
    if ccurve_TOTAL_READS[-1] < value:
        print("WARNING: value is set to a value higher than the highest extrapolated point by preseq (value={}, ccurve_TOTAL_READS[-1]={}).".format(value, ccurve_TOTAL_READS[-1]))
    first_point = 0
    last_point  = len(ccurve_TOTAL_READS)
    iterations = 0
    while first_point != last_point:
        middle_point = int((first_point + last_point)/2)
        middle_value = ccurve_TOTAL_READS[middle_point]
        if middle_value == value or iterations >= 10000:
            return middle_point
        elif middle_value >= value:
            last_point = middle_point -1
        else:
            first_point = middle_point +1
        iterations += 1

    return first_point


if __name__ == '__main__':
    parser = argparse.ArgumentParser("plot_complexity_curves.py", description=plot_complexity_curves.__doc__)
    parser.add_argument('ccurves', metavar='<preseq file>', nargs='+',
        help="List of input files generated by preseq")
    parser.add_argument('-c', '--coverage', dest='coverage', type=int, default=0,
        help="Use coverage on axes instead of read counts. Enter number of base pairs of refernce. Human GRCh38.p2 = 3221487035, GRCh37 = 3137144693")
    parser.add_argument('-l', '--read-length', dest='read_length', type=int, default=0,
        help="Sequence read length, for use in coverage calculations")
    parser.add_argument('-r', '--real-counts', dest='real_counts_path', type=str, default=False,
        help="File name for file with three columns - preseq filename, total number reads, number of unique reads (unique optional, whitespace delimited)")
    parser.add_argument('-u', '--ignore_unique', dest='use_unique', action='store_false',
        help="Ignore any information about unique read counts found in --real-counts file.")
    parser.add_argument('-o', '--output-name', dest='output_name', type=str, default='complexity_curves',
        help="Output name (.png will be automatically added)")
    parser.add_argument('-m', '--x-min', dest='x_min', type=int, default=0,
        help="Lower x-limit (default 0)")
    parser.add_argument('-x', '--x-max', dest='x_max', type=int, default=500000000,
        help="Upper x-limit (default 500 million)")
    kwargs = vars(parser.parse_args())
    plot_complexity_curves(**kwargs)
