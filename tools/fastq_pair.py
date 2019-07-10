#!/usr/bin/env python
"""
An alternate method to pair fastq files

From: https://edwards.sdsu.edu/research/sorting-and-paring-fastq-files/
"""

from argparse import ArgumentParser
import os
import sys
from pararead import add_logging_options
from pararead import logger_via_cli
import gzip


def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quality scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rb')
    else:
        qin = open(fqfile, 'r')

    while True:
        header = qin.readline()
        if not header:
            break
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seq = qin.readline()
        seq = seq.strip()
        qualheader = qin.readline()
        qualscores = qin.readline()
        qualscores = qualscores.strip()
        header = header.replace('@', '', 1)
        yield seqid, header, seq, qualscores


def fast_pair(read1, read2, unmapped):
    _LOGGER.info("Read the first file into data structure...")
    if read1.endswith(('.gz', '.gzip')):
        name, ext = os.path.splitext(read1)
        r1_name, r1_ext = os.path.splitext(name)
    else:
        r1_name, r1_ext = os.path.splitext(read1)

    if read2.endswith(('.gz', '.gzip')):
        name, ext = os.path.splitext(read2)
        r2_name, r2_ext = os.path.splitext(name)
    else: 
        r2_name, r2_ext = os.path.splitext(read2)

    seqs = {}
    for (seqid, header, seq, qual) in stream_fastq(read1):
        seqid = seqid.replace('.1', '')
        seqs[seqid] = [header, seq, qual]

    if unmapped:
        r1_un = open(r1_name + "_unpaired.fastq", "w")
        r2_un = open(r2_name + "_unpaired.fastq", "w")

    r1_pe = open(r1_name + "_paired.fastq", "w")
    r2_pe = open(r2_name + "_paired.fastq", "w")

    _LOGGER.info("Read the second file into data structure...")
    seen = set()
    for (seqid, header, seq, qual) in stream_fastq(read2):
        seqid = seqid.replace('.2', '')
        seen.add(seqid)
        if seqid in seqs:
            r1_pe.write("@" + seqs[seqid][0] + "\n" + seqs[seqid][1] +
                        "\n+\n" + seqs[seqid][2] + "\n")
            r2_pe.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
        else:
            if unmapped:
                r2_un.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

    for seqid in seqs:
        if seqid not in seen and unmapped:
            r1_un.write("@" + seqs[seqid][0] + "\n" + seqs[seqid][1] +
                        "\n+\n" + seqs[seqid][2] + "\n")

    r1_pe.close()
    r2_pe.close()
    
    if unmapped:
        r1_un.close()
        r2_un.close()


def lowmem_pair(read1, read2, unmapped):
    _LOGGER.info("Read the first file into data structure...")
    if read1.endswith(('.gz', '.gzip')):
        name, ext = os.path.splitext(read1)
        r1_name, r1_ext = os.path.splitext(name)
    else:
        r1_name, r1_ext = os.path.splitext(read1)

    if read2.endswith(('.gz', '.gzip')):
        name, ext = os.path.splitext(read2)
        r2_name, r2_ext = os.path.splitext(name)
    else: 
        r2_name, r2_ext = os.path.splitext(read2)

    seqs = []
    for (seqid, header, seq, qual) in stream_fastq(read1):
        seqid = seqid.replace('.1', '')
        seqs.append(seqid)

    if unmapped:
        r2_un = open(r2_name + "_unpaired.fastq", "w")

    r2_pe = open(r2_name + "_paired.fastq", "w")

    _LOGGER.info("Read the second file into data structure...")
    seen = set()
    for (seqid, header, seq, qual) in stream_fastq(read2):
        seqid = seqid.replace('.2', '')
        if seqid in seqs:
            seqs.remove(seqid)
            r2_pe.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
        else:
            if unmapped:
                r2_un.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

    r2_pe.close()
    if unmapped:
        r2_un.close()
        r1_un = open(r1_name + "_unpaired.fastq", "w")

    r1_pe = open(r1_name + "_paired.fastq", "w")

    for (seqid, header, seq, qual) in stream_fastq(read1):
        seqid = seqid.replace('.1', '')
        if seqid not in seqs:
            # we have seen this so we deleted it
            r1_pe.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
        else:
            if unmapped:
                r1_un.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

    r1_pe.close()

    if unmapped:
        r1_un.close()


def parse_args(cmdl):
    parser = ArgumentParser(description="Re-pair fastq files.")
    parser.add_argument('-r1', '--read1', dest='read1', required=True,
                        help='Pair read1 file.')
    parser.add_argument('-r2', '--read2', dest='read2', required=True,
                        help='Pair read2 file.')
    parser.add_argument('-un', '--unmapped', dest='unmapped', required=False,
                        default=False, action='store_true',
                        help='Write the unpaired reads to separate files.')
    parser.add_argument('-w', '--lowmem', dest='lowmem', required=False,
                        default=False, action='store_true',
                        help='Conserve memory. Slow!.')
    parser = add_logging_options(parser)
    return parser.parse_args(cmdl)



if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    _LOGGER = logger_via_cli(args)

    if args.lowmem:
        lowmem_pair(args.read1, args.read2, args.unmapped)
    else:
        fast_pair(args.read1, args.read2, args.unmapped)
