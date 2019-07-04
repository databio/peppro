#!/usr/bin/env python

"""
Splits a paired-end fastq file into 2 files, one for each read. It can handle 2
formats: the reads are specified as with /1 or /2 in the read name; or with the
other notation: "1:N:0:NNNNNNNN" vs "2:N:0:NNNNNNNN" (which essentially splits
on "whitespace#")

By default it uses a named pipe which offers a balance between speed and memory.
The default, named-pipe method uses 2 processors.
You can opt for an in-memory solution (a direct pipe) which is faster, or a
single-processor version which is slowest and uses no memory (--gzip)
"""

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"
__credits__ = []
__license__ = "GPL3"
__version__ = "0.1"
__status__ = "Development"

from argparse import ArgumentParser
import gzip
import re
import os

parser = ArgumentParser(description='Fastq split')

# Add any pipeline-specific arguments
parser.add_argument('-p', '--prefix', dest='prefix', default=None,
	help="Read names are identified with a 1 or 2; What string is used as a prefix to identify 1 or 2? Default: autodetect.",
	required=False)

parser.add_argument('-i', '--infile', dest='infile',
		help="Input file (in fastq or fastq.gz)", required=True)

parser.add_argument('-o', '--outfile', dest='outfile',
		help="Output file basename (optional)", required=False)

parser.add_argument('--gzip', dest='gzip',
		help="Use gzip library for reading (slower, but uses no memory; uses only a single processor) (Default:False)",
		default=False, action="store_true", required=False)

parser.add_argument('--mem', dest='mem',
		help="Don't open a FIFO; use memory (fastest speed, but uses memory; uses 2 processors) (Default:False)",
		default=False, action="store_true", required=False)

parser.add_argument('--temp-dir', dest='temp_dir',
		help="Where to store temporary FIFO? (Default:/tmp)",
		default="/tmp", required=False)


args = parser.parse_args()

if args.infile.endswith(".fastq"):
	fin = open(args.infile, "r")
	fin_head = open(args.infile, "r")
elif args.infile.endswith(".gz"):
	if args.gzip:
		fin = gzip.open(args.infile, "r")
		fin_head = gzip.open(args.infile, "r")
	elif args.mem:
		# zcat pipe is 4 times faster than gzip.open (?)
		# Credits: Codebright (https://codebright.wordpress.com/2011/03/25/139/)
		import sys
		import subprocess
		if sys.version.startswith("3"):
			import io
			io_method = io.BytesIO
		else:
			import cStringIO
			io_method = cStringIO.StringIO

		p = subprocess.Popen(["zcat", args.infile], stdout = subprocess.PIPE)
		fin = io_method(p.communicate()[0])
		fin_head = io_method(p.communicate()[0])
	else:
		# Another way to do this, could be to use a named pipe (FIFO),
		# but this looks slower than the zcat direct method below.
		import io
		import subprocess
		import tempfile
		import shutil # to remove temporary folder
		import atexit
		prefix = "tmp_fqsplit_"
		tmp_fifo_dir = tempfile.mkdtemp(prefix=prefix, dir=args.temp_dir)
		tmp_fifo = os.path.join(tmp_fifo_dir, "fifo")
		def exit_handler():
			shutil.rmtree(tmp_fifo_dir)
		atexit.register(exit_handler)
		
		os.mkfifo(tmp_fifo)

		p = subprocess.Popen("gzip --stdout -d " + args.infile + " > %s" % tmp_fifo, shell=True)
		fin = io.open(tmp_fifo, "r")
		fin_head = io.open(tmp_fifo, "r")

# Which regex code should we use?
slash = 0  # \1 or \2 identifiers
space = 0  # " 1:N" or " 2:N" identifiers
# Test for \1 or " 1:" style notation
for i, line in enumerate(fin_head):
	if i % 4 == 0:
		if re.search("/1", line) or re.search("/2", line):
			slash += 1
		elif re.search(" 1", line) or re.search(" 2", line):
			space += 1
	if i * 4 > 11:
		break

fin_head.seek(0)  # Reset file pointer to beginning
fin_head.close()
 
if space > slash:
	args.prefix = " "
elif slash > space:
	args.prefix = "/"
else:
	raise Exception("Unknown regex to split reads.")

print("Identified regex split code: '" + args.prefix + "'")
if args.outfile == "" or args.outfile is None:
	outfile = os.path.basename(args.infile).split('.')[0]
else:
	outfile = args.outfile

fastq_file1 = outfile + "_R1.fastq"
fastq_file2 = outfile + "_R2.fastq"
fastq_1 = open(fastq_file1, "w")
fastq_2 = open(fastq_file2, "w")
print("Outfiles: '" + fastq_file1 + "' and '" + fastq_file2 + "'")

for line in fin:
	if re.search(args.prefix + "1", line):
		# Write 4 lines.
		fastq_1.write(line)
		fastq_1.write(next(fin))
		fastq_1.write(next(fin))
		fastq_1.write(next(fin))
	elif re.search(args.prefix + "2", line):
		fastq_2.write(line)
		fastq_2.write(next(fin))
		fastq_2.write(next(fin))
		fastq_2.write(next(fin))

fastq_1.close()
fastq_2.close()
fin.close()