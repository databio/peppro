#!/usr/bin/env bash
# Regenerate small test FASTQ files from the existing SE test read.
# Requires: seqtk, awk, gzip  (or run generate_test_data.py for a Python-only version)
# Run from the repository root: bash tests/scripts/generate_test_data.sh

set -euo pipefail

SRC="examples/data/test_r1.fq.gz"
OUTDIR="tests/data"

echo "Generating test data from ${SRC} ..."

# R1 copy (rename convention to fastq.gz)
cp "${SRC}" "${OUTDIR}/test_R1.fastq.gz"

# R2: reverse complement of R1
seqtk seq -r "${SRC}" | gzip > "${OUTDIR}/test_R2.fastq.gz"

# UMI R1: prepend 8-nt UMI (ACGTACGT / IIIIIIII) to every read
zcat "${SRC}" | awk '
    NR%4 == 1 { print; next }
    NR%4 == 2 { print "ACGTACGT" $0; next }
    NR%4 == 3 { print; next }
    NR%4 == 0 { print "IIIIIIII" $0 }
' | gzip > "${OUTDIR}/test_R1_umi.fastq.gz"

echo "Done. Files written to ${OUTDIR}/"
ls -lh "${OUTDIR}/"
