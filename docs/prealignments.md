# Use prealignments with <img src="../../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px">

One feature of the pipeline is *prealignments*, which siphons off reads by aligning to small genomes before the main alignment to the primary reference genome.

Ideas for common prealignment references are provided by [ref_decoy](https://github.com/databio/ref_decoy).

Prealignments can be added via the command-line or included in a *PEP* project configuration file.

## Using prealignments on the command-line

In this example, we'll align sequentially to human mitochondrial sequence (e.g. the Revised Cambridge Reference Sequence of Human Mitochondrial DNA) and then a "human_repeats" genome which is the combination of alu elements, centromeric or alpha-satellite DNA, and ribosomal sequence.

```console
/pipelines/peppro.py \
  --sample-name test \
  --genome hg38 \
  --prealignments human_rDNA rCRSd \
  --input examples/data/test_r1.fq.gz \
  --single-or-paired single \
  -O $HOME/peppro_example/

```

## Adding prealignments to a project configuration file

See the included [`peppro_test.yaml`](https://github.com/databio/peppro/tree/master/examples/meta/peppro_test.yaml) for a simple example of setting `prealignments` in a project configuration file.
