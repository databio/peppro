# Running the pipeline script directly

It's easiest to run PEPPRO using `looper`, as described in the [install](install.md) and [tutorial](tutorial.md) guides. This greatly simplifies the process of running the pipeline across many sampels in a project. But really, the pipeline at its core is just a python script, and you can run it on the command line for a single sample (see [command-line usage](usage.md) by running `pipelines/peppro.py --help`). You just need to pass a few command-line parameters to specify sample name, reference genome, input files, etc. Here's the basic command to run the included small test example through the pipeline:

```console
cd peppro
./pipelines/peppro.py \
  --sample-name test \
  --genome hg38 \
  --input examples/data/test_r1.fq.gz \
  --single-or-paired single \
  -O $HOME/peppro_example/
```

This test example takes less than 5 minutes to complete.

