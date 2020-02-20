# Configure UMI settings

By default, the pipeline assumes there is *not* a UMI. In other words, the parameter `umi_len` is set to 0. See the [pipeline usage documentation](usage.md) for additional parameter settings.

## Specify a UMI length 

There are three approaches for specifying the `umi_len` parameter for your samples.

### 1: Pass the `--umi-len` parameter at the command line

If you're running `PEPPRO` at the command line for a single sample, you may specify the UMI length using the `--umi-len` argument.  
For example: 
```
./pipelines/peppro.py \
  --sample-name test \
  --genome hg38 \
  --input examples/data/test_r1.fq.gz \
  --single-or-paired single \
  --umi-len 8 \
  -O $HOME/peppro_example/
```

### 2: Pass the `--umi-len` parameter to the pipeline using `looper`

If you're running `PEPPRO` with `looper`, you can also pass any number of additional arguments to `looper` that will be automatically passed to the pipeline.  
For example:
```
looper run examples/meta/peppro_test.yaml -d \
  --compute slurm \
  --umi-len 8
```

In this case, `looper` will automatically pass the `--umi-len 8` argument to each sample in the `peppro_test.yaml` file.

### 3: Specify a `--umi-len` argument in the project configuration file

If you're using `looper` and you'd like to set the `--umi-len` for individual samples that is entirely possible with some customization to the configuration and annotation files.  For a real life example, check out the [`peppro_paper.yaml`](https://github.com/databio/ppqc/blob/master/peppro_paper.yaml) and [`peppro_paper.csv`](https://github.com/databio/ppqc/blob/master/peppro_paper.csv) project files.

Below we'll go over two examples of customization in the project configuration files.

#### 1: Set a universal `--umi_len` in the project configuration file

```
# Run test sample through PEPPRO
name: test

metadata:
  sample_annotation: "peppro_test.csv"
  output_dir: "$PROCESSED/peppro/peppro_test/"
  pipeline_interfaces: "$CODE/peppro/pipeline_interface.yaml" 
  
derived_columns: [read1]

data_sources:
  R1: "$CODE/peppro/examples/data/{sample_name}_r1.fq.gz"
  
implied_columns:
  organism:
    human:
      genome: hg38
      prealignments: human_rDNA
      adapter: cutadapt  # Default
      dedup: seqkit      # Default
      trimmer: seqtk     # Default
      protocol: pro      # Default
      umi_len: 8         # Custom --umi-len that will be passed to **all** samples
      max_len: -1        # Disable length trimming 

pipeline_args:
#  peppro.py:
#    "--prioritize": null # Default is FALSE. Pass flag to prioritize features by the order they appear in the feat_annotation asset when calculating FRiF/PRiF
#    "--sob": null        # Default is FALSE. Pass flag to use seqOutBias for signal track generation and to incorporate mappability
#    "--scale": null      # Default is FALSE. Pass flag to scale seqOutBias signal tracks
#    "--coverage": null   # Default is FALSE. Pass flag to use coverage when producing library complexity plots.
#    "--keep": null       # Default is FALSE. Pass flag to keep prealignment BAM files.
#    "--noFIFO": null     # Default is FALSE. Pass flag to NOT use named pipes during prealignments.
#    "--complexity": null # Default is TRUE.  Pass flag to disable library complexity calculation. Faster.
```

#### 2: Set custom `--umi_len` arguments for individual samples

What do you do if your project doesn't have the same UMI for all samples? This requires a bit more complexity.  Our paper's samples do exactly this, and you can [go check out those configuration files specifically for complete detail](https://github.com/databio/ppqc/).  Here we'll highlight the relevant components.

In our [paper' project configuration file](https://github.com/databio/ppqc/blob/master/peppro_paper.yaml) we use the `implied_columns` feature to create custom arguments based on our [corresponding annotation file](https://github.com/databio/ppqc/blob/master/peppro_paper.csv).  We can include any column of choice, in this case we call it `umi_status`, that includes multiple possible matches. If a sample includes `true_8` in the `umi_status` column, we will pass a `--umi-len 8` to the pipeline.  Conversely, if a sample includes `true_6` in the `umi_status` column, we will pass a `--umi-len 6` to the pipeline. If neither argument is present, we use the default setting of length 0.

Here's a snippet of the relevant portion of the configuration files:

- configuration file
```
 umi_status:
    true_8:
      umi_len: 8
    true_6:
      umi_len: 6
```
- annotation file
```
umi_status, umi_length
FALSE, 0
FALSE, 0
true_6, 6
true_8, 8
true_8, 8
```

