# Using custom adapters

A custom adapter file can be set in the pipeline configuration file, [`pipelines/peppro.yaml`](https://github.com/databio/peppro/blob/master/pipelines/peppro.yaml).  If left `null`, the default, then it uses the [included adapter file](https://github.com/databio/peppro/blob/master/tools/adapter.fa).

Alternatively, you may specify a path to your own adapter file.
For example: `adapters: /home/userid/my_custom_adapters.fa`

The only requirement here is that it expects specific headers for 5' and 3' adapters.
For example:
```
>5prime
TGGAATTCTCGGGTGCCAAGG
>3prime
GATCGTCGGACTGTAGAACTCTGAAC
```

## Using a custom pipeline configuration file with looper

If you want to use a custom `peppro.yaml` without modifying the original, specify the path via the `config_file` attribute in your project configuration file or sample sheet. This is passed to the pipeline as the `-C` argument.

### Apply to all samples (project-wide)

Add `config_file` under `sample_modifiers.append` in your project config:

```yaml
sample_modifiers:
  append:
    config_file: /path/to/my_peppro.yaml
```

### Apply to individual samples

Add a `config_file` column to your sample sheet CSV:

```
sample_name,genome,protocol,read1,config_file
my_sample,hg38,PRO,/path/to/reads.fq.gz,/path/to/my_peppro.yaml
```
