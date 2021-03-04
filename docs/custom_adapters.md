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
