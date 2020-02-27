# FAQ

## Do I have to use PEPPRO with looper?

No. `PEPPRO` by itself does not specify any cluster resources, so you could just roll your own and submit individual jobs to a cluster however you choose. But because `PEPPRO` is already `looper`-compatible, the easier way is to use `looper's` built-in template system, which `looper` uses to build flexible shell scripts for job submission. These templates can be used to run jobs in a container, to submit to a cluster resource manager, or both.

## Will PEPPRO run on a mac?

The pipeline has been successfully run in both a Linux and MacOS environment.

## Can I run the project used in the paper?

Of course! Here is the exact [project configuration file](https://github.com/databio/ppqc/blob/master/peppro_paper.yaml) and [project annotation file](https://github.com/databio/ppqc/blob/master/peppro_paper.csv).

To run using `looper`:
```
looper run peppro_paper.yaml
```

## How do I determine if my samples are good quality?

To assist in determining whether you are working with a good quality nascent RNA library, we recommend several cutoffs for `PEPPRO` quality metrics that can be used to help make this decision.  If your sample meets or exceeds these values, it suggests it is a high quality sample.

| Metric                        | Recommended value |
|-------------------------------|-------------------|
| Degradation ratio             | < 1               |
| rDNA alignment rate           | < 20%             |
| Pause index                   | > 10              |
| mRNA contamination            | 1 - 1.5           |
| % uninformative adapter reads | < 25%             |
| TSS enrichment (coding)       | > 10              |
| TSS enrichment (non-coding)   | > 5               |
| % unique at 10M reads         | > 75%             |