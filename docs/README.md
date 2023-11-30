# <img src="img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:80px; margin-top:10px; margin-bottom:-10px" align="left">  

<br clear="all">

[![PEP compatible](https://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

`PEPPRO` is a pipeline for nascent RNA sequencing data. It can process PRO-seq, GRO-seq, and ChRO-seq data and is optimized on unique features of nascent RNA to be fast and accurate. It performs variable-length UMI adapter removal, read deduplication, trimming, mapping, QC, and signal tracks (bigWig) for plus and minus strands using mappability-scaled or unscaled read counts. 

## Outputs

`PEPPRO` produces quality control plots, statistics, and data formats to set the stage for project-specific analysis. We have produced an [interactive display of the output folder structure](browse_output/), which includes:

- **HTML report**: an easily-navigable HTML report with pretty plots: [HTML summary report demo](files/examples/paper/PEPPRO_summary.html).
- **Stats**: An easily parsable stats file: [Summary statistics demo file](files/examples/paper/PEPPRO_stats_summary.tsv).
- **Processed data**: Several bigWig signal tracks (plus and minus stranded), with options to produce: smoothed signal; exact (nucleotide-resolution) RNA polymerase position signal; or nucleotide-resolution signal corrected for enzymatic sequence bias.

## User interface

`PEPPRO` is a python script that runs on the command line (See [usage](usage)). It can also read projects in [PEP format](https://pepkit.github.io/). This means that `PEPPRO` projects are also compatible with other PEP tools, and output can be conveniently read into `R` using [the `pepr` package](http://code.databio.org/pepr/) or into `Python` using the [`peppy` package](https://peppy.readthedocs.io/en/latest/). The pipeline itself is customizable, enabling a user to adjust individual command settings or even swap out specific software by editing a few lines of human readable configuration files.

## Availability

You can download the latest version from the [releases](https://github.com/databio/peppro/releases) page, or find changes listed in the [CHANGELOG](changelog).

## Citing

If you find PEPPRO useful in your research, please cite:

Smith, J.P., Dutta, A.B., Sathyan, K.M. et al. PEPPRO: quality control and processing of nascent RNA profiling data. Genome Biol 22, 155 (2021). [https://doi.org/10.1186/s13059-021-02349-4](https://doi.org/10.1186/s13059-021-02349-4)
