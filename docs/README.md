# <img src="img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:200px; margin-top:10px; margin-bottom:-10px" align="left">  

<br clear="all">

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

`PEPPRO` is a pipeline designed to process PRO-seq data. It is optimized on unique features of PRO-seq to be fast and accurate. It performs adapter removal, including UMI of variable length, read deduplication, trimming, mapping, and signal tracks (bigWig) for plus and minus strands using scaled (based on mappability information) or unscaled read count patterns. 

## Outputs

`PEPPRO` produces quality control plots, summary statistics, and several data formats to set the stage for project-specific analysis. 

- PEPPRO produces an easily-navigable HTML report when used with [`Looper`](http://looper.databio.org/en/latest/): View this [HTML Summary report demo](files/examples/K562/summary.html)
- We have produced an [interactive display of the output folder structure](browse_output/), which includes:
	- [Easily parsable summary statistics file](files/examples/K562/results_pipeline/K562/stats.tsv)
	- BigWig signal tracks (plus and minus stranded):
	    - nucleotide-resolution ("exact cut") signal
	    - smoothed signal

## User interface

`PEPPRO` is a python script that runs on the command line (See [usage](usage)). It can also read projects in [PEP format](https://pepkit.github.io/). This means that `PEPPRO` projects are also compatible with other PEP tools, and output can be conveniently read into `R` using [the `pepr` package](http://code.databio.org/pepr/) or into `Python` using the [`peppy` package](https://peppy.readthedocs.io/en/latest/). The pipeline itself is customizable, enabling a user to adjust individual command settings or even swap out specific software by editing a few lines of human readable configuration files.

## Availability

You can download the latest version from the [releases](https://github.com/databio/peppro/releases) page, or find changes listed in the [CHANGELOG](changelog).

