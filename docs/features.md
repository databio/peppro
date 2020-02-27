# <img src="../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:50px; margin-top:-15px; margin-bottom:-10px"> features at-a-glance

Here are a few of the highlights that make `PEPPRO` valuable.

- **Scalability.** Run the pipeline easily on a project with a single sample or a thousand. This pipeline is compatible with [`looper`](https://looper.databio.org), so it can run locally, in a cloud container engine, or with any cluster resource manager (e.g. SLURM, SGE, or LFS).
- **Restartability.** The pipeline is built using [`pypiper`](https://pypiper.databio.org), so it automatically picks up where it left off in case of preemption or crash.
- **Copious logging.** The pipeline produces a detailed log file recording all output from every command run, and also records the time and memory use of every process, the version of the pipeline and other software, and other useful run information.
- **Flexibility.** The pipeline provides options for multiple peak callers, multiple adapter trimmers, and fully configurable parameterization for many underlying tools.
- **Portability.** Run it with [`bulker`](https://bulker.databio.org/en/latest/) which may use either `docker` or `singularity`, or it can be [run natively without containers](install.md). The choice is yours.
- **Standardized user interface.** The pipeline reads sample metadata formatted in [standard PEP format](http://pepkit.github.io/), so you can use the same sample annotation sheets for your downstream R or python analysis using tools from [pepkit](http://pepkit.github.io/).
- **Standardized reference genome assembly.** The pipeline uses standard reference genome assemblies produced by [refgenie](http://refgenie.databio.org), which provides a scripted way to produce a compatible reference assembly for any custom genome. For common genomes, you can either download pre-indexed assemblies or build your own.
- **Prealignments.** The pipeline can (optionally) first align to any number of reference assemblies separately before the primary genome alignment. This increases both speed and accuracy and can be used, for example, to align sequentially to rDNA, mtDNA, repeats, or spike-ins.
- **Nascent sequencing quality control metrics** The pipeline includes measures of nascent RNA purity, RNA integrity, library complexity, run-on efficiency and read feature distributions.
- **Beautiful `HTML` reports.** Your results include an easy-to-navigate `HTML` report with a sample table, job status, summary statistics, and QC plots at your fingertips when run using our [pipeline submission engine, looper](https://looper.databio.org).

