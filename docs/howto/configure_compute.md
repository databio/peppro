# Configuring <img src="../../img/peppro_logo_black.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> computing settings

## Default computing options

When you run your PEPPRO project using `looper run`, by default it will simply run each sample locally. You can adjust the computing settings using the `--compute` argument to `looper run`. This enables you to adjust your computing preferences on-the-fly when you run a project. You have several built-in options to change this. Some common examples are these:

- `--compute slurm`. Submit the jobs to a SLURM cluster using `sbatch`.
- `--compute sge`. Submit the jobs to a SGE cluster using `sbatch`.
- `--compute docker`. Submit the jobs locally using the `databio/peppro` docker image.
- `--compute singularity`. Submit the jobs locally using the `peppro` singularity image.
- `--compute singularity_slurm`. Submit the jobs using `sbatch`, but then run them using the `peppro` singularity image.

To use the docker or singularity options, you'll need to make sure you're [set up for using PEPPRO containers](use_container.md). These available computing options are actually using a standardized computing system called [divvy](https://divvy.databio.org), and you can view a list of all your available options with this command:

```console
divvy list
```

## Customizing compute options

Divvy also allows you to very easily change these templates or add your own, so you can run PEPPRO in any possible computing environment. Changing these computing configuration options will work for any software that relies on `divvy`. For complete instructions, you should consult the [divvy documentation](https://divvy.databio.org). In a nutshell, you will first create a `divvy` computing configuration file (compute_config.yaml) and point an environment variable (`DIVCFG`) to that file. You then have access to any configured computing packages by using `looper --compute <package>`, where `package` can be any computing system you configure.  


