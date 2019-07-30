# Running on a cluster

## Default computing options

When you run your PEPPRO project using `looper run`, by default it will simply run each sample locally. You can adjust the computing settings using the `--compute` argument to `looper run`. This enables you to adjust your computing preferences on-the-fly when you run a project. You have several built-in options to change this. Some common examples are these:

- `--compute slurm`. Submit the jobs to a SLURM cluster using `sbatch`.
- `--compute sge`. Submit the jobs to a SGE cluster using `qsub`.
- `--compute docker`. Submit the jobs locally using the `databio/peppro` docker image.
- `--compute singularity`. Submit the jobs locally using the singularity image.
- `--compute singularity_slurm`. Submit jobs using `sbatch`, but run them using the singularity image.

These available computing options are actually using a standardized computing system called [divvy](https://divvy.databio.org), and you can view a list of all your available options with this command:

```console
divvy list
```

To show how this works, let's run the example project using the `slurm` compute package. Used `-d` for a dry run to create the submits scripts but not run them:

```
cd peppro
looper run examples/meta/peppro_test.yaml -d \
  --compute slurm
```

This will give us a script produced, which we can look at:

```
cat peppro_test/submission/peppro_test.sub
```

If all looks well, run looper without `-d` to actually submit the jobs. To use the docker or singularity options, see [running PEPPRO in containers](container.md). 

## Customizing compute options

Divvy also allows you to very easily change these templates or add your own, so you can run PEPPRO in any possible computing environment. Changing these computing configuration options will work for any software that relies on `divvy`. For complete instructions, you should consult the [divvy documentation](https://divvy.databio.org). In a nutshell, you will first create a `divvy` computing configuration file and point an environment variable (`DIVCFG`) to that file. You then have access to any configured computing packages by using `looper --compute <package>`, where `package` can be any computing system you configure. The instructions for changing these computing configuration options are universal for any software that relies on `divvy`. 