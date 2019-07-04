# Running individual samples in a container

We have produced both docker and singularity containers that hold all the necessary software for `PEPPRO`. You can run `PEPPRO` as an individual pipeline on a single sample using these containers by directly calling `docker run` or `singularity exec`. Or, you can rely on `looper`, which is already set up to run any pipeline in existing containers using the `divvy` templating system.

Individual jobs can be run in a container by simply running the `peppro.py` command through `docker run` or `singularity exec`. You can run containers either on your local computer, or in an HPC environment, as long as you have `docker` or `singularity` installed. For example, run it locally in singularity like this:
```
singularity exec --bind $GENOMES $SIMAGES/peppro pipelines/peppro.py --help
```

With `docker`, you can use:
```
docker run --rm -it databio/peppro pipelines/peppro.py --help
```
Be sure to mount the volumes you need with `--volume`. If you're using any environment variables (e.g. `$GENOMES`), don't forget to include those in your docker command with the `-e` option.

### Detailed example using `docker`
The pipeline has been successfully run in both a Linux and MacOS environment. With `docker` you need to bind mount your volume that contains the pipeline and your `$GENOMES` location, as well as provide the container the same environment variables your host environment is using.

In the first example, we're mounting our home user directory (`/home/jps3ag/`) which contains the parent directories to our `$GENOMES` folder and to the pipeline itself. We'll also provide the pipeline two environment variables, `$GENOMES` and `$HOME`.

Here's that example command in a Linux environment to run the test example through the pipeline:
```
docker run --rm -it --volume /home/jps3ag/:/home/jps3ag/ \
  -e GENOMES='/home/jps3ag/genomes/' \
  -e HOME='/home/jps3ag/' \
  databio/peppro \
  /home/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignments rCRSd human_repeats \
  --genome hg38 \
  --sample-name test1 \
  --input /home/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  -O $HOME/peppro_test
```

In this second example, we'll perform the same command in a Mac environment using [Docker for Mac](https://docs.docker.com/v17.12/docker-for-mac/install/). 

This necessitates a few minor changes to run that same example:

- replace `/home/` with `/Users/` format
- e.g. `--volume /Users/jps3ag/:/Users/jps3ag/`

Remember to [allocate sufficient memory](https://docs.docker.com/docker-for-mac/#advanced) (6-8GB should generally be more than adequate) in Docker for Mac.

```
docker run --rm -it --volume /Users/jps3ag/:/Users/jps3ag/ \
  -e GENOMES="/Users/jps3ag/genomes" \
  -e HOME="/Users/jps3ag/" \
  databio/peppro \
  /Users/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignments rCRSd human_repeats \
  --genome hg38 \
  --sample-name test1 \
  --input /Users/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  -O $HOME/peppro_test
```

### Detailed example using `singularity`

First, build a singularity container from the docker image and create a running instance (be sure to mount your directories containing your `$GENOMES` folder and pipeline.
```
singularity build peppro docker://databio/peppro:latest
singularity instance.start -B /home/jps3ag/:/home/jps3aq/ peppro peppro_instance
```

Second, run your command.
```
singularity exec instance://peppro_instance \
  /home/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignments rCRSd human_repeats \
  --genome hg38 \
  --sample-name test1 \
  --input /home/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  -O $HOME/peppro_test
```

Third, close your instance when finished.
```
singularity instance.stop peppro_instance
```
