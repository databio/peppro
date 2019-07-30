# Running individual samples in a container

Individual jobs can be run in a container by simply running the `peppro.py` command through `docker run` or `singularity exec`. For example, run it locally in singularity like this:
```
singularity exec --bind $REFGENIE $SIMAGES/peppro pipelines/peppro.py --help
```

With `docker`, you can use:
```
docker run --rm -it databio/peppro pipelines/peppro.py --help
```
Be sure to mount the volumes you need with `--volume`. If you're utilizing any environment variables (e.g. `$REFGENIE`), don't forget to include those in your docker command with the `-e` option.

### Container details 

#### Using `docker`
With `docker` you need to bind mount your volume that contains the pipeline and your `$REFGENIE` location, as well as provide the container the same environment variables your host environment is using.

In the first example, we're mounting our home user directory (`/home/jps3ag/`) which contains the parent directories to our `$REFGENIE` folder and to the pipeline itself. We'll also provide the pipeline two environment variables, `$REFGENIE` and `$HOME`.

Here's that example command in a Linux environment to run the test example through the pipeline:
```
docker run --rm -it --volume /home/jps3ag/:/home/jps3ag/ \
  -e REFGENIE='/home/jps3ag/genomes/genome_config.yaml' \
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
  -e REFGENIE="/Users/jps3ag/genomes/genome_config.yaml" \
  -e HOME="/Users/jps3ag/" \
  databio/peppro \
  /Users/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignments rCRSd human_repeats \
  --genome hg38 \
  --sample-name test1 \
  --input /Users/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  -O $HOME/peppro_test
```

#### Using `singularity`

First, build a singularity container from the docker image and create a running instance (be sure to mount your directories containing your `$REFGENIE` configuration file, corresponding `refgenie` genomes, and the pipeline.
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




