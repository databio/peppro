<<<<<<< HEAD
=======
# Run <img src="../../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> in a container
>>>>>>> master

## Setting up PEPPRO containers

First, make sure your environment is set up to run either docker or singularity containers. Then, pull the container image (choose either docker or singularity). If your containers are set up correctly, then you won't need to install any additional software. 

### Docker

You can pull the docker [databio/peppro image](https://hub.docker.com/r/databio/peppro/) from dockerhub like this:

```
docker pull databio/peppro
```

Or build the image using the included Dockerfile (you can use a recipe in the included Makefile):
```
cd peppro/
make docker
```

### Singularity

You can [download the singularity image](http://big.databio.org/simages/peppro) or build it from the docker image using the Makefile:
```
cd peppro/
make singularity
```

Now you'll need to tell the pipeline where you saved the singularity image. You can either create an environment variable called `$SIMAGES` that points to the folder where your image is stored, or you can tweak the `pipeline_interface.yaml` file so that the `compute.singularity_image` attribute is pointing to the right location on disk.

<<<<<<< HEAD
 You can now use the containerized [compute packages](configure_compute.md), *e.g.*, `looper run --compute docker` or `looper run --compute singularity`.
=======
If your containers are set up correctly, then won't need to install any additional software. 

## Running individual samples in a container

Individual jobs can be run in a container by simply running the `peppro.py` command through `docker run` or `singularity exec`. You can run containers either on your local computer, or in an HPC environment, as long as you have `docker` or `singularity` installed. For example, run it locally in singularity like this:
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
The pipeline has been successfully run in both a Linux and MacOS environment. With `docker` you need to bind mount your volume that contains the pipeline and your `$REFGENIE` location, as well as provide the container the same environment variables your host environment is using.

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

## Running multiple samples in a container with looper

To run multiple samples in a container, you simply need to configure `looper` to use a container-compatible template. The looper documentation has detailed instructions for [how to run pipelines in containers](http://code.databio.org/looper/containers/).
>>>>>>> master
