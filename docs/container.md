# Running in a container

We have produced both docker and singularity containers that hold all the necessary software for `PEPPRO`. If your containers are set up correctly, then you won't need to install any additional software. It's easy to run your jobs in a container by configuring `looper` to use a container-compatible template. Follow the instructions below for either `docker` or `singularity` as you wish:

## Run PEPPRO using docker

You can pull the docker [databio/peppro image](https://hub.docker.com/r/databio/peppro/) from dockerhub like this:

```
docker pull databio/peppro
```

Or build the image using the included Dockerfile (you can use a recipe in the included Makefile):
```
cd peppro/
make docker
```

Next, just add `--compute docker` to your `looper run` command:

```
cd peppro
looper run examples/meta/peppro_test.yaml --compute docker
```

## Run PEPPRO using singularity

You can [download the singularity image](http://big.databio.org/simages/peppro) or build it from the docker image using the Makefile:
```
cd peppro/
make singularity
```

Now you'll need to tell the pipeline where you saved the singularity image. By default PEPPRO expects you to put your singularity image in a folder referred to with an environment variable called `$SIMAGES`:

```
export SIMAGES=path/to/singularity/folder/
```

You could also tweak the `pipeline_interface.yaml` file so that the `compute.singularity_image` attribute is pointing to the right location on disk. Run it like this:

```
cd peppro
looper run examples/meta/peppro_test.yaml \
  --compute singularity
```

## Run PEPPRO using singularity with SLURM

```
cd peppro
looper run examples/meta/peppro_test.yaml \
  --compute singularity_slurm
```

## More details on containers

You may need to adjust the built-in templates to fit with how you run docker or singularity in your environment. For example, you'll need to make sure to mount any filesystems you need. To do this pretty easy, you will just need to tweak the default templates to fit your environment. Here are some resources to get you started:

- [Divvy documentation on container templates](http://divvy.databio.org/en/latest/containers/)
- Looper documentation detailed instructions for [how to run pipelines in containers](http://looper.databio.org/en/latest/containers/).
