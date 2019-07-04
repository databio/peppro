
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

 You can now use the containerized [compute packages](configure_compute.md), *e.g.*, `looper run --compute docker` or `looper run --compute singularity`.