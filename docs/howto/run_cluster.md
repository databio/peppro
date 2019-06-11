# Run <img src="../../img/peppro_logo_black.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> on a cluster

`PEPPRO` by itself does not specify any cluster resources, so you could just roll your own and submit individual jobs to a cluster however you choose. But because `PEPPRO` is already `looper`-compatible, the easier way is to use `looper's` built-in template system, which `looper` uses to build flexible shell scripts for job submission. These templates can be used to run jobs in a container, to submit to a cluster resource manager, or both.

To use `looper` templates, we must create a `divvy` computing configuration file (compute_config.yaml) and point an environment variable (`DIVCFG`) to that file. You then have access to any configured computing packages by using `looper --compute <package>`, where `package` can be any computing system you configure.  In short, you will need to:

- Set up a compute configuration file that includes a containerized or cluster compute template (or both).
- Point the environment variable `DIVCFG` to the location of this file.
- Run the pipeline with `looper run --compute PACKAGE` (where `PACKAGE` is specified in your `DIVCFG` file).

This enables you to adjust your computing preferences on-the-fly when you run a project.

The complete description of setting up `looper` to use `DIVCFG` is generic to any pipeline. If you want to use looper with containers or clusters, you should consult the complete docs in the looper documentation on [configuring looper to use a cluster](http://code.databio.org/looper/cluster-computing/).
