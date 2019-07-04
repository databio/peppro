# FAQ

## Do I have to use PEPPRO with looper?

No. `PEPPRO` by itself does not specify any cluster resources, so you could just roll your own and submit individual jobs to a cluster however you choose. But because `PEPPRO` is already `looper`-compatible, the easier way is to use `looper's` built-in template system, which `looper` uses to build flexible shell scripts for job submission. These templates can be used to run jobs in a container, to submit to a cluster resource manager, or both.
