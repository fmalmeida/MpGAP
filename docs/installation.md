# Installation

## Dependencies

The pipeline require only a UNIX system, [Nextflow](https://www.nextflow.io/docs/latest/index.html#) and either [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/docs/) or [conda](https://conda.io/). Please, for installing these tools refer to their manual.

!!! note "About NF profiles"

    Please read more about how to [proper select NF profiles](profiles.md#) to better understand it.

## Downloading the pipeline

You can easily get a copy of the pipeline with:

```bash
# nextflow pull
nextflow pull fmalmeida/mpgap
```

!!! warning
    
    The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the [Linux subsystem for window](https://docs.microsoft.com/pt-br/windows/wsl/install-win10). Nextflow team has made available a [nice tutorial](https://www.nextflow.io/blog.html) about this issue.

## Downloading docker images

> All images can be downloaded on the fly, automatically by nextflow, and this is the recommended way to do it.

!!! info "If using singularity"

    **Docker and singularity images are downloaded on the fly**. Be sure to properly set `NXF_SINGULARITY_LIBRARYDIR` env variable to a writable directory if using Singularity. This will make that the downloaded images are reusable through different executions. Read more at: https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub

    For example, you would:

    ```bash
    export NXF_SINGULARITY_LIBRARYDIR=MY_SINGULARITY_IMAGES    # Set a path to your singularity storage dir
    export NXF_SINGULARITY_CACHEDIR=MY_SINGULARITY_CACHE       # Set a path to your singularity cache dir

    # run
    nextflow run fmalmeida/mpgap -profile singularity [options]
    ```

!!! info "If using conda"

    You would need to first download the environment file and create the pipeline's conda environment. Example:

    ```bash
    wget https://github.com/fmalmeida/mpgap/raw/master/environment.yml
    conda env create -f environment.yml # advice: use mamba
    ```

Now, everything is set up and ready to run. Remember to always keep your Docker images up to date (Docker pull will always download the latest).
