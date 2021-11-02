# Running earlier releases

## How?

By default, nextflow will execute the code that is available in master/main branch of a repository. This is the recommended action because the main branch will always hold the most up-to-date code with bug fixes from other versions.

However, if for some reason you need to execute a different version of a pipeline, you can use the parameter `-r branch/tag`. For example:

```bash
# to run the code in tag/branch v2.0
nextflow run fmalmeida/mpgap -r v2.0 --help
```

One can check all available versions with the command:

```bash
# get information on pipeline
nextflow info fmalmeida/mpgap
```

## However

Please note that even though I try to maintain all the docker images available in dockerhub, docker automatically deactivates and deletes images that are not used for a certain amount of time.

If you try to execute a earlier release, and the image is not available anymore, you will have to access the Dockerfile of that specific image and built it locally. However, if a bug ar any problem arise when using an earlier release, I cannot provide much assistance and maintainance.

That's why I recommend always using the latest version available.
