# MpGAP pipeline changelog

The tracking for changes started in v2.1

## v2.1

This version have no additions to the pipeline workflow. It has additions in the modes of configuring and executing the pipeline, which are highlighted below.

### nf-core schema

We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility, also it is possible to track the pipeline execution with [nextflow tower](https://tower.nf/).

```bash
# It is triggered as
nf-core launch fmalmeida/mpgap
```

Checkout the paremeters `--use_tower` and `--tower_token` to activate pipeline execution in nextflow tower.
