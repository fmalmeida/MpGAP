# MpGAP
General multi-platform genome assembly pipeline

MpGAP is a nextflow docker-based pipeline that wraps up [Canu](https://github.com/marbl/canu), [Flye](https://github.com/fenderglass/Flye), [Unicycler](https://github.com/rrwick/Unicycler), [Spades](https://github.com/ablab/spades), [Nanopolish](https://github.com/jts/nanopolish), [QUAST](https://github.com/ablab/quast) and [Pilon](https://github.com/broadinstitute/pilon).

This is an easy to use pipeline that uses well known software for genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes.

This pipeline has only two dependencies: [Docker](https://www.docker.com) and [Nextflow](https://github.com/nextflow-io/nextflow).
