# MpGAP
General multi-platform genome assembly pipeline

MpGAP is a nextflow docker-based pipeline that wraps up [Canu](https://github.com/marbl/canu), [Flye](https://github.com/fenderglass/Flye), [Unicycler](https://github.com/rrwick/Unicycler), [Spades](https://github.com/ablab/spades), [Nanopolish](https://github.com/jts/nanopolish), [QUAST](https://github.com/ablab/quast) and [Pilon](https://github.com/broadinstitute/pilon).

This is an easy to use pipeline that uses well known software for genome assembly of Illumina, Pacbio and Oxford Nanopore sequencing data through illumina only, long reads only or hybrid modes.

This pipeline has only two dependencies: Docker and Nextflow.

This pipeline requires only the installation of Docker and Nextflow. It is able to execute long reads only, illumina only or hybrid assemblies using these tools cited above.
