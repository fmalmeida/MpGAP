.. _installation:

Installation
============

Dependencies
------------

This pipeline requires only `Docker <https://www.docker.com/>`_ or `Singularity <https://sylabs.io/docs/>`_ and `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ to run.

1. Installing Docker or Singularity

   + Read more in the `docker <https://docs.docker.com/>`_ and `singularity <https://sylabs.io/docs/>`_ manuals

2. Installing Nextflow

   ``curl -s https://get.nextflow.io | bash``

3. Download the pipeline

   ``./nextflow pull fmalmeida/mpgap``

4. Test your installation

   ``./nextflow run fmalmeida/mpgap --help``

5. Download required Docker images

   .. code-block:: bash

      # for docker
      docker pull fmalmeida/mpgap:v3.1

      # for singularity
      # remember to properly set NXF_SINGULARITY_LIBRARYDIR
      # read more at https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub
      export NXF_SINGULARITY_LIBRARYDIR=MY_SINGULARITY_IMAGES
      singularity pull --dir MY_SINGULARITY_IMAGES docker://fmalmeida/mpgap:v3.1
   
   * To understand how to change between profiles please refer to https://github.com/fmalmeida/MpGAP/tree/master#selecting-between-profiles

      * You can also run it with conda

6. (Optional) Install nf-core

   ``pip install nf-core>=1.10``

.. tip::

	The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the `Linux subsystem for windows <https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html>`_.
