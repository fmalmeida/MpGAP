.. _installation:

Installation
============

Dependencies
------------

This pipeline requires only `Docker <https://www.docker.com/>`_ (and its Docker images) and
`Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ to run.

1. Installing Docker
  * Read more in their `manual <https://docs.docker.com/>`_
  * Or give this `in-house script <https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh>`_ a try.
2. Installing Nextflow

    ``curl -s https://get.nextflow.io | bash``

3. Download the pipeline

    ``./nextflow pull fmalmeida/ngs-preprocess``

4. Test your installation

    ``./nextflow run fmalmeida/ngs-preprocess --help``

5. Download required Docker images

    ``docker pull fmalmeida/mpgap``

.. note::

  Now, everything is set up and ready to run.
  Remember to always keep your Docker images up to date
  (Docker pull will always download the latest).
