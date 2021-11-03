.. _installation:

Installation
============

Dependencies
------------

This pipeline requires only `Docker <https://www.docker.com/>`_ (and its Docker images) and `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ to run.

1. Installing Docker

   + Read more in their `manual <https://docs.docker.com/>`_

2. Installing Nextflow

   ``curl -s https://get.nextflow.io | bash``

3. Download the pipeline

   ``./nextflow pull fmalmeida/mpgap``

4. Test your installation

   ``./nextflow run fmalmeida/mpgap --help``

5. Download required Docker images

   ``docker pull fmalmeida/mpgap:v3.0``

6. (Optional) Install nf-core

   ``pip install nf-core>=1.10``

.. note::

  Now, everything is set up and ready to run. Remember to always keep your Docker images up to date (Docker pull will always download the latest).

.. tip::

	The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the `Linux subsystem for windows <https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html>`_.
