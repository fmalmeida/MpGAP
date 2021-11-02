.. _config:

Configuration File
==================

.. code-block:: bash

  # To download a configuration file template users just need to run:
  nextflow run fmalmeida/mpgap [--get_config]

  # Using a config file your code is a lot more clean and concise:
  nextflow run fmalmeida/mpgap -c [path-to-config]

Main config file
----------------

.. include:: ../nextflow.config
   :code: groovy