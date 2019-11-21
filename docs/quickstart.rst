.. _quickstart:

Quickstart
**********

Overview
--------

We provide users a few test cases for evaluating the pipeline's commands and workflow.
We've made available two datasets:

* `Dataset 1 <https://drive.google.com/file/d/1xm3R97HXcfhsjSyhoTnvbA8HDK4Ij8ws/view?usp=sharing>`_.

    * Oxford Nanopore data (FAST5 and FASTQ);
    * Illumina paired end reads;

* `Dataset 2 <https://drive.google.com/file/d/1JmBPucItax7t2lH1XpqPxkdX9e7JLR6E/view?usp=sharing>`_

    * Pacbio data (subreads.bam);
    * Illumina paired end reads;

Getting the data
================

Users can directly download data through the link given or by this command line method below.

Add to your ``bashrc`` or ``bash_aliases``
""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

  function gdrive_download () {
  CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate \
  "https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
  wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
  rm -rf /tmp/cookies.txt
  }

.. tip::

  The function is used as: ``gdrive_download [gdrive id] [output name]``

Then, you can download the datasets as follows:

* **Dataset 1**

    * ``gdrive_download 1xm3R97HXcfhsjSyhoTnvbA8HDK4Ij8ws dataset_1.tar.gz``

* **Dataset 2**

    * ``gdrive_download 1JmBPucItax7t2lH1XpqPxkdX9e7JLR6E dataset_2.tar.gz``

Preprocessing the data
----------------------

Dataset 1
=========

After downloaded. the dataset shall be available as ``dataset_1`` directory. The first step, right after installing
the pipeline and downloading the docker image is to download the configuration file templates.

Download config files
"""""""""""""""""""""

.. code-block:: bash

  ## Get configuration for illumina data
  nextflow run fmalmeida/ngs-preprocess --get_illumina_config && mv illumina_data.config 01_illumina_data.config

  ## Get configuration for nanopore data
  nextflow run fmalmeida/ngs-preprocess --get_ont_config && mv ont_data.config 01_ont_data.config

After properly configuration of the files, they might look as this:

* `01_illumina_data.config <https://drive.google.com/file/d/1misoPDB66ai2J9cKhEyUKSO--H-933xv/view?usp=sharing>`_
* `01_ont_data.config <https://drive.google.com/file/d/16A3Uc6Ixqj-jYniSXPOSwNNzthKL3Ucz/view?usp=sharing>`_

Running the pipeline
""""""""""""""""""""

.. code-block:: bash

  ## Run for illumina
  nextflow run fmalmeida/ngs-preprocess -c 01_illumina_data.config &> 01_illumina_preprocess.log

  ## Run for nanopore
  nextflow run fmalmeida/ngs-preprocess -c 01_ont_data.config &> 01_ont_preprocess.log

Outputs will be at ``dataset_1/preprocessed``

Dataset 2
=========

After downloaded. the dataset shall be available as ``dataset_2`` directory. The first step, right after installing
the pipeline and downloading the docker image is to download the configuration file templates.

Download config files
"""""""""""""""""""""

.. code-block:: bash

  ## Get configuration for illumina data
  nextflow run fmalmeida/ngs-preprocess --get_illumina_config && mv illumina_data.config 02_illumina_data.config

  ## Get configuration for pacbio data
  nextflow run fmalmeida/ngs-preprocess --get_pacbio_config && mv pacbio_data.config 02_pacbio_data.config

After properly configuration of the files, they might look as this:

* `02_illumina_data.config <https://drive.google.com/file/d/17_lipuPHWOHUKj9TcW9ouDUpuzb7h3gQ/view?usp=sharing>`_
* `02_pacbio_data.config <https://drive.google.com/file/d/1gEsZ5KglbW-uYpYHnBrKIX7V5oTcMQuO/view?usp=sharing>`_

Running the pipeline
""""""""""""""""""""

.. code-block:: bash

  ## Run for illumina
  nextflow run fmalmeida/ngs-preprocess -c 02_illumina_data.config &> 02_illumina_preprocess.log

  ## Run for pacbio (we will use subreads.bam as input)
  nextflow run fmalmeida/ngs-preprocess -c 02_pacbio_data.config &> 02_pacbio_preprocess.log

Outputs will be at ``dataset_2/preprocessed``
