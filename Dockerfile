FROM mambaorg/micromamba
LABEL authors="Felipe Almeida" \
      description="Docker image containing all software requirements for the fmalmeida/mpgap pipeline"

# Install the conda environment
COPY environment.yml /
RUN micromamba env create --quiet -f /environment.yml --yes && micromamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/mpgap-3.2/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN micromamba env export --name mpgap-3.2 > mpgap-3.2.yml

# check problematic installation
RUN medaka --help

# download busco dbs
ENV CONDA_PREFIX=/opt/conda
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/busco/
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/busco/bacteria.tar.gz https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/busco/eukaryota.tar.gz https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/busco/fungi.tar.gz https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2024-01-08.tar.gz
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/augustus3.2.3 && wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/augustus3.2.3/augustus.tar.gz http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.3.tar.gz

# fix permissions
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/busco
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/quast_libs/augustus3.2.3
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/medaka

# install ps
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
USER mambauser