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

# fix permissions
USER root
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.9/site-packages/medaka && \
      chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.9/site-packages/medaka

# pre-download BUSCO bacteria database
RUN mkdir -p /opt/busco_db/ && \
      wget -O /opt/busco_db/bacteria.tar.gz https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz && \
      cd /opt/busco_db/ && \
      tar zxvf bacteria.tar.gz && \
      rm bacteria.tar.gz && \
      chmod -R 777 /opt/busco_db/

# install ps
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

# return user
USER mambauser