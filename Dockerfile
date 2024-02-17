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
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/medaka && \
      chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.8/site-packages/medaka

# install ps
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

# return user
USER mambauser