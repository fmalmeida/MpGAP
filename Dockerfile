FROM mambaorg/micromamba
LABEL authors="Felipe Almeida" \
      description="Docker image containing all software requirements for the fmalmeida/mpgap pipeline"

# Install the conda environment
COPY environment.yml /
RUN micromamba env create --quiet -f /environment.yml --yes && micromamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/mpgap-3.3/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN micromamba env export --name mpgap-3.3 > mpgap-3.3.yml

# fix permissions
USER root
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.3/lib/python3.9/site-packages/medaka/data && \
      chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.3/lib/python3.9/site-packages/medaka

# install packages
RUN apt-get update && apt-get install -y procps zlib1g-dev build-essential && rm -rf /var/lib/apt/lists/*

# create dir
RUN mkdir -p /opt/busco_db/ && chmod -R 777 /opt/busco_db && chmod 777 /opt

# return user
USER mambauser

# download medaka pipeline's default model
RUN medaka tools download_models --model r941_min_high_g360

# pre-download BUSCO bacteria / fungi / eukaryota database
RUN wget -O /opt/busco_db/bacteria.tar.gz \
      https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz && \
      cd /opt/busco_db/ && \
      tar zxvf bacteria.tar.gz && \
      rm bacteria.tar.gz
RUN wget -O /opt/busco_db/fungi.tar.gz \
      https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2024-01-08.tar.gz && \
      cd /opt/busco_db/ && \
      tar zxvf fungi.tar.gz && \
      rm fungi.tar.gz
RUN wget -O /opt/busco_db/eukarya.tar.gz \
      https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz && \
      cd /opt/busco_db/ && \
      tar zxvf eukarya.tar.gz && \
      rm eukarya.tar.gz

USER root
RUN chmod -R 777 /opt/busco_db
USER mambauser

# build haslr tool due problems with conda version
# conda only used to install dependencies
RUN cd /opt && \
      git clone https://github.com/vpc-ccg/haslr.git && \
      cd haslr && \
      sed -i 's/-march=native/-msse4.1/g' src/haslr_assemble/lib/spoa.make && \
      /usr/bin/make
