FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing all software requirements for the fmalmeida/mpgap pipeline"

# Install the conda environment
RUN conda install -y -c conda-forge mamba
COPY environment.yml /
RUN mamba env create --quiet -f /environment.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/mpgap-3.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name mpgap-3.1 > mpgap-3.1.yml

# download busco dbs
RUN quast-download-busco || true
ENV CONDA_PREFIX=/opt/conda
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.1/lib/python3.6/site-packages/quast_libs/busco/
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.1/lib/python3.6/site-packages/quast_libs/busco/bacteria.tar.gz https://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.1/lib/python3.6/site-packages/quast_libs/busco/eukaryota.tar.gz https://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.1/lib/python3.6/site-packages/quast_libs/busco/fungi.tar.gz https://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz

# fix permissions
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.1/lib/python3.6/site-packages/quast_libs/busco
RUN chmod -R 777 /opt/conda/envs/mpgap-3.1/lib/python3.6/site-packages/medaka

