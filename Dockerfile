FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing all software requirements for the fmalmeida/mpgap pipeline"

# Install the conda environment
RUN conda install -y -c conda-forge mamba
COPY environment.yml /
RUN mamba env create --quiet -f /environment.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/mpgap-3.2/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name mpgap-3.2 > mpgap-3.2.yml

# download busco dbs
ENV CONDA_PREFIX=/opt/conda
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/busco/
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/busco/bacteria.tar.gz https://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/busco/eukaryota.tar.gz https://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
RUN wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/busco/fungi.tar.gz https://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
RUN mkdir -p $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/augustus3.2.3 && wget -O $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/augustus3.2.3/augustus.tar.gz http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.3.tar.gz

# fix permissions
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/busco
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/quast_libs/augustus3.2.3
RUN chmod -R 777 $CONDA_PREFIX/envs/mpgap-3.2/lib/python3.6/site-packages/medaka

