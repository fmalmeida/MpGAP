FROM ubuntu:20.04

MAINTAINER Felipe Marques de Almeida <200055801@aluno.unb.br>
SHELL ["/bin/bash", "-c"]

# SET DIR
WORKDIR /work

# Dependencies
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && \
			apt-get install -y git libc6-dev build-essential gcc g++ make dpkg-dev python-pip-whl samtools libidn11 cmake wget curl python3 python3-pip libpbcopper-dev libpbcopper1.3.0

# Set miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
			bash /miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# Set ncbi-blast+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz && \
        tar -zxvf ncbi-blast-2.11.0+-x64-linux.tar.gz && \
				cp ncbi-blast-2.11.0+/bin/* /usr/local/bin/

# INSTALL Racon
RUN git clone --recursive https://github.com/isovic/racon.git racon && \
			cd racon && mkdir build && cd build && \
			cmake -DCMAKE_BUILD_TYPE=Release .. && \
			make && make install

# INSTALL Medaka
RUN conda create -y -n MEDAKA python=3.8 && \
			conda install -y -n MEDAKA -c conda-forge -c bioconda medaka nomkl && \
			conda clean -afy

# Set Pilon
RUN conda install -y -c bioconda pilon nomkl && \
			conda clean -afy

# Set Canu assembler
RUN wget https://github.com/marbl/canu/releases/download/v2.1.1/canu-2.1.1.Linux-amd64.tar.xz && \
			tar -xJf canu-2.*.tar.xz && \
			ln -rs canu-2.*/*/bin/* /usr/local/bin

# Set SPAdes assembler
RUN conda install -y -c bioconda spades nomkl && \
			conda clean -afy

# Set Unicycler
RUN pip install --upgrade pip && \
			pip3 install --upgrade pip && \
			hash -d pip && \
			hash -d pip3 && \
			pip install setuptools && \
			pip3 install setuptools && \
			apt install zlib1g-dev && \
			wget https://github.com/rrwick/Unicycler/archive/v0.4.8.tar.gz && \
			tar zxvf v0.4.8.tar.gz && \
			mv Unicycler-0.4.8 Unicycler && \
			cd Unicycler && python3 setup.py install

## Set Unicycler complementary tools
RUN apt-get install -y mummer && \
			#conda install -y -c bioconda samtools && \
			conda install -y -c bioconda bowtie2 nomkl && \
			conda install -y -c bioconda freebayes nomkl && \
			#conda create -y -n legacy_python python=2.7 && conda install -y -c bioconda -n legacy_python pbalign && \
			#conda install -y -n legacy_python -c bioconda bax2bam && \
			#conda install -y -c conda-forge arrow && \
			git clone https://github.com/sc932/ALE.git && cd ALE/src && make && ln -rs ALE /usr/local/bin && \
			conda install -y -c bioconda minimap nomkl && \
			curl -L https://github.com/lh3/minimap2/releases/download/v2.10/minimap2-2.10_x64-linux.tar.bz2 | tar -jxvf - ; ln -rs minimap2-2.10_x64-linux/minimap2 /usr/local/bin && \
			chmod a+x /miniconda/bin/* && \
			conda clean -afy

# Set Nanopolish
RUN apt-get install -y libhdf5-dev parallel && \
    	conda config --add channels defaults && \
      conda config --add channels conda-forge && \
      conda config --add channels bioconda && \
      conda create -y -n NANOPOLISH && \
      conda install -y -n NANOPOLISH nanopolish nomkl && \
			conda clean -afy

# Set PACBIO GENOMICONSENSUS
RUN conda create --name pacbio python=2 && \
			conda install -y -c bioconda -n pacbio genomicconsensus openssl=1.0 pbalign blasr bax2bam pbcore pbcommand pbcopper nomkl && \
			chmod a+x /miniconda/bin/* && \
			conda clean -afy

# Set Flye
RUN conda create -y -n flye && \
			conda install -y -n flye -c bioconda flye nomkl && \
			conda clean -afy

# Set Quast
RUN python -m pip install -U matplotlib && \
			apt-get install -y pkg-config libfreetype6-dev libpng-dev bwa && \
			git clone https://github.com/ablab/quast.git && \
			cd quast && pip install --upgrade setuptools pip && \
			./setup.py install_full

# Set miniasm
RUN conda install -y -c bioconda miniasm nomkl && \
			conda clean -afy

# Custom Permissions
RUN chmod -R 777 /miniconda/envs/MEDAKA/lib/python3.8/site-packages/medaka/data

# Check envs
RUN conda info -e

# Fix Canu links
RUN ln -frs canu-2.1.1/bin/* /usr/local/bin/

# Install seqtk
RUN git clone https://github.com/lh3/seqtk.git && \
				cd seqtk && \
				make && \
				ln -rs seqtk /usr/local/bin

# Install shovill
RUN conda create --name shovill -y -c conda-forge -c bioconda -c defaults shovill nomkl && \
			conda clean -afy

# Install haslr
RUN git clone https://github.com/vpc-ccg/haslr.git && \
			cd haslr && \
			make
ENV PATH="$PATH:/work/haslr/bin"

# Install raven
RUN conda create -y -n RAVEN && \
			conda install -y -n RAVEN -c bioconda raven-assembler nomkl

# Install multiqc
RUN conda install -y -c bioconda -c conda-forge multiqc nomkl && \
			conda clean -afy

# fix tbb libraries
RUN conda install -y -c conda-forge 'tbb=2020.2' && \
			conda clean -afy

# Remove tars
RUN rm -rf *.tar.gz *.tar.xz

# Set WorkDir
WORKDIR /work
