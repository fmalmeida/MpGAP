FROM ubuntu:20.04

LABEL MAINTAINER "Felipe Marques de Almeida <200055801@aluno.unb.br>"
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
RUN export NCBI=$(curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ | grep "x64-linux.tar.gz$" | rev | cut -d' ' -f1 | rev) && \
 			wget "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/${NCBI}" -O blast.tar.gz
RUN 			tar -zxvf blast.tar.gz && \
			cp ncbi-blast-*/bin/* /usr/local/bin/

# INSTALL Racon
RUN git clone --recursive https://github.com/isovic/racon.git racon && \
			cd racon && mkdir build && cd build && \
			cmake -DCMAKE_BUILD_TYPE=Release .. && \
			make && make install

# INSTALL Medaka
RUN conda create -y -n MEDAKA "python=3.8" && \
			conda install -y -n MEDAKA -c conda-forge -c bioconda medaka nomkl && \
			conda clean -afy

# Set Pilon
RUN conda install -y -c bioconda pilon nomkl && \
			conda clean -afy

# Set Canu assembler
RUN curl -s https://api.github.com/repos/marbl/canu/releases/latest | grep "browser_download_url.*.Linux-amd64.tar.xz" | cut -d ":" -f 2,3  | tr -d '"' | wget -qi - -O canu.tar.xz && \
			tar -xJf canu.tar.xz && \
			ln -rs canu-*/*/bin/* /usr/local/bin

# Set SPAdes assembler
## install compatible spades for unicycler
## unicycler edge sorting functions of illumina graph
## only works until spades 3.13.0
RUN conda install -y -c bioconda "spades==3.13.0"	 nomkl && \
			conda clean -afy

# Set Unicycler
RUN pip install --upgrade pip && \
			pip3 install --upgrade pip && \
			hash -d pip && \
			hash -d pip3 && \
			pip install setuptools && \
			pip3 install setuptools && \
			apt install zlib1g-dev && \
			export UNICYCLER_TAG=$(curl -s https://api.github.com/repos/rrwick/Unicycler/releases/latest | grep "tag_name" | cut -d ":" -f 2,3  | tr -d '"' | tr -d "," | tr -d " ") && \
			wget "https://github.com/rrwick/Unicycler/archive/refs/tags/${UNICYCLER_TAG}.tar.gz" -O unicycler.tar.gz && \
			tar zxvf unicycler.tar.gz && \
			mv Unicycler-* Unicycler && \
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
			curl -s https://api.github.com/repos/lh3/minimap2/releases/latest | grep "browser_download_url.*x64-linux.tar.bz2" | cut -d ":" -f 2,3  | tr -d '"' | wget -qi - -O minimap.tar.bz2 && \
			tar -jxvf minimap.tar.bz2 && \
			ln -rs minimap2-*-linux/minimap2 /usr/local/bin && \
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

# Set PACBIO polisher (newest gcpp and mm2)
RUN conda create --name pacbio -y -c bioconda pbgcpp pbmm2 nomkl && \
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
RUN ln -frs canu-*/bin/* /usr/local/bin/

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
