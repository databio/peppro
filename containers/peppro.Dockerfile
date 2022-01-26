# Pull base image
FROM phusion/baseimage:master

# Who maintains this image
LABEL maintainer Jason Smith "jasonsmith@virginia.edu"

# Version info
LABEL version 0.10.0

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install dependencies
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes \    
    curl \
    cmake \
    default-jre \
    default-jdk \
    git \
    gsl-bin \
    libbz2-dev \
    libgsl-dbg \
    libgsl-dev \
    libcommons-math3-java \
    libcurl4-gnutls-dev \ 
    libjbzip2-java \
    libpng-dev \
    liblua5.1-0-dev \
    libisal-dev \
    libdeflate-dev \
    libssl-dev \
    libtbb2 \
    libtbb-dev \
    libbam-dev \
    libssl-dev \
    libtbb2 \
    libtbb-dev \
    lua-filesystem-dev \
    lua-lpeg-dev \
    lua-md5-dev \
    libexpat1-dev \
    libtre-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libsqlite3-dev \
    libxml2-dev \
    openssl \
    pigz=2.4-1 \
    python3.8 \
    python3-pip \
    python3-dev \
    software-properties-common \
    build-essential \
    rustc \
    wget \
    zlib1g \
    zlib1g-dev  

# Install MySQL server
RUN DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes mysql-server \
    mysql-client \
    libmysqlclient-dev
    
# Install python tools
RUN python3.8 -m pip install -U pip
RUN pip install attmap>=0.12.9 && \
    pip install cython>=0.29 && \
    pip install cykhash && \
    pip install jinja2>=2.11.2 && \
    pip install jsonschema>=3.0.1 && \
    pip install logmuse>=0.2.5 && \
    pip install numpy>=1.17 && \
    pip install https://github.com/pepkit/looper/zipball/master && \
    pip install pararead && \
    pip install pandas>=0.20.2 && \
    pip install peppy>=0.31.0 && \
    pip install piper>=0.12.1 && \
    pip install psutil>=5.6.3 && \
    pip install pysam>=0.13 && \
    pip install pyyaml>=3.13 && \
    pip install refgenconf>=0.7.0 && \
    pip install refgenie>=0.9.3 && \
    pip install ubiquerg>=0.6.1 && \
    pip install yacman>=0.6.7 && \
    pip install cutadapt

# Install R
RUN apt update -qq && \
    DEBIAN_FRONTEND=noninteractive apt --assume-yes install --no-install-recommends dirmngr
RUN apt-key adv --keyserver  hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN DEBIAN_FRONTEND=noninteractive apt-get --assume-yes install r-base r-base-dev r-base-core r-recommended && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN Rscript -e "install.packages('argparser')" && \
    Rscript -e "install.packages('optigrab')" && \
    Rscript -e "install.packages('data.table')" && \
    Rscript -e "install.packages('xml2')" && \
    Rscript -e "install.packages('roxygen2')" && \
    Rscript -e "install.packages('rversions')" && \
    Rscript -e "install.packages('callr')" && \
    Rscript -e "install.packages('pkgbuild')" && \
    Rscript -e "install.packages('rcmdcheck')" && \
    Rscript -e "install.packages('testthat')" && \
    Rscript -e "install.packages('devtools')"
    
RUN Rscript -e "devtools::install_github('pepkit/pepr')" && \    
    Rscript -e "install.packages('data.table')" && \
    Rscript -e "install.packages('BiocManager')" && \
    Rscript -e "BiocManager::install('GenomicRanges')" && \
    Rscript -e "BiocManager::install('BSgenome')" && \
    Rscript -e "BiocManager::install('GenomicFeatures')" && \
    Rscript -e "BiocManager::install('ensembldb')" && \
    Rscript -e "BiocManager::install('ExperimentHub')" && \
    Rscript -e "devtools::install_github('databio/GenomicDistributions')" && \
    Rscript -e "install.packages('http://big.databio.org/GenomicDistributionsData/GenomicDistributionsData_0.0.2.tar.gz', repos=NULL)" &&\
    Rscript -e "install.packages('ggrepel')" && \
    Rscript -e "install.packages('ggplot2')" && \  
    Rscript -e "install.packages('gplots')" && \
    Rscript -e "install.packages('grid')" && \    
    Rscript -e "install.packages('gtable')" && \
    Rscript -e "install.packages('scales')" && \
    Rscript -e "install.packages('stringr')" && \
    Rscript -e "devtools::install_github('databio/peppro/PEPPROr/', ref = 'master')"

# Install htslib
WORKDIR /home/src/
RUN wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
    tar xf htslib-1.12.tar.bz2 && \
    cd /home/src/htslib-1.12 && \
    ./configure --prefix /home/src/ && \
    make && \
    make install

# Install samtools
WORKDIR /home/src/
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar xf samtools-1.12.tar.bz2 && \
    cd /home/src/samtools-1.12 && \
    ./configure && \
    make && \
    make install && \
    ln -s /home/src/samtools-1.12/samtools /usr/bin/

# Install bedtools
RUN DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes \
    ant \
    bedtools>=2.29.2 

# Install bowtie2
WORKDIR /home/src/
RUN wget -O bowtie2-2.4.2-source.zip 'https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.4.2/bowtie2-2.4.2-source.zip?ts=gAAAAABgfxZxKMUjBU0A0XjfO55q36KUoO9RRemjzTT_WCDpSSZCy8NtKrFODKV4xS_135KTiIdnBSaqmvHuQw9l6nqM2EULvw%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.4.2%2Fbowtie2-2.4.2-source.zip%2Fdownload' && \
    unzip bowtie2-2.4.2-source.zip && \
    cd /home/src/bowtie2-2.4.2 && \
    make && \
    make install && \
    ln -s /home/src/bowtie2-2.4.2/bowtie2 /usr/bin/

# Install seqkit
WORKDIR /home/src/
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.10.1/seqkit_linux_amd64.tar.gz && \
    tar -zxvf seqkit_linux_amd64.tar.gz && \
    ln -s /home/src/seqkit /usr/bin/

# Install fastp
WORKDIR /home/src/
RUN git clone https://github.com/OpenGene/fastp.git && \
    cd fastp && \
    make && \
    make install && \
    ln -s /home/src/fastp/fastp /usr/bin/

# Install seqtk
WORKDIR /home/src/
RUN git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    make && \
    ln -s /home/src/seqtk/seqtk /usr/bin/

# Install preseq
WORKDIR /home/src/
RUN wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz && \
    tar xf preseq-3.1.2.tar.gz && \
    cd preseq-3.1.2 && \
    mkdir build && cd build && \
    ../configure --enable-hts \
        CPPFLAGS=-I"/home/src/include" \
        LDFLAGS="-L/home/src/lib -Wl,-R/home/src/lib" && \
    make && \
    make install

# Install UCSC tools
WORKDIR /home/tools/
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat && \
    chmod +x /home/tools/wigToBigWig && \
    chmod +x /home/tools/bigWigCat && \
    ln -s /home/tools/wigToBigWig /usr/bin/ && \
    ln -s /home/tools/bigWigCat /usr/bin/

# Install FLASH
WORKDIR /home/src/
RUN wget -O FLASH-1.2.11.tar.gz http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-Linux-x86_64.tar.gz && \
    tar xf FLASH-1.2.11.tar.gz && \
    ln -s /home/src/FLASH-1.2.11-Linux-x86_64/flash /usr/bin/

# Install fastq_pair
WORKDIR /home/src/
RUN git clone https://github.com/linsalrob/fastq-pair.git && \
    cd fastq-pair/ &&\
    mkdir build && cd build && \
    cmake /home/src/fastq-pair/ && \
    make && \
    make install

# OPTIONAL REQUIREMENTS
# Install fastqc
WORKDIR /home/tools/
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    cd /home/tools/FastQC && \
    chmod 755 fastqc && \ 
    ln -s /home/tools/FastQC/fastqc /usr/bin/

# Install fqdedup
WORKDIR /home/tools/
RUN git clone https://github.com/guertinlab/fqdedup.git && \
    cd fqdedup && \
    cargo build --release && \
    ln -s /home/tools/fqdedup/target/release/fqdedup /usr/bin/

# Install fastx_toolkit
WORKDIR /home/tools/
RUN git clone https://github.com/agordon/libgtextutils.git && \
    cd libgtextutils && \
    ./reconf && \
    ./configure && \
    make && \
    make install && \
    cd /home/tools/ && \
    git clone https://github.com/agordon/fastx_toolkit && \
    cd fastx_toolkit && \
    ./reconf && \
    sed -i 's/-Werror//g' configure.ac && \
    ./configure && \
    make && \
    make install

# Install genometools
WORKDIR /home/tools/
RUN wget http://genometools.org/pub/genometools-1.6.1.tar.gz && \
    tar xf genometools-1.6.1.tar.gz && \
    cd /home/tools/genometools-1.6.1 && \
    make useshared=yes && \
    make install

# Install seqOutBias
WORKDIR /home/tools/
RUN wget -O seqOutBias-v1.3.0.tar.gz 'https://github.com/guertinlab/seqOutBias/archive/refs/tags/v1.3.0.tar.gz' && \
    tar xf seqOutBias-v1.3.0.tar.gz && \
    cd seqOutBias-1.3.0 && \
    cargo build --release && \
    ln -s /home/tools/seqOutBias-1.3.0/target/release/seqOutBias /usr/bin/

# Set environment variables
ENV PATH=/home/tools/bin:/home/tools/:/home/tools/bin/kentUtils/:/home/src/bowtie2-2.4.2:/home/src/skewer:/home/src/samtools-1.12:/home/src/htslib-1.12:$PATH \
    R_LIBS_USER=/usr/local/lib/R/site-library/ \
    PYTHONPATH=/usr/local/lib/python3.8/dist-packages:$PYTHONPATH

# Set Python3 as python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# Define default command
WORKDIR /home/
CMD ["/bin/bash"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

