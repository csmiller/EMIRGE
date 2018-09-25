#
# VERSION               0.62.1

FROM      ubuntu:18.04
# For singularity on the hutch cluster
RUN export DEBIAN_FRONTEND=noninteractive
RUN mkdir /working
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
wget \
build-essential \
libtool \
automake \
zlib1g-dev \
libbz2-dev \
libncurses5-dev \
liblzma-dev \
unzip \
cython \
pkg-config
RUN mkdir -p /src/
WORKDIR /src/


# Pandaseq compile / install
RUN wget https://github.com/neufeld/pandaseq/archive/v2.11.tar.gz && \
tar xzvf v2.11.tar.gz
WORKDIR /src/pandaseq-2.11
RUN ./autogen.sh && ./configure && make && make install && ldconfig
WORKDIR /src
RUN rm -r /src/pandaseq* && rm -r /src/v2.11.tar.gz 

# VSEARCH
RUN wget https://github.com/torognes/vsearch/releases/download/v2.8.4/vsearch-2.8.4-linux-x86_64.tar.gz && \
tar xzvf vsearch-2.8.4-linux-x86_64.tar.gz && \
cp /src/vsearch-2.8.4-linux-x86_64/bin/vsearch /usr/local/bin/ && \
rm -r /src/vsearch*

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar xjvf samtools-1.9.tar.bz2
WORKDIR /src/samtools-1.9
RUN ./configure --prefix=/usr/local && make && make install
WORKDIR /src/
RUN rm -r /src/samtools*

# bowtie
WORKDIR /src/
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download -O bowtie.zip && \
unzip bowtie.zip && \
cp bowtie-1.2.2-linux-x86_64/bowtie* /usr/local/bin && \
rm -r /src/bowtie*

# USEARCH
#RUN wget https://golob.org/software/usearch10.0.240_i86linux32
#RUN cp /src/usearch10.0.240_i86linux32 /usr/local/bin/usearch
#RUN chmod +x /usr/local/bin/usearch

# EMIRGE (finally!)
RUN pip install numpy scipy pysam biopython cython
ADD . /src/emirge/
WORKDIR /src/emirge
RUN /usr/bin/python setup.py build
RUN /usr/bin/python setup.py install
RUN mkdir -p /emirge/db/
#WORKDIR /emirge/db/
#RUN emirge_makedb.py  --silva-license-accepted

# Cleanup
WORKDIR /
RUN rm -R /src