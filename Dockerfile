#
# VERSION               golob/emirge:0.62.1D

FROM      --platform=amd64 ubuntu:18.04


RUN export DEBIAN_FRONTEND=noninteractive
RUN mkdir /working
RUN apt-get update && apt-get install -y \
python-dev \
python-pip \
wget \
build-essential \
libtool \
automake \
zlib1g-dev \
libbz2-dev \
libncurses5-dev \
liblzma-dev \
unzip \
pkg-config \
python-pip && apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
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

# EMIRGE (finally!)
RUN pip install \
numpy \
scipy \
pysam \
biopython==1.76 \
setuptools \
cython

RUN mkdir -p /src/emirge/ && mkdir -p /src/emirge/utils && mkdir -p /src/emirge/pykseq
ADD emirge_makedb.py /usr/local/bin/
ADD ./*.c /src/emirge/
ADD ./*.py /src/emirge/
ADD ./*.pyx /src/emirge/
ADD ./*.h /src/emirge/
ADD ./utils/* /src/emirge/utils/
ADD ./pykseq /src/emirge/pykseq/
WORKDIR /src/emirge
RUN python setup.py build
RUN python setup.py install

# Finally the DB
RUN mkdir -p /emirge/db/arf_c100/
WORKDIR /emirge/db/arf_c100/
ADD SSU_candidates.fasta.gz /emirge/db/arf_c100/
RUN  gunzip -c SSU_candidates.fasta.gz > SSU_candidates.fasta && \
bowtie-build SSU_candidates.fasta SSU_candidates_btindex --threads 4 && \
rm SSU_candidates.fasta
RUN gunzip SU_candidates.fasta.gz
WORKDIR /emirge/

# Cleanup
WORKDIR /
RUN rm -R /src