# Use a LTS Ubuntu version as parent image
FROM ubuntu:20.04
# NOTE: without interactive dialogue needed to install r-base
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# metadata
LABEL base.image="ubuntu:20.04"
LABEL container.version="0.1"
LABEL software="assembly_pipeline"
LABEL software.version="1.0"
LABEL description="assembly_pipeline.py is a computational pipeline to produce polished bacterial de novo assemblies from paired Illumina data using Spaces and improve_assembly pipeline"
LABEL website="https://github.com/francesccoll/assembly_pipeline"
LABEL lisence="https://github.com/francesccoll/assembly_pipeline/blob/main/LICENSE"
LABEL maintainer="Francesc Coll"
LABEL maintainer.email="francesc.coll@lshtm.ac.uk"


# Install general dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    curl \
    build-essential \
    automake \
    pkg-config \
    zlib1g-dev \
    unzip \
    autoconf \
    check \
    libtool \
    libsubunit-dev \
    git \
    wget \
    python3.8 \
    python3.8-dev \
    python3-setuptools \
    python3-pip


# Installing fastqcheck
# Docker commands extracted from https://hub.docker.com/r/sangerpathogens/fastqcheck/dockerfile 
RUN git clone https://github.com/sanger-pathogens/fastqcheck.git
RUN cd fastqcheck && autoreconf -i -f && ./configure && make && make install
ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH


# Downloading SPAdes Linux binaries
RUN wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz
RUN tar -xzf SPAdes-3.15.3-Linux.tar.gz
ENV PATH="${PATH}:/opt/SPAdes-3.15.3-Linux/bin"
# NOTE: 'ln -s /usr/bin/python3 /usr/bin/python' used so that spades can run with header '#!/usr/bin/env python' using python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# assembly improvement required dependancies
# See https://github.com/sanger-pathogens/assembly_improvement
RUN apt-get -y install cpanminus
RUN cpanm -f Bio::AssemblyImprovement

# Installing quast
# docker commands extracted from: https://github.com/ablab/quast/blob/master/Dockerfile
RUN apt-get update && apt-get install -y pkg-config libfreetype6-dev \
    libpng-dev python3-matplotlib

# COPY . quast
RUN git clone https://github.com/ablab/quast.git
# RUN pip3 install --upgrade setuptools pip && \
RUN cd quast && \
    ./setup.py install
    # python setup.py develop && \
    # cd ..


# And finally, make assembly_pipeline.py script available and executable on image

WORKDIR /usr/local/bin
COPY assembly_pipeline.py .
RUN chmod +x assembly_pipeline.py