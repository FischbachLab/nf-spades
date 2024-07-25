# Base Image
FROM continuumio/miniconda
#FROM jcrist/alpine-conda
#:4.6.8
LABEL authors="Xiandong Meng" \
      description="Docker image for Hybrid assembly for VITAL project" \
      software.name="metaspades" \
      software.description="SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines." \
      software.website="http://cab.spbu.ru" \
      container.category="tool"

# Update conda to latest version.
RUN /opt/conda/bin/conda update -n base -c defaults conda

# Install software packages
RUN /opt/conda/bin/conda install -c conda-forge -c bioconda -y spades=3.15.4 samtools bwa bbmap fastqc qualimap
RUN /opt/conda/bin/conda install -c conda-forge -y awscli

# Install time and memory usage tool (timem)
RUN /opt/conda/bin/conda install -c jrmadsen -y timemory

RUN mkdir -p /work
WORKDIR /work

# Get Repo
COPY . .

RUN chmod -R +rx ./
