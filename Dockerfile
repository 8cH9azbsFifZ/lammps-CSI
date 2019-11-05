FROM debian:stretch-slim
MAINTAINER Gerolf Ziegenhain "gerolf.ziegenhain@gmail.com"

RUN apt-get update && \
    apt-get install -y wget build-essential ssh zlib1g-dev\
                       libfftw3-dev libopenblas-dev libopenmpi-dev && \
    rm -rf /var/lib/apt/lists/*

ADD . /lammps
RUN cd /lammps/STUBS ; make ; cd .. ; make serial
