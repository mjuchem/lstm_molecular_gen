FROM debian:unstable-slim

#####
# apt

ENV DEBIAN_PRIORITY=critical
ENV DEBIAN_FRONTEND=noninteractive

##############
# repositories

RUN truncate --size=0 /etc/apt/sources.list
COPY debian-testing.list /etc/apt/sources.list.d/
COPY debian-unstable.list /etc/apt/sources.list.d/

##########
# python 3

RUN apt-get update \
  && apt-get install -y \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-venv

RUN python3 -m venv /srv/venv

######################
# jupyter-lab and deps

RUN . /srv/venv/bin/activate \
  && pip install --upgrade \
    jupyterlab \
    matplotlib \
    numpy \
    pandas \
    rdkit \
    tensorflow

#############
# openbabel 3

RUN apt-get update \
  && apt-get install -y \
    libopenbabel-dev \
    swig

RUN . /srv/venv/bin/activate \
  && pip install --upgrade \
    --global-option="build_ext" \
    --global-option="-I/usr/include/openbabel3" \
    openbabel

############
# entrypoint

COPY entrypoint.sh /srv
ENTRYPOINT "/srv/entrypoint.sh"

#########
# cleanup

RUN apt-get upgrade -y
RUN apt autoremove --purge -y
RUN apt-get clean -y
RUN rm -rf /var/lib/apt/lists/*
