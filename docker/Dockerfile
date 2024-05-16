FROM debian:unstable-slim

#####
# apt

ENV DEBIAN_PRIORITY=critical
ENV DEBIAN_FRONTEND=noninteractive

##############
# repositories

RUN truncate --size=0 /etc/apt/sources.list
RUN rm -rf /etc/apt/sources.list.d
RUN mkdir -p /etc/apt/sources.list.d
COPY debian-testing.list /etc/apt/sources.list.d/
COPY debian-unstable.list /etc/apt/sources.list.d/

RUN apt-get update && apt-get upgrade -y

#############
# GPU drivers

RUN apt-get update \
  && apt-get install -y \
    build-essential dkms \
    vainfo \
      mesa-va-drivers \
    vdpauinfo \
      mesa-vdpau-drivers \
    vulkan-tools \
      mesa-vulkan-drivers \
    mesa-opencl-icd \
    intel-media-va-driver-non-free \
      intel-gpu-tools \
      intel-opencl-icd \
    nvidia-driver \
      nvidia-driver-full \
      nvidia-kernel-dkms \
      nvidia-vdpau-driver \
      nvidia-cuda-toolkit \
      nvidia-tesla-driver \
        nvidia-tesla-kernel-dkms \
        nvidia-tesla-opencl-icd \
        nvidia-tesla-vulkan-icd

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
  && CFLAGS="-march=native" \
    CPPFLAGS="-march=native" \
    pip install --upgrade \
      ipykernel \
      ipython \
      jupyterlab \
      matplotlib \
      numpy \
      pandas \
      rdkit \
      tensorflow \
      tensorrt

#############
# openbabel 3

RUN apt-get update \
  && apt-get install -y \
    libopenbabel-dev \
    swig

RUN . /srv/venv/bin/activate \
  && CFLAGS="-march=native" \
    CPPFLAGS="-march=native" \
    pip install --upgrade \
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
