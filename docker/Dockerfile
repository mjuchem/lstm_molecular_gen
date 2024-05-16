FROM debian:testing-slim

ENV DEBIAN_PRIORITY=critical
ENV DEBIAN_FRONTEND=noninteractive

RUN truncate --size=0 /etc/apt/sources.list
COPY debian-testing.list /etc/apt/sources.list.d/
COPY debian-unstable.list /etc/apt/sources.list.d/

RUN apt-get update \
  && apt-get install -y \
    python3 \
    python3-matplotlib \
    python3-numpy \
    python3-openbabel \
    python3-pandas \
    python3-rdkit \
    python3-setuptools \
    python3-venv

RUN cd /srv \
  && python3 -m venv venv \
  && . venv/bin/activate \
  && pip install --upgrade \
    jupyterlab \
    tensorflow

RUN apt-get upgrade -y
RUN apt autoremove --purge -y
RUN apt-get clean -y
RUN rm -rf /var/lib/apt/lists/*

COPY entrypoint.sh /srv
ENTRYPOINT "/srv/entrypoint.sh"
