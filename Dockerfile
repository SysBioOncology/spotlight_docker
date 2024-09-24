FROM ubuntu:20.04
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get install -y python3.8 && \
    ln -s /usr/bin/python3.8 /usr/bin/python && \
    apt-get install -y python3-pip && \
    apt-get install -y openslide-tools && \
    apt-get install -y python3-openslide && \
    apt-get install -y libgl1-mesa-glx && \
    apt-get install -y python3.8-dev && \ 
    apt-get install -y build-essential && \ 
    apt-get install -y pkg-config && \
    apt-get install -y python-dev && \
    apt-get install -y libhdf5-dev && \
    apt-get install -y libblosc-dev
    

# Set up python environment #
RUN apt install python3.8-venv

# Add nf-bin with all Python/R scripts
ENV VIRTUAL_ENV=/spotlight_venv
RUN python3 -m venv ${VIRTUAL_ENV}
ENV PATH="${VIRTUAL_ENV}/bin:$PATH"

# RUN . spotlight_venv/bin/activate
COPY ./env_requirements.txt ./
RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install --default-timeout=900 -r env_requirements.txt

# Set up directories #
# WORKDIR /nf-bin

ENV PATH="nf-bin:${PATH}"


