FROM ubuntu:20.04
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get install -y python3.8 && \
    ln -s /usr/bin/python3.8 /usr/bin/python && \
    apt-get install -y python3-pip && \
    apt-get install -y openslide-tools && \
    apt-get install -y python3-openslide && \
    apt-get install -y libgl1-mesa-glx

# Set up python environment #
RUN apt install python3.8-venv
RUN python3 -m venv /spotlight_venv
RUN . spotlight_venv/bin/activate
COPY ./env_requirements.txt ./
RUN pip3 install -r env_requirements.txt

# Set up directories #
RUN mkdir -p /project/Python/libs
WORKDIR /

ENV PYTHONPATH /project/Python/libs
