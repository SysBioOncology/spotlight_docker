FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt-get install -y python3.8
RUN ln -s /usr/bin/python3.8 /usr/bin/python
RUN apt-get install -y python3-pip
RUN apt-get install -y openslide-tools
RUN apt-get install -y python3-openslide
RUN apt-get install -y libgl1-mesa-glx

# Set up python environment #
RUN apt install python3.8-venv
RUN python3 -m venv /spotlight_venv
RUN . spotlight_venv/bin/activate
COPY ./env_requirements.txt ./
RUN pip3 install -r env_requirements.txt

# Set up directories #
WORKDIR /

# Set up required data #
COPY Python/ Python
COPY run_scripts/ run_scripts

# ------------- #
# Run spotlight:
# ------------- #

# 1. Extract 1,536 histopathological features computed using PC-CHiP from Fu et al.

# CMD ["bash", "/run_scripts/1_extract_histopatho_features.sh", "/data_example/images", "/data_example", "/data/checkpoint/Retrained_Inception_v4/model.ckpt-100000" ,"FFPE", "/output_example"]

# 2. Compute cell-type probability maps using transfer learning model (TCGA-SKCM)

CMD ["bash", "/run_scripts/2_tile_level_cell_type_quantification.sh", "/output_example/features_format_parquet", "FFPE", "/output_example"]

# 3. Compute spatial graph-based interpretable features

# CMD ["bash", "run_scripts/3_compute_spatial_features.sh", "path.to.tile.quant", "path.to.output"]

