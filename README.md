# Integrating histopathology and transcriptomics for spatial profiling of the tumor microenvironment: a melanoma case study

Our pipeline (SPoTLIghT) to derive spatial graph-based interpretable features from H&E (fresh-frozen, FF) tissue slides is packaged as a docker container.

![](spotlight_a.jpg)
![](spotlight_b.jpg)

## Run SPoTLIghT

1. Build the docker image as follows:

```bash
docker build -t run_spotlight:vfinal . 
```
2. Create the structure of the repository (in your working directory):
```bash
mkdir -p data_example/images
mkdir -p data/checkpoint/Retrained_Inception_v4
```
2. Add your FF histopathology slides into the `data_example/images` folder.
3. Download retrained models to extract the histopathological features, available from Fu et al., Nat Cancer, 2020 ([Retrained_Inception_v4](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BSST292)). 
Once you unzip the folder, extract the files to the `data/checkpoint/Retrained_Inception_v4/` folder.

4. Run the docker to execute the pipeline.

```bash
docker run \
-v $(pwd)/data/:/data:ro \
-v $(pwd)/data_example/:/data_example:ro \
-v $(pwd)/output_example/1_histopathological_features:/output_example/1_histopathological_features:rw \
-v $(pwd)/output_example/2_tile_level_quantification:/output_example/2_tile_level_quantification:rw \
-v $(pwd)/output_example/3_spatial_features:/output_example/3_spatial_features:rw \
run_spotlight:vfinal
```

The pipeline comprises the following steps:
1. Extract 1,536 histopathological features from Inception V4 model.
2. Predict tile level abundances for the different cell types.
3. Compute spatial features.

## Output documentation

SPoTLIghT generates the following output directory structure:

