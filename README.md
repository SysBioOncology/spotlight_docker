# Spatial Profiling of Tumors by Leveraging Imaging and Transcriptomics (SPoTLIghT)

<img align="left" width=100% src="https://github.com/olapuentesantana/spotlight_docker/blob/main/PipelineA.jpg">
<img align="left" width=100% src="https://github.com/olapuentesantana/spotlight_docker/blob/main/PipelineB.jpg">
<img align="left" width=100% src="https://github.com/olapuentesantana/spotlight_docker/blob/main/PipelineC.jpg">

First, build the docker image as follows:

```bash
docker build -t run_spotlight_example:v1 . 
```

Then, run the docker to execute the whole pipeline or execute separate steps as follows:

```bash
docker run \                                                  
-v $(pwd)/data/:/data:ro \
-v $(pwd)/data_example/:/data_example:ro \
-v $(pwd)/output_example/:/output_example:rw \
run_spotlight_example:v1
```

The pipeline comprises the following steps:
1. Extract 1,536 histopathological features from Inception V4 model.
2. Predict tile level abundances for the different cell types.
3. Compute spatial features.
