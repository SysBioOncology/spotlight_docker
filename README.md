# Spatial Profiling of Tumors by Leveraging Imaging and Transcriptomics (SPoTLIghT)

<img align="left" width=100% src="https://github.com/joan-yanqiong/spotlight/blob/master/PipelineA.png">
<img align="left" width=100% src="https://github.com/joan-yanqiong/spotlight/blob/master/PipelineB.png">
<img align="left" width=100% src="https://github.com/joan-yanqiong/spotlight/blob/master/PipelineC.png">

First, build the docker image as follows:

```bash
docker build -t run_spotlight_example:v1 . 
```

Then, run the docker to execute the whole pipeline or execute separate steps as follows:

```bash
docker run \                                                  
-v $(pwd)/spotlight_docker/data/checkpoint/Retrained_Inception_v4/:/data/checkpoint/Retrained_Inception_v4:ro \
-v $(pwd)/spotlight_docker/data_example/:/data_example:ro \
-v $(pwd)/spotlight_docker/output_example/:/output_example:rw \
run_spotlight_example:v1
```

The pipeline comprises the following steps:
1. Extract 1,536 histopathological features from Inception V4 model.
2. Predict tile level abundances for the different cell types.
3. Compute spatial features.
