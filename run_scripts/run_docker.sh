#!/bin/bash

docker run \                              
-v /Users/Oscar/ownCloud2/SystemsImmunoOncology/spatial_localization/spotlight/spotlight_repo/spotlight/data_example/:/data_example:ro \
-v /Users/Oscar/ownCloud2/SystemsImmunoOncology/spatial_localization/spotlight/spotlight_repo/spotlight/output_example/:/output_example:rw \
run_spotlight_example:v1
