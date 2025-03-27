#!/bin/bash
module load StdEnv/2020
module load python/3.9.6
module load gcc/9.3.0
module load cuda/11.4
module load opencv/4.6.0
module load scipy-stack/2022a
source .objslam/bin/activate
pip install -r requirements.txt
