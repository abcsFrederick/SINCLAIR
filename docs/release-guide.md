# Release Guide

## How to test a pre-release on biowulf

Install the development version of sinclair.

```sh
# activate the conda env for development
. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
conda activate py311

# go to the source on biowulf and update
cd /data/CCBR_Pipeliner/Pipelines/SINCLAIR/dev
git pull
# optionally switch to different branch if needed

# install the version to a hidden path (e.g. .v0.2.0-dev) in /data/CCBR_Pipeliner/Pipelines/SINCLAIR
cd ..
pip install ./dev -t ./.v0.2.0-dev
# add it to your PATH and PYTHONPATH with:
export PATH="$PATH:/data/CCBR_Pipeliner/Pipelines/SINCLAIR/.v0.2.0-dev/bin/"
export PYTHONPATH="$PYTHONPATH:/data/CCBR_Pipeliner/Pipelines/SINCLAIR/.v0.2.0-dev/"
```
