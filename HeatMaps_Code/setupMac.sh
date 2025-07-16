#!/bin/bash
# Initialize conda
eval "$(conda shell.bash hook)"

echo "Installing required environments"
# Create and set up napari environment
conda env remove --name napari-env --all
conda create -y -n napari-env -c conda-forge python=3.11
conda activate napari-env
conda install -y -c conda-forge napari pyqt
conda update -y napari
conda install -c conda-forge brainreg -y
conda install -c conda-forge niftyreg -y
pip install brainglobe
conda deactivate

# Create and set up allensdk environment
conda create -y --name allensdk python=3.10
conda activate allensdk
pip install allensdk
conda install -y -c conda-forge tqdm
conda deactivate