# napari_viewer_with_auto_install.py
"""
napari_viewer_with_auto_install.py

Opens downsampled projection and reference mask images in Napari, 
with automatic installation of required Python packages.

@author: Hemanth Mohan and Shreyas Suryanarayana
Duke University School of Medicine
Date: July 14, 2025
"""
# %% Install required packages directly (no subprocess)
import sys
import pkg_resources

required = {
    "numpy",
    "PyYAML",
    "scikit-image",
    "napari[all]"
}

installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print(f"Installing missing packages: {missing}")
    try:
        from pip._internal import main as pipmain
    except ImportError:
        from pip import main as pipmain  # Fallback for older pip versions
    pipmain(['install', *missing])

# %% Import packages
from os.path import join, split
import os
import numpy as np
import skimage.io as io
import napari
import yaml

# %% Load config and images
configPath = 'config.yaml'
with open(configPath, 'r') as f:
    config = yaml.safe_load(f)

file_path = config["file_path"]
dataPath = join(split(file_path)[0], 'downsampled_standard_binary.tiff')
data = io.imread(dataPath).astype(np.uint8)

areaMaskPath = join(split(file_path)[0], 'downsampled_standard_areaMask.tiff')
areaMask = io.imread(areaMaskPath).astype(np.uint8)

# %% Start Napari
viewer = napari.Viewer()
viewer.theme = 'light'

viewer.add_image(data, name='projection', colormap='I Forest', rendering='mip', contrast_limits=(0, 1))
viewer.add_image(areaMask, name='refMap', colormap='plasma', rendering='mip', opacity=0.5, contrast_limits=(0, 5))

napari.run()
