#!/bin/bash

set -e  # Exit on error

eval "$(conda shell.bash hook)"
TARGET_ENV="allensdk"

echo "Activating '${TARGET_ENV}' environment..."
conda activate "$TARGET_ENV"

echo "Using Python at: $(which python)"
echo "Python version: $(python --version)"

# List of required Python modules
REQUIRED_MODULES=("PIL" "numpy" "scipy" "matplotlib.pyplot" "tqdm" "skimage.io" "tifffile" "yaml" "allensdk")

# Map import names to pip install names
declare -A MODULE_MAP=(
    ["PIL"]="Pillow"
    ["numpy"]="numpy"
    ["scipy"]="scipy"
    ["matplotlib.pyplot"]="matplotlib"
    ["tqdm"]="tqdm"
    ["skimage.io"]="scikit-image"
    ["tifffile"]="tifffile"
    ["yaml"]="pyyaml"
    ["allensdk"]="allensdk"
)

for MODULE in "${!MODULE_MAP[@]}"; do
    echo "Checking for Python module: ${MODULE}"
    python -c "import ${MODULE}" 2>/dev/null || {
        echo "⚠️  Module '${MODULE}' not found. Installing '${MODULE_MAP[$MODULE]}'..."
        pip install "${MODULE_MAP[$MODULE]}"
    }
done


echo "✅ All required packages are installed."

echo "Running Python script: binarizeAndGenerateAreaMask.py"
python binarizeAndGenerateAreaMask.py

echo "Deactivating environment..."
conda deactivate

echo "✅ Script completed successfully."
