# -*- coding: utf-8 -*-
"""
Created on Sat Apr  5 15:28:57 2025
To be executed in napari environment with brainreg installed

@author: Hemanth Mohan and Shreyas Suryanarayana, Duke University School of Medicine

"""
import os
from PIL import Image
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import tkinter as tk
from tkinter import filedialog
from tqdm import tqdm
import subprocess
import argparse

# Function to downsample and save image
def downsample_image(input_folder, output_folder, filename, dsXY):
    try:
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, filename)
        with Image.open(input_path) as img:
            new_size = (img.width // dsXY, img.height // dsXY)
            downsampled_img = img.resize(new_size, Image.LANCZOS)
            downsampled_img.save(output_path)
        return f"Processed: {filename}"
    except Exception as e:
        return f"Error processing {filename}: {e}"

def run_downsampling(downsampleXY,downsampleZ,input_folder, output_folder):
    tiff_extensions = ('.tif', '.tiff')
    all_files = sorted(f for f in os.listdir(input_folder) if f.lower().endswith(tiff_extensions))

    tiff_files = all_files[::downsampleZ]  # Z downsampling

    func = partial(downsample_image, input_folder, output_folder, dsXY=downsampleXY)

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(func, f) for f in tiff_files]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Downsampling"):
            result = future.result()
            # print(future.result())
    return

def allenRegistration(output_folder,input_folder,downsampleXY,downsampleZ):
    input_data_folder = output_folder
    output_data_folder = os.path.join(input_folder,'downsampled_reg')
    if not os.path.exists(output_data_folder):
        os.makedirs(output_data_folder, exist_ok=True)
        
        voxel_x = 1.8 * downsampleXY # e.g., pixel size in X dimension
        voxel_y = 1.8 * downsampleXY # e.g., pixel size in Y dimension
        voxel_z = 4.0 * downsampleZ
        
        atlas_name = "allen_mouse_25um"
        
        data_orientation = "sal" ### sal
        
        n_free_cpus = 2 # Example: leave 2 cores free
        
        command = [
            "brainreg",
            str(input_data_folder), # Convert Path object to string for subprocess
            str(output_data_folder),     # Convert Path object to string for subprocess
            "-v",                   # Flag for voxel sizes
            str(voxel_z),           # Voxel Z (needs to be string)
            str(voxel_y),           # Voxel Y (needs to be string)
            str(voxel_x),           # Voxel X (needs to be string)
            "--atlas",              # Flag for atlas name
            atlas_name,
            "--orientation",        # Flag for data orientation
            data_orientation,
            "--n-free-cpus",        # Flag for CPU control
            str(n_free_cpus)        # Number of cores (needs to be string)
        ]
        
        print("-" * 30)
        print(f"Running Brainreg with command:")
        # Print the command in a readable format
        print(" ".join(command))
        print("-" * 30)
        
        result = subprocess.run(command, check=True, capture_output=True, text=True, encoding='utf-8')
        
        # Print the standard output from brainreg
        print("\n--- Brainreg Output ---")
        print(result.stdout)
        print("-----------------------\n")
        print(f"Brainreg registration completed successfully!")
        print(f"Results saved in: {output_folder}")
    else:
        print("registration already done. Delete reg folder to rerun.")
    return

# GUI + processing logic
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Downsample and register brain images')
    parser.add_argument('--downXY', type=int, default=4, help='Downsample rate for each image (XY dimensions)')
    parser.add_argument('--downZ', type=int, default=5, help='Downsample rate in depth (Z dimension)')
    # parser.add_argument('--voxel_size', type=float, default=1.8, help='Original voxel size in microns')
    # parser.add_argument('--voxel_z', type=float, default=4.0, help='Original voxel z size in microns')
    # parser.add_argument('--atlas', type=str, default='allen_mouse_25um', help='Atlas name for registration')
    # parser.add_argument('--orientation', type=str, default='sal', help='Data orientation')
    # parser.add_argument('--n_free_cpus', type=int, default=2, help='Number of CPUs to leave free')
    
    args = parser.parse_args()
     
    root = tk.Tk()
    root.withdraw()
    
    input_folder = filedialog.askdirectory(title="Select INPUT folder")
    if not input_folder:
        print("No folder selected.")
        
    output_folder = os.path.join(input_folder, 'downsampled')
    # downXY = 4  # Downsample rate for reach image
    # downZ = 5   # Downsample rate in depth
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
        run_downsampling(downsampleXY=args.downXY, downsampleZ=args.downZ,
                         input_folder=input_folder, output_folder=output_folder)
    else:
        print("Downsample folder already exists. Delete to rerun.")
    
    
    # registration to atlas
    allenRegistration(output_folder=output_folder,input_folder=input_folder,downsampleXY=args.downXY,downsampleZ=args.downZ)
    
if __name__ == "__main__":
    main()