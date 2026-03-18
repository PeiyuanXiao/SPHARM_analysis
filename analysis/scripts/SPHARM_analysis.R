library(reticulate) # Load bridge
library(tidyverse)

# Force-connect the F-drive engine has been ignited
use_python("F:/ANACONDA/envs/spharm_local/python.exe", required = TRUE)

py_run_string(r"(
import sys, numpy as np, trimesh, os, glob, pandas as pd
import pyshtools as pysh
from unittest.mock import MagicMock

# Interference shielding library
for m in ['pyvista', 'pyvista.core._vtk_core', 'open3d', 'igl']:
    sys.modules[m] = MagicMock()

# --- Path Configuration --- 
base_dir = 'H:/SDG_Lithic_Analysis'
module_dir = os.path.join(base_dir, 'SPHARM_modules')
model_dir = os.path.join(base_dir, '3D_models_cores')
output_csv = os.path.join(base_dir, 'SDG_Cores_Full_Metrics.csv')

if module_dir not in sys.path:
    sys.path.append(module_dir)

import mesh_processing, pca_align, spherical_harmonics

def analyze_full_metrics(f_path):
    # A. Load
    mesh = trimesh.load(f_path)
    if not mesh.is_watertight:
        mesh.fill_holes()
    
    # B. Compute Sphericity
    vol, area = mesh.volume, mesh.area
    sphericity = (np.pi**(1/3) * (6*vol)**(2/3)) / area
    
    # C. Run SPHARM pipeline
    rng = np.random.default_rng(seed=1999)
    indices = rng.choice(len(mesh.vertices), 12000, replace=False)
    norm_v = mesh_processing.normalize_mesh(mesh.vertices[indices])
    align_v = pca_align.robust_pca_alignment(norm_v, enforce_direction=True)
    
    # D. Extract spherical Harmonics (SH)
    sph_coords = spherical_harmonics.cartesian_to_spherical(align_v)
    grid_r = spherical_harmonics.spherical_interpolate(sph_coords[:,0], sph_coords[:,1], sph_coords[:,2], 256)
    grid_dh = pysh.SHGrid.from_array(grid_r)
    clm = grid_dh.expand()
    
    # E.  Extract power spectrum
    power_spectrum = clm.spectrum()
    
    # F. Compute SHE
    she = np.sum(power_spectrum[1:16]) / np.sum(power_spectrum[0:16])
    
    # G. Extract spherical harmonic Coefficients
    coeffs_flattened = clm.coeffs.flatten()
    
    # Construct return dictionary for convenient Pandas processing 
    data_row = {
        'ID': os.path.basename(f_path),
        'Sphericity': round(sphericity, 4),
        'SHE': round(she, 4)
    }
    
    # Expand power spectrum into multiple columns
    for i, p_val in enumerate(power_spectrum):
        data_row[f'Spec_{i}'] = p_val
        
    # Expand coefficients into multiple columns
    for j, c_val in enumerate(coeffs_flattened):
        data_row[f'Coeff_{j}'] = c_val
        
    return data_row

# --- Batch Execution ---
files = glob.glob(os.path.join(model_dir, "*.stl"))

print(f">>> Preparing to process {len(files)} specimens")

if os.path.exists(output_csv):
    os.remove(output_csv)

header_written = False

for i, f in enumerate(files):
    name = os.path.basename(f)
    try:
        row_data = analyze_full_metrics(f)
        
        # Write immediately and release after writing
        df_row = pd.DataFrame([row_data])
        df_row.to_csv(
            output_csv,
            mode='a',
            header=not header_written,
            index=False
        )
        header_written = True
        
        # Explicitly release memory
        del df_row, row_data
        
        print(f"[{i+1}/{len(files)}] {name} - Written")
    except Exception as e:
        print(f"[{i+1}/{len(files)}] Skipping {name}，reason: {e}")

import gc
gc.collect()

print(f"\n>>> Mission completed! Data saved to: {output_csv}")
)")






