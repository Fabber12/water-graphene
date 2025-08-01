# Water-Graphene
Data and analysis scripts associated with the publication “Role of surface oxidation in enhancing heat transfer across graphene/water interface via Thermal Boundary Resistance modulation”. Includes LAMMPS input files and post-processing scripts (Python, MATLAB)

## Contact Angle (CA)

This section describes the workflow to compute the water–graphene contact angle via LAMMPS simulations and MATLAB post-processing.

### Overview

A series of three replica simulations (`0-replica`, `1-replica`, `2-replica`) are performed to improve statistical reliability. Each replica reads the final structure of the previous run, appends water molecules, and dumps atom coordinates at fixed intervals.

### Directory Structure

```python
lammps/CA/
        ├── relaxed_graph/              # Pre-simulated relaxed graphene structures at various oxidation degree
        |
        ├── 0-replica/                  # Base simulation folder
        |   ├── TERSOFF_forcefield.ff       # Forcefield
        │   ├── Wet.in*                     # LAMMPS input files
        │   └── water_box.data              # Water box to merge with relaxed graphene structure
        |
        ├── 1-replica/                  # Simulation folder
        |   ├── TERSOFF_forcefield.ff       # Forcefield
        |   └── Wet.in*                     # LAMMPS input files
        |
        └── 2-replica/                  # Simulation folder
            ├── TERSOFF_forcefield.ff       # Forcefield
            └── Wet.in*                     # LAMMPS input files
```

### Usage

1. Navigate into `lammps/CA/0-replica/` and run LAMMPS simulations (replace `X` with number of MPI processes):

   ```bash
   mpirun -np X lmp_mpi -in Wet.in
   ```
2. Repeat for `1-replica` and `2-replica` in sequence.
3. Post-process .dump output files by running MATLAB scripts:
    
    ```python
    post-processing/CA/
                    ├── CA_data_parser.m         # Reads dump files from lammps/CA/*-replica/, produces wet_*.mat 
                    └── CA_trend.m               # Loads wet_*.mat, compute CA mean ± SE, and plot trend with error bars

    ```
* Ensure that each simulation completes before invoking MATLAB scripts to guarantee all data is available for analysis.




## Density Profile (DP)

This section describes how to compute the water density profile across the graphene interface using a single LAMMPS simulation followed by MATLAB and python post‑processing.

### Overview

The trajectory (.dump) is first binned along the surface-normal axis—using molecular surfaces generated with EDTSurf to define the bin boundaries—to tally oxygen and hydrogen atoms in 0.1 Å layers, after which these counts are converted into a mass-density profile across the graphene–water interface.

### Directory Structure

```python
lammps/DP/
        ├─ TERSOFF_forcefield.ff            # Forcefield
        └─ Water-Graph_density.in*          # LAMMPS input files

post-processing/DP/
                ├─ MS/                      # Molecular surface generation
                |   ├─ remove_dump_lines.py             # Removes water molecules from .dump file
                |   ├─ PDB_conversion.tcl               # Converts .dump file single frames PDB files
                |   ├─ EDTSurf                          # Software used in surface-generator.sh (see https://zhanggroup.org/EDTSurf/)
                |   ├─ surface-generator.sh             # Generates .ply 3D molecular surfaces 
                |   └─ (surface-generator_parallel.sh)  # Generates .ply 3D molecular surfaces (with `GNU parallel`)
                │
                ├─ data_parser.m            # Counts particles between bins for each frame
                ├─ density_profile.m        # Plots water density profile
                └─ fastPlyRead.m            # Matlab function to read ply files
```

### Usage

1. **Run the LAMMPS simulation** (replace `X` with MPI ranks):

   ```bash
   cd lammps/DP
   mpirun -np X lmp_mpi -in Water-Graph_density.in
   ```

   *Tip*: it's possible to start from a relaxed data structure or generate your own structure by simulating `lammps/equilibration/`. Change the `read_data` line in `Water-Graph_density.in` accordingly.

2. **Generate Molecular surface**:
    ```bash
    cd post-processing/DP/MS
    python remove_dump_lines.py N           # N : oxidation degree (e.g., 20)
    vmd -dispdev text PDB_conversion.tcl    
    ./surface-generator.sh ms 1501          # ms: molecular surface | 1501: total number of frames (i.e. total number of PBD or ply files)
    ```
3. **Density profile**:
    ```bash
    cd post-processing/DP
    ```
    Open MATLAB and run: 
    ```
    data_parser.m
    ```
    ``` 
    density_profile.m
    ```
    
    
    Open matlab and run `data_parser.m` firstly and then `density_profile.m`.

* EDTSurf is a binary file for linux systems.
* If you want to speed up `surface-generator.sh` use `surface-generator_parallel.sh`. Be sure to have `GNU parallel` installed.
* Running `data_parser.m` **before** `density_profile.m` is mandatory, otherwise the counters file will be missing.
* Variables `maxZheight` and `thickness` in `data_parser.m` and `density_profile.m` **must match**. Default values are 15 Å and 0.1 Å respectively, producing 150 bins.



## Phonon Density of States (PDOS)

This section explains how to obtain and compare the phonon density of states of water and graphene through a velocity‑autocorrelation (VACF) analysis based on a LAMMPS simulation followed by a Python spectral post‑processing.

### Overview

A microcanonical production run records atomic velocities every femtosecond; the resulting VACF files for water and graphene are Fourier‑transformed to yield their single‑sided amplitude spectra, which are then averaged over 100 blocks and used to compute the spectral overlap factor.

### Directory Structure

```python
lammps/PDOS/
        ├─ TERSOFF_forcefield.ff        # Forcefield
        └─ Water‑Graph_pdos.in*         # LAMMPS input files

post‑processing/PDOS/               
                  └─ pdos.ipynb                   # Computes DOS and overlap factor
```

### Usage

1. **Run the LAMMPS simulation** (replace `X` with MPI ranks):

   ```bash
   cd lammps/PDOS
   mpirun -np X lmp_mpi -in Water-Graph_pdos.in
   ```

   *Tip*: it's possible to start from a relaxed data structure or generate your own structure by simulating `lammps/equilibration/`. Change the `read_data` line in `Water-Graph_density.in` accordingly.

2. **Compute the PDOS and spectral overlap** in Python:

   ```bash
    cd post-processing/PDOS
    ```

    Open the Jupyter notebook and run:

    ```text
    pdos.ipynb
    ```

   The notebook loads `vacf_water.{1..100}` and `vacf_graph.{1..100}`, performs the FFT, plots the averaged single‑sided spectra, and prints the overlap factor $S_{graph/water}$.


* Required Python libraries: `numpy`, `pandas`, `scipy`, `matplotlib`.
* Ensure all 100 VACF files are present before launching the notebook; otherwise adjust the loop range inside the first cell.
