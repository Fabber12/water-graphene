<h1 align="center">Water - Graphene Oxide</h1>

<p align="center">
<a href="https://www.nature.com/articles/sdata201618" target="_blank">
  <img src="https://custom-icon-badges.demolab.com/badge/data-FAIR-blue?logo=database&logoColor=white" alt="FAIR data" />
</a>
<a href="https://www.linux.org/" target="_blank">
  <img src="https://custom-icon-badges.demolab.com/badge/OS-Linux-orange?logo=linux&logoColor=white" alt="Linux" />
  </a>
<a href="https://python.org" target="_blank">
  <img src="https://custom-icon-badges.demolab.com/badge/Python-3.11%2B-blue?logo=python&logoColor=white" alt="Python 3.11+" />
</a>
<a href="https://www.matlab.com/" target="_blank">
   <img src="https://img.shields.io/badge/MATLAB-R2024b-orange" alt="Matlab R2024b" />
</a>
<a href="https://moltemplate.org/" target="_blank">
  <img src="https://custom-icon-badges.demolab.com/badge/Moltemplate-2025.3.18-lightblue?logo=gear&logoColor=white" alt="Moltemplate 2025-3-18" />
</a>
<a href=".github/CONTRIBUTING.md" target="_blank">
  <img src="https://custom-icon-badges.demolab.com/badge/contributions-open-4cb849?logo=code-of-conduct&logoColor=white" alt="Contributions welcome" />
</a>
</p>
<p align="center">
<a href="LICENSE" target="_blank">
  <img src="https://custom-icon-badges.demolab.com/badge/license-CC--BY%204.0-lightgray?logo=law&logoColor=white" alt="License CC-BY 4.0" />
</a>
</p>

Data and analysis scripts associated with the publication "*Role of surface oxidation in enhancing heat transfer across graphene/water interface via Thermal Boundary Resistance modulation*". Includes LAMMPS input files (some generated via Moltemplate) and post-processing scripts (Python, MATLAB)
 

## Contents
- [Contact Angle (CA)](#contact-angle-ca)
- [Density Profile (DP)](#density-profile-dp)
- [Phonon Density of States (PDOS)](#phonon-density-of-states-pdos)
- [Thermal Boundary Resistance (TBR)](#thermal-boundary-resistance-tbr)


<p style="color: orange; font-size: small;"><em>For clarity, all paths used with the <code>cd</code> command in this guide should be interpreted as absolute.</em></p>

## Contact Angle (CA)

This section describes the workflow to compute the water–graphene contact angle via LAMMPS simulations and MATLAB post-processing.

### Overview

A series of three replica simulations (`0-replica`, `1-replica`, `2-replica`) are performed to improve statistical reliability. Each replica starts from the final structure of the previous run and provides the trajectory files.

### Directory Structure

```bash
lammps/CA/
        ├── relaxed_graph/              # Pre-simulated relaxed graphene structures at various oxidation degrees
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

1. Navigate into `lammps/CA/0-replica/` and modify `Wet.in` by adjusting the `Atom Definition Section`. Then run LAMMPS simulations (replace `X` with MPI ranks):

   ```bash
   mpirun -np X lmp_mpi -in Wet.in
   ```
2. Repeat for `1-replica` and `2-replica` in sequence.
3. Post-process trajectory output files by running MATLAB scripts:
    
    ```bash
    post-processing/CA/
                    ├── CA_data_parser.m         # Reads dump files from lammps/CA/*-replica/, produces wet_*.mat 
                    └── CA_trend.m               # Loads wet_*.mat, computes the mean contact angle (CA) ± standard error, and plots trend with error bars
    ```
> Notes
> - Ensure that each simulation completes before running MATLAB scripts to guarantee all data is available for analysis.
> - It's possible to leverage `lammps/TBR/add_OH.m` script to create a new graphene oxide and then run an equilibration simulation to generate your own relaxed graphene oxide configuration. (Paths and settings in the MATLAB code have to be adjusted).



## Density Profile (DP)

This section describes how to compute the water density profile across the graphene interface using a single LAMMPS simulation followed by MATLAB and python post-processing.

### Overview

The trajectory (.dump) is first binned along the surface-normal axis—using molecular surfaces generated with EDTSurf to define the bin boundaries—to tally oxygen and hydrogen atoms in 0.1 Å layers, after which these counts are converted into a mass-density profile across the graphene–water interface.

### Directory Structure

```bash
lammps/DP/
        ├─ TERSOFF_forcefield.ff                 # Forcefield
        └─ Water-Graph_density.in*               # LAMMPS input files

post-processing/DP/
                ├─ MS/                           # Molecular surface generation
                |   ├─ remove_dump_lines.py             # Removes water molecules from .dump file
                |   ├─ PDB_conversion.tcl               # Converts single frames of a .dump file into PDB files
                |   ├─ EDTSurf                          # Software used in surface-generator.sh (see https://zhanggroup.org/EDTSurf/)
                |   ├─ surface-generator.sh             # Generates .ply 3D molecular surfaces 
                |   └─ (surface-generator_parallel.sh)  # Generates .ply 3D molecular surfaces (with `GNU parallel`)
                │
                ├─ data_parser.m                 # Counts the number of particles within each bin for every frame
                ├─ density_profile.m             # Plots water density profile
                └─ fastPlyRead.m                 # MATLAB function
```

### Usage

1. **Run the LAMMPS simulation** (replace `X` with MPI ranks):

   ```bash
   cd lammps/DP
   mpirun -np X lmp_mpi -in Water-Graph_density.in
   ```

   *Note*: it's possible to use a pre-relaxed system (available in `lammps/transient/systems_relaxed/`) or generate your own structure by simulating `lammps/equilibration/` (see *Thermal Boundary Resistance (TBR)* section). Change the `read_data` command in `Water-Graph_density.in` accordingly.

2. **Generate Molecular surface**:
    ```bash
    cd post-processing/DP/MS
    python remove_dump_lines.py ../../../lammps/DP/water-graph_density.dump water-graph_reduced.dump N           # N : oxidation degree (e.g., 20)
    vmd -dispdev text PDB_conversion.tcl    
    ./surface-generator.sh ms 1501          # ms: molecular surface | 1501: total number of frames (i.e. total number of PDB or ply files)
    ```
3. **Density profile**:
    ```bash
    cd post-processing/DP
    ```
    Open MATLAB and launch: 
    ```
    data_parser.m
    ```
    ``` 
    density_profile.m
    ```

> Notes
> - EDTSurf is a binary file for Linux systems.
> - If you want to speed up `surface-generator.sh` use `surface-generator_parallel.sh`. Be sure to have `GNU parallel` installed.
> - Running `data_parser.m` **before** `density_profile.m` is mandatory, otherwise the particle count file will be missing.
> - Variables `maxZheight` and `thickness` in `data_parser.m` and `density_profile.m` **must match**. Default values are 15 Å and 0.1 Å respectively, producing 150 bins.



## Phonon Density of States (PDOS)

This section explains how to obtain and compare the phonon density of states of water and graphene through a velocity-autocorrelation (VACF) analysis based on a LAMMPS simulation followed by a Python spectral post-processing.

### Overview

A microcanonical production run records atomic velocities every femtosecond; the resulting VACF files for water and graphene are Fourier-transformed to yield their single-sided amplitude spectra, which are then averaged over 100 blocks and used to compute the spectral overlap factor.

### Directory Structure

```bash
lammps/PDOS/
        ├─ TERSOFF_forcefield.ff        # Forcefield
        └─ Water-Graph_pdos.in*         # LAMMPS input files

post-processing/PDOS/               
                  └─ pdos.ipynb         # Computes PDOS and overlap factor
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

   The notebook loads `vacf_water.{1..100}` and `vacf_graph.{1..100}`, performs the FFT, plots the averaged single-sided spectra, and prints the overlap factor $S_{graph/water}$.

> Notes
> - Required Python libraries: `numpy`, `pandas`, `scipy`, `matplotlib`.
> - Ensure all 100 VACF files are present before launching the notebook; otherwise adjust the loop range inside the first cell.




## Thermal Boundary Resistance (TBR)

This section explains how to set up LAMMPS runs that yield the raw quantities required to determine the Kapitza resistance at the water-graphene oxide interface.

### Overview

A pristine graphene sheet is functionalised with a user-defined fraction of hydroxyl groups using `add_OH.m`. The resulting functionalised sheet is combined with a water box and equilibrated to yield a relaxed configuration (`lammps/TBR/equilibration` folder). A temperature difference is then imposed between the graphene oxide and the water, allowing the graphene's energy and system's temperature evolution to be tracked for TBR evaluation (`lammps/TBR/transient` folder).

### Directory Structure

```bash
lammps/TBR/
        ├─ graphene.lt                     # Base graphene lattice (pristine)
        ├─ add_OH.m                        # Adds -OH groups to graphene sheet
        │
        ├─ equilibration/                  # Equilibrates water-graphene oxide system
        │   ├─ H2O_Compass.lt                   # COMPASS water model
        │   ├─ Water-Graph_equilibration.lt     # Generates LAMMPS input files to run the simulation 
        │   └─ TERSOFF_forcefield.ff             # Forcefield
        │
        └─ transient/                      # Production run to evaluate the Kapitza resistance
            ├─ systems_relaxed/                 
            │   └─ relaxed_water-graph_*.data         # Pre-relaxed data files (oxidation (%): 0,5,10,20,40,60,80)
            │
            ├─ Water-Graph_transient.in*         # LAMMPS input files
            └─ TERSOFF_forcefield.ff             # Forcefield
```

### Usage

1. **Build the system**:

   ```bash
   cd lammps/TBR
   ```
   Open MATLAB, adjust the `oxid` variable to specify the desired percentage of –OH group coverage on the graphene sheet, and run the script:
   ```
   add_OH.m
   ```
2. **Equilibrate the system** (replace `X` with MPI ranks):
   ```bash
   cd lammps/TBR/equilibration
   ```
      ```bash
   moltemplate.sh Water-Graph_equilibration.lt
   ```
   ```bash
   mpirun -np X lmp_mpi -in Water-Graph_equilibration.in
   ```

3. **Launch the production run**:
   ```bash
   cd lammps/TBR/transient
   # edit the read_data command in `Water-Graph_transient.in`
   mpirun -np X lmp_mpi -in Water-Graph_transient.in
   ```
> Notes
> - A Moltemplate installation is required.
> - It's possible to skip `1. Build the system` by using one of the ready-made files in `lammps/transient/systems_relaxed/`, where all the systems analysed in this study are available.
