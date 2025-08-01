# Water-Graphene
Data and analysis scripts associated with the publication “Role of surface oxidation in enhancing heat transfer across graphene/water interface via Thermal Boundary Resistance modulation”. Includes LAMMPS input files and post-processing scripts (Python, MATLAB)

## Contact Angle (CA)

This section describes the workflow to compute the water–graphene contact angle via LAMMPS simulations and MATLAB post-processing.

### Overview

A series of three replica simulations (`0-replica`, `1-replica`, `2-replica`) are performed to improve statistical reliability. Each replica reads the final structure of the previous run, appends water molecules, and dumps atom coordinates at fixed intervals.

### Directory Structure

```
lammps/CA/
        ├── relaxed_graph/         # Relaxed graphene structures at various oxidation degree
        |
        ├── 0-replica/             # Base simulation folder
        │   ├── Wet.in*              # LAMMPS input files
        │   └── water_box.data       # Water box to merge with relaxed graphene structure
        |
        ├── 1-replica/             # Simulation folder
        |   └── Wet.in*              # LAMMPS input files
        |
        └── 2-replica/             # Simulation folder
            └── Wet.in*              # LAMMPS input files

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
                    └── CA_trend.m               # Load wet_*.mat, compute CA mean ± SE, and plot trend with error bars

    ```
Ensure that each simulation completes before invoking MATLAB scripts to guarantee all data is available for analysis.



