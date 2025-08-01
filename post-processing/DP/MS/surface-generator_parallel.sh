#!/bin/bash

# Created by F. Tarulli, Politecnico di Torino, Italy – February 27, 2025
#       Process PDB frames with EDTSurf in parallel (using 'GNU parallel') for the given mode (sas,vws,ms).

# Usage: ./surface-generator_parallel.sh [sas|vws|ms] <tot frames>
# Default mode: ms

usage="Usage: $0 [sas|vws|ms] <total_frames>. Example: $0 ms 1501"

mode=${1:-ms}

if [ -z "$2" ]; then
    echo "Error: you must specify the total number of frames in the dump"
    echo "$usage"
    exit 1
fi
frames=$2

if [ "$mode" = "sas" ]; then
    s_val=2
    mkdir -p "ply_sas"
    output_prefix="ply_sas/sas_"
elif [ "$mode" = "vws" ]; then
    s_val=1
    mkdir -p "ply_vws"
    output_prefix="ply_vws/vws_"
elif [ "$mode" = "ms" ]; then
    s_val=3
    mkdir -p "ply_ms"
    output_prefix="ply_ms/ms_"
else
    echo "Mode unknown: $mode. Use 'sas', 'vws' or 'ms'."
    exit 1
fi

echo "Mode: $mode, option -s: $s_val"
echo "Output: ${output_prefix}N"

N=$(( frames - 1 ))
# 36 job in parallel (adjust if needed)
seq 0 $N | parallel -j36 '
  input_file="PDB/frame_{}.pdb";
  output_file="'"$output_prefix"'{}";
  echo "Processing ${input_file} -> ${output_file}";
  ./EDTSurf -i "${input_file}" -o "${output_file}" -s '"$s_val"' -p 1.4 -f 20.0
'

