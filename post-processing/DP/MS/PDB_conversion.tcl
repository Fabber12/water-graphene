mol new water-graph_reduced.dump type lammpstrj \
    first 0 last -1 step 1 waitfor -1

mkdir -p PDB
cd PDB

set nf [molinfo top get numframes]

for {set i 0} {$i < $nf} {incr i} {

set carb [atomselect top "type 3 or type 4" frame $i]
set oxy [atomselect top "type 5" frame $i]
set hyd [atomselect top "type 6" frame $i]

$carb set name C
$oxy set name O
$hyd set name H

set sel [atomselect top "all" frame $i]

# Write PDB file "frame_i.pdb"
$sel writepdb frame_$i.pdb

$sel delete

}

puts "Conversion completed! Saved $nf PDB files.

