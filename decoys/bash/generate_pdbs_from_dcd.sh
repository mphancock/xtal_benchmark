#!/bin/bash
# Author: Matthew Hancock
# Date: 10/28/2021
# Description: Read in the MD trajectory from the DCD file and generate individual pdb files corresponding to the desired frames.


PROJ_DIR="$HOME/xtal_benchmark/decoys"
DATA_DIR="$PROJ_DIR/data"
SRC_DIR="$PROJ_DIR/src"

# Activate the md environment.
source activate md

# Generate the atom id indexes for mdconvert.
python "$SRC_DIR/generate_atom_selection.py" >"$DATA_DIR/processed/atoms_ids.txt"

mdconvert --output "$DATA_DIR/processed/1ubq.decoy.all.tmp.pdb" --chunk 1000 --force --stride 10 --atom_indices "$DATA_DIR/processed/atoms_ids.txt" --topology "$DATA_DIR/input/topology.prmtop" "$DATA_DIR/output/1ubq_st.dcd"

tail -n +2 "$DATA_DIR/processed/1ubq.decoy.all.tmp.pdb" >"$DATA_DIR/processed/1ubq.decoy.all.pdb"
rm "$DATA_DIR/processed/1ubq.decoy.all.tmp.pdb"

# Split the condensed pdb file from the trajectory into pdb individual files.
NUM_LINES_PER_ENTRY=663
NUM_LINES_TOTAL=$(wc -l "$DATA_DIR/processed/1ubq.decoy.all.pdb" | awk '{ print $1 }')
NUM_ENTRY=$((NUM_LINES_TOTAL / NUM_LINES_PER_ENTRY))
for i in $(seq 0 $((NUM_ENTRY-1)))
  do
    sed -n '1,1p;2q' "$DATA_DIR/processed/1ubq.decoy.all.pdb" >"$DATA_DIR/processed/1ubq.decoy.$i.pdb"
    START=$((2+i*NUM_LINES_PER_ENTRY))
    END=$((2+(i+1)*NUM_LINES_PER_ENTRY))
    END_MINUS=$((END-1))

    sed -n "$START,$END_MINUS p;$END q" "$DATA_DIR/processed/1ubq.decoy.all.pdb" >>"$DATA_DIR/processed/1ubq.decoy.$i.pdb"
  done

rm "$DATA_DIR/processed/1ubq.decoy.all.pdb"


