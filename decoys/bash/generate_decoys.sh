#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/xtal_benchmark/decoys/bash
#$ -e /wynton/home/sali/mhancock/xtal_benchmark/decoys/bash
#$ -j y
#$ -l h_rt=32:00:00
#$ -l mem_free=16G


source activate md
python "$HOME"/xtal_benchmark/decoys/src/md_simulation.py