# Author: Matthew Hancock
# Date: 10/28/2021
# Description: Generate the atom indices from the MD trajectory corresponding to the main chain heavy atoms as well as 58 hetero atoms recorded in the native ubiquitin structure.
from pathlib import Path
import mdtraj as md


data_dir = Path(Path.home(), "xtal_benchmark/decoys/data")
traj_file = Path(data_dir, "output/1ubq_st.dcd")
pdb_file = Path(data_dir, "input/topology.prmtop")
traj = md.load(
    filename_or_filenames=str(traj_file),
    top=str(pdb_file)
)

topology = traj.topology
# Select the main chain
main_chain_selection = topology.select("protein and type != H")

# Only keep first 58 water atoms (the # of hetero-atoms in pdb file).
hetero_atom_selection = topology.select("water and type == O")
hetero_atom_selection = hetero_atom_selection[:58]

selection = list(main_chain_selection)
selection.extend(list(hetero_atom_selection))
selection = [str(id) for id in selection]
print("\t".join(selection))


