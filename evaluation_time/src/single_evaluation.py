# Author: Matthew Hancock
# Date: 10/18/2021
# Description: This script runs a single evaluation of the xtal restraint using the Monte Carlo optimizer and 1UBQ.

import sys
import time
import numpy as np

import IMP
import IMP.xtal


sys.path.append("/home/matthew/imp_tools/model_tools/src")
sys.path.append("/home/matthew/imp_development/imp/modules/xtal/utility")
from models import generate_imp_model_from_pdb
from read_cif import get_structure_factor_table_from_cif

IMP.set_log_level(IMP.SILENT)
pdb_file = "/home/matthew/xtal_benchmark/data/pdb/1ubq.pdb"
m, h = generate_imp_model_from_pdb(pdb_file=pdb_file,
                                   pdb_selector=IMP.atom.AllPDBSelector())

uc = IMP.xtal.UnitCell([50.84,42.77,28.95,90,90,90])
sg = IMP.xtal.SpaceGroup("P 2ac 2ab")
calc = IMP.xtal.StructureFactorCalc(uc, sg, "IT1992")

cif_file = "/home/matthew/xtal_benchmark/data/reflections/1ubq_fmodel.cif"
Fo_table = get_structure_factor_table_from_cif(cif_file=cif_file)

r = IMP.xtal.Restraint(m, Fo_table, calc)
r.evaluate(False)

sf = IMP.core.RestraintsScoringFunction([r])
mc = IMP.core.MonteCarlo(m)

movers = list()
ps = m.get_particle_indexes()
for p in ps:
    if IMP.atom.Atom.get_is_setup(m, p):
        movers.append(IMP.core.BallMover(m, p, 1.0))

mc.set_scoring_function(sf)
mc.set_kt(1)
mc.add_movers(movers)

times = list()
for i in range(10):
    t1 = time.time()
    mc.optimize(1)
    t2 = time.time()
    times.append(t2-t1)

print(np.mean(times), np.std(times))