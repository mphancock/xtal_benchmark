import random
import pickle 
import sys 
import os     
import time 
import numpy as np 

import IMP 
import IMP.atom
import IMP.algebra 
import IMP.core 
import IMP.xtal

sys.path.append("/root/imp_devel/imp/modules/xtal/utility")
from structure_manipulation import add_waters
from structure_manipulation import perturb_model
from structure_manipulation import align_imp_models
from models import generate_imp_model_from_pdb
from models import generate_nth_imp_model_from_pdb
from models import generate_imp_model_copy
from model_statistics import compute_rmsd 
from model_statistics import compute_xtal_score


## generate pickle file containing score and rmsd of all models analyzed 
if __name__ == '__main__':
    IMP.set_log_level(0)
    IMP.atom.add_atom_type("HET: O  ", IMP.atom.ElementTable().get_element("O"))
 
    data_dir = "/wynton/home/sali/mhancock/job_output/"
    # data_dir = "/root/imp_devel/imp/modules/xtal/benchmark/data/"
    data_file = os.path.join(data_dir, "1ubq_md_score_vs_rmsd.p")

    cif_file = "/root/imp_devel/imp/modules/xtal/benchmark/data/reflections/1ubq_fexp.cif" 
    
    pdb_cryst = "/root/imp_devel/imp/modules/xtal/benchmark/data/pdb/1ubq.pdb" 
    pdb_nmr = "/root/imp_devel/imp/modules/xtal/benchmark/data/1d3z_aligned.pdb"  
    pdb_md = os.path.join(data_dir, "md_1")

    m_cryst_h20, h_cryst_h20 = generate_imp_model_from_pdb(pdb_cryst, IMP.atom.AllPDBSelector())
    m_cryst, h_cryst = generate_imp_model_from_pdb(pdb_cryst)

    data = {"score": [], "rmsd": []}

    for i in range(0,10000,10):
    # for i in range(0,10000,50): 
        t0 = time.time() 

        ## randomly perturb atoms (no bonds)
        # m_rand = generate_imp_model_copy(m_true)
        # perturb_model(m_rand, 10)

        ## use nmr structures 
        # m_nmr = generate_nth_imp_model_from_pdb(pdb_nmr, i)

        ## read in md structures 
        md_file = os.path.join(pdb_md, "1ubq_md_{}.pdb".format(i))
        ## md_file = "/root/imp_devel/imp/modules/xtal/benchmark/data/pdb/1ubq_md_sample.pdb"
        m_md, h_md = generate_imp_model_from_pdb(md_file, IMP.atom.AllPDBSelector())
        m_md, h_md = align_imp_models(m_cryst, m_md, h_md)


        ## temporary just to deal with missing atom 
        p = m_md.add_particle("OT")
        IMP.core.XYZ.setup_particle(m_md, p, IMP.algebra.Vector3D(40.862, 39.575, 36.251))
        a = IMP.atom.Atom.setup_particle(m_md, p, IMP.atom.AtomType("O"))
        a.set_temperature_factor(36.27)
        a.set_occupancy(.25)
        b = IMP.xtal.IsoBFactor.setup_particle(m_md, p, 36.27)

        # IMP.atom.write_pdb(h_md, "/root/imp_devel/imp/modules/xtal/benchmark/data/pdb/1ubq_md_sample_aligned.pdb")
    
        add_waters(m_md, m_cryst_h20)

        rmsd = compute_rmsd(m_cryst_h20, m_md) 
        score = compute_xtal_score(m_md, cif_file)
        t1 = time.time() 
        print(i, score, rmsd, t1-t0)

        data["score"].append(score) 
        data["rmsd"].append(rmsd)

    with open(data_file, 'wb') as f: 
        pickle.dump(data, f) 
