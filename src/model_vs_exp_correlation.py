import matplotlib.pyplot as plt 
import os 
import numpy as np 

import IMP.xtal 


def plot(sft_exp, sft_calc, plots_dir): 
    f_exp = [] 
    f_c = []
    hkls = sft_exp.get_hkl() 

    for hkl in hkls: 
        sf_exp = sft_exp.get_structure_factor(hkl)
        sf_calc = sft_calc.get_structure_factor(hkl) 

        f_exp.append(sf_exp.amplitude())
        f_c.append(sf_calc.amplitude())

    m, b = np.polyfit(np.array(f_exp), np.array(f_c), 1)

    plt.xlabel("F_exp (R-factor: 0.1856)")
    plt.ylabel("IMP F_c (R-factor: 0.0031)")

    plt.scatter(f_exp, f_c)
    plt.plot(np.array(f_exp), m*np.array(f_exp) + b, "r")
    plt.text(100,1750, "m: {0:.5f}".format(m))
    plt.savefig(os.path.join(plots_dir, "calc_vs_exp_correlation.png"))


if __name__ == '__main__':
    input_dir = "/root/imp_devel/imp/modules/xtal/test/input"
    plots_dir = "/root/imp_devel/imp/modules/xtal/benchmark/plots"

    sft_exp = IMP.xtal.StructureFactorTable()  
    sft_calc = IMP.xtal.StructureFactorTable() 

    sft_exp.read_cif(os.path.join(input_dir, "1ubq_fexp.cif"))
    sft_calc.read_cif(os.path.join(input_dir, "1ubq_imp.cif"))
    plot(sft_exp, sft_calc, plots_dir)