#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 16:23:35 2022

@author: mingxue
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from tools import read_tpa,heatmap_IS,print_contributions, get_M_entangled,plot_contributions,get_sigma,get_M_sos_entangled,get_sigma_entangled
import matplotlib.pyplot as plt
from itertools import product
import matplotlib.pyplot as plt
from scipy.stats import norm
import statistics

plt.rcParams.update({
    "font.family": "serif",
#    "font.sans-serif":['Helvetica'],
    "font.size": "7",
    "figure.facecolor":"white",
    "text.usetex":True,
    "pgf.texsystem" : "pdflatex",
    'savefig.dpi':300,
    'figure.dpi':100,
#    'figure.figsize':[6,5]
#    "text.latex.preamble":r"\usepackage{amsmath}",
})

mpl.style.use('seaborn')

def converter(instr):
     # when datafram is saved to csv, numpy array becomes string
    return np.fromstring(instr[1:-1],sep =' ') 

def etpa_time(T_e,pd_tpa,considered_state,tot_states,polarization = "parallel"):
    """
    input: list of entanglement time
    ouput: list- etpa value
    """
    M_matrix = {}
    sigma_entangled = []
    for t in T_e:    
        #print(t)
        cart = product(range(3), repeat = 2)
        M = np.zeros((3,3),dtype = complex)
        for item in cart:
            get_M_sos_entangled(item,pd_tpa,M,t,considered_state, tot_states)
            #get_M_entangled(item,pd_tpa,M,t,considered_state, 1)
        M_matrix[t] = M
        sig = get_sigma_entangled(M,polarization)
        #print(sig)
        sigma_entangled.append(sig)
    sigma_entangled = np.around(np.real(sigma_entangled),decimals=2)
    print('finished the last one:',sigma_entangled[-1])
    return sigma_entangled

def parse_pdTPA(root,system,iso,basis=""):
    path = os.path.join(root,system,basis,iso,'tpa.csv')
    pd_path = pd.read_csv(path,converters={'Transition dipole moment':converter,\
                                           'dipole moment':converter},index_col = 'Excited state')
    return pd_path

# Define parameters
root = os.getcwd()
polarization = "parallel"
T_e = np.arange(0.,100.,0.05)
pd_tpa_c4h4n2_c2v_10 = read_tpa(os.path.join(root,'c4h4n2-c2v','10es'))
