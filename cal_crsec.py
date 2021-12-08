#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:14:17 2019

Example for studing C2H4_h2o tpa cross section
@author: fu
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from tools import read_tpa, print_contributions, print_contributions_entangled,plot_contributions,plot_contributions_entangled, get_sigma_IScontri,get_M_sos_randomOmega,get_sigma,get_M_sos_entangled,get_sigma_entangled
import matplotlib.pyplot as plt
from itertools import product
import matplotlib.pyplot as plt

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

root = os.getcwd()
filename=[f for f in os.listdir() if os.path.isdir(os.path.join(root,f))]   
filename.remove('__pycache__')
basis = ["acdz","actz","acqz","dacdz","dactz","dacqz"]
considered_state = {'A_MP2':[6,6,6,5,5,5],'emb_ME':[6,6,5,5,5,5],'emb_SE':[5,5,5,5,5,5],'AB_MP2':[6,6,6,6,6,6]}
intermediate_state2 = {'A_MP2':[3,3,3,4,3,3],'emb_ME':[3,3,3,4,3,3],'emb_SE':[4,3,3,4,3,3],'AB_MP2':[4,4,3,5,4,3]}
intermediate_state1 = {'A_MP2':[1,1,1,1,1,1],'emb_ME':[1,1,1,1,1,1],'emb_SE':[1,1,1,1,1,1],'AB_MP2':[2,2,1,2,1,1]}
method = "A_MP2"
tot_states = 10
# =============================================================================
# User difined parameters
# =============================================================================
#ex_specific = input("which excited state for TPA cross section:")
#num_states = input("The totol number of excited states for SOS expression:")
#intermediate_state = 3
# =============================================================================
# read output and generate a dataframe
# =============================================================================    
 #read from a csv file
#pd_tpa = pd.read_csv('tpa.csv',converters={'Transition dipole moment':converter},index_col = 'Excited state')
 #read from the an output file
# =============================================================================
# test session
# ============================================================================= 
considered_state = 6
#T_e = 5740 #(140 fs, 1fs= 41 a.u.)
pd_tpa = pd.read_csv(os.path.join(root,"acdz","AB_MP2","tpa.csv"),converters={'Transition dipole moment':converter,'dipole moment':converter},index_col = 'Excited state')
#pd_tpa = read_tpa(os.path.join(root,'acqz','emb_SE'))
#print_contributions(pd_tpa,considered_state, tot_states)
#print_contributions_entangled(pd_tpa, 5740, considered_state, tot_states) 
#plot_contributions(pd_tpa, "acdz", considered_state, tot_states)#, ETPA = False)
#plot_contributions_entangled(pd_tpa, "dacdz",T_e, considered_state, tot_states)#, ETPA = False)
# =============================================================================
# For random input photon frequency
# =============================================================================
#Omega = np.arange(2,10,0.5)
#omega_depend = []
#for omega in Omega:
#    cart = product(range(3), repeat = 2)
#    M = np.zeros((3,3))
#    for item in cart:
#        get_M_sos_randomOmega(item,omega,pd_tpa,M,considered_state, tot_states)
#    sigma_ran = get_sigma(M)
#    omega_depend.append(sigma_ran)
#plt.scatter(Omega,omega_depend) 
# =============================================================================
# For entablged photon cases
# =============================================================================
sigma_entangled = []
#T_e = 5740 #(140 fs, 1fs= 41 a.u.)
T_e = np.arange(0,1000,1)
M_matrix = {}
for t in T_e:    
    print(t)
    cart = product(range(3), repeat = 2)
    M = np.zeros((3,3),dtype = complex)
    for item in cart:
        get_M_sos_entangled(item,pd_tpa,M,t,considered_state, tot_states)
    M_matrix[t] = M
    sig = get_sigma_entangled(M)
    print(sig)
    sigma_entangled.append(sig)
T_e = T_e/41
plt.plot(T_e,np.real(sigma_entangled)) 
plt.title(r"ETPA cross sections for $C_2H_4$ vary with $T_e$")
plt.ylabel("TPA cross section (a.u.)")
plt.xlabel(r"$T_e$ (fs)")    
# =============================================================================
# plot the IS develops with the basis set size
# =============================================================================
#basis_evlo_IS1 = []
#i = 0
#for bas in basis:
#    path = os.path.join(root,bas,method,"tpa.csv")
#    pd_tpa = pd.read_csv(path,converters={'Transition dipole moment':converter,'dipole moment':converter},index_col = 'Excited state')
#    basis_evlo_IS1.append(get_sigma_IScontri(pd_tpa,considered_state[method][i],intermediate_state1[method][i]))
#    i += 1
#basis_evlo_IS2 = []
#i = 0
#for bas in basis:
#    path = os.path.join(root,bas,method,"tpa.csv")
#    pd_tpa = pd.read_csv(path,converters={'Transition dipole moment':converter,'dipole moment':converter},index_col = 'Excited state')
#    basis_evlo_IS2.append(get_sigma_IScontri(pd_tpa,considered_state[method][i],intermediate_state2[method][i]))
#    i += 1
#
#fig, axs = plt.subplots(2,1,sharex=True)
#axs[0].bar(basis, basis_evlo_IS1, width=0.2, color='tab:blue')
#axs[0].set_title("The Intermediate state 1")
#axs[1].bar(basis, basis_evlo_IS2, width=0.2, color='tab:blue')
#axs[1].set_title("The Intermediate state 2")
##plt.xticks(np.arange(1,7),basis)
##plt.tight_layout(rect=[0.02, 0.01, 0, 0.9])
##plt.ylabel("TPA cross section contribution (a.u.)")
#fig.text(0.35, 0.95, "Frozen B embedding (SE)", fontsize =15 , va='center')
#fig.text(0.04, 0.5, 'TPA cross section contribution (a.u.)', fontsize =15 , va='center', rotation='vertical')
#plt.savefig(os.path.join(root,'three-state-model',method,'IS_contribute'))
# =============================================================================
# =============================================================================
# plot one considered state develops with basis set
# =============================================================================
#i = 0
#fig = plt.figure()
#for bas in basis:
#    path = os.path.join(root,bas,method)
#    #pd_tpa = read_tpa(path)
#    pd_tpa = pd.read_csv(os.path.join(path,"tpa.csv"),converters={'Transition dipole moment':converter,'dipole moment':converter},index_col = 'Excited state')
#    fig.add_subplot(2,3,i+1)
#    plot_contributions(pd_tpa,bas,considered_state[method][i],tot_states)
#    i += 1
#    plt.tight_layout(rect=[0.07, 0, 1, 0.9]) # left, bottom,right,top default(0,0,1,1)
#fig.text(0.4, 0.95, r"$C_2H_4$ alone ADC2", fontsize =15 , va='center')    
#fig.text(0.04, 0.5, 'TPA cross section contribution (a.u.)', fontsize =15 , va='center', rotation='vertical')
#plt.savefig(os.path.join(root,'three-state-model',method,'sum_final_state'))
# =============================================================================
# =============================================================================
# generate tpa.csv through output
# =============================================================================
#for bas in basis:
#    path = os.path.join(root,bas,'A_MP2')
#    pd_tpa = read_tpa(path)


    