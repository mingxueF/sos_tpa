#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 16:19:51 2022

@author: mingxue
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from tools import read_tpa,heatmap_IS, read_tpa_adc1,print_contributions, get_M_entangled,print_contributions_entangled,plot_contributions,plot_contributions_entangled, get_sigma_IScontri,get_M_sos_randomOmega,get_sigma,get_M_sos_entangled,get_sigma_entangled
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

# =============================================================================
# add subplots in the figure for the comparison
#fig = plt.figure()  
#ax1 = fig.add_subplot(221)
#plt.yscale("log")  
#plt.plot(T_e, sigma_c2h4_3h2o_A,label='c2h4-parallel')
#plt.legend()
#ax2 = fig.add_subplot(222,sharey=ax1)#,sharey=ax1)
#plt.yscale("log")  
#plt.plot(T_e, sigma_c2h4_3h2o_AB,color='C2',label='c2h4-3h2o-parallel')
#plt.legend()
#ax3 = fig.add_subplot(223,sharey=ax1)#,sharey=ax1)
#plt.yscale("log")  
#plt.plot(T_e, sigma_c2h4_5h2o_A,label='c2h4-parallel')
#plt.legend()
#ax4 = fig.add_subplot(224,sharey=ax1)#,sharey=ax1)
#plt.yscale("log")  
#plt.plot(T_e, sigma_c2h4_5h2o_AB,color='C2',label='c2h4-5h2o-parallel')
#plt.legend()
## Set common labels
#fig.text(0.5, 0.04, r"$T_e$ (fs)", fontsize =12,ha='center', va='center')
#fig.text(0.06, 0.5, "ETPA cross section (a.u.)",fontsize =12, ha='center', va='center', rotation='vertical') 
## to set the top title
#fig.text(0.45, 0.95, r"$C_2H_4$\_$H_2O$", fontsize =15 , va='center')
# =============================================================================

# =============================================================================
# Plot normalization distribution
# =============================================================================
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
## Calculating mean and standard deviation
#x_axis = np.around(np.real(sigma_c2h4_h2o_A),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),0.1), norm.pdf(np.arange(0.,max(x_axis),0.1), mean, sd),\
#         color='C0',label=r"c2h4-parallel,$\sigma= $"+ str(round(sd,0)))
#
#x_axis = np.around(np.real(sigma_c2h4_h2o_A_p),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         linestyle='dashed',color='C0',label=r"c2h4-perpendicular,$\sigma= $"+ str(round(sd,0)))
#plt.xlabel("ETPA cross section (a.u.)")
#plt.ylabel("density of probability")
#x_axis = np.around(np.real(sigma_c2h4_h2o_AB),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),0.1), norm.pdf(np.arange(0.,max(x_axis),0.1), mean, sd),\
#         color='C2',label=r"c2h4-h2o-parallel,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_h2o_AB_p),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C2',linestyle='dashed',label=r"c2h4-h2o-perpendicular,$\sigma= $"+ str(round(sd,0)))
#plt.legend()
##------------------------------
#ax2 = fig.add_subplot(312)
#x_axis = np.around(np.real(sigma_c2h4_3h2o_A),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C0',label=r"c2h4-parallel,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_3h2o_A_p),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         linestyle='dashed',color='C0',label=r"c2h4-perpendicular,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_3h2o_AB),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C2',label=r"c2h4-3h2o-parallel,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_3h2o_AB_p),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C2',linestyle='dashed',label=r"c2h4-3h2o-perpendicular,$\sigma= $"+ str(round(sd,0)))
#plt.legend()
##------------------------------
#ax3 = fig.add_subplot(313,sharey= ax2)
#x_axis = np.around(np.real(sigma_c2h4_5h2o_A),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C0',label=r"c2h4-parallel,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_5h2o_A_p),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         linestyle='dashed',color='C0',label=r"c2h4-perpendicular,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_5h2o_AB),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C2',label=r"c2h4-5h2o-parallel,$\sigma= $"+ str(round(sd,0)))
#x_axis = np.around(np.real(sigma_c2h4_5h2o_AB_p),decimals=0) 
#mean = statistics.mean(x_axis)
#sd = statistics.stdev(x_axis)
#plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),\
#         color='C2',linestyle='dashed',label=r"c2h4-5h2o-perpendicular,$\sigma= $"+ str(round(sd,0)))
#plt.legend()
#fig.text(0.5, 0.03, "ETPA cross section (a.u.)", fontsize =12,ha='center', va='center')
#fig.text(0.06, 0.5, "density of probability",fontsize =12, ha='center', va='center', rotation='vertical')
##-----------------------------------------------------------------------
#fig, axes = plt.subplots(1, 2)     
#axes[0].plot(T_e, sigma_entangled)
#axes[0].set_title("Normal scale")
#
#axes[1].plot(T_e, sigma_entangled)
#axes[1].set_yscale("log") #the log transformation
#axes[1].set_title("Logarithmic scale")
#
# Set common labels
#fig.text(0.5, 0.04, r"$T_e$ (fs)", fontsize =12,ha='center', va='center')
#fig.text(0.06, 0.5, "ETPA cross section (a.u.)",fontsize =12, ha='center', va='center', rotation='vertical') 
#fig.text(0.45, 0.95, r"$C_2H_4$\_$H_2O$", fontsize =15 , va='center')
#fig.text(0.45, 0.95, r"$C_2H_4$", fontsize =15 , va='center')    
#fig.text(0.45, 0.95, "rhodamin", fontsize =15 , va='center')    
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
# plot intermediate states contribution
# =============================================================================
#i = 1
#fig = plt.figure()
#while i <= tot_states:
#    #path = os.path.join(root,bas,method)
#    #pd_tpa = read_tpa(path)
#    #pd_tpa = pd.read_csv(path_c10h8,converters={'Transition dipole moment':converter,'dipole moment':converter},index_col = 'Excited state')
#    fig.add_subplot(5,2,i)
#    plot_contributions(pd_tpa_c10h7cl,'acdz',i,tot_states)
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