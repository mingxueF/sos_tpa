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
from tools import read_tpa,exc_ladder,print_contributions, etpa_time,get_weight
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.signal import find_peaks
import matplotlib.gridspec as gridspec
import statistics

plt.rcParams.update({
    "font.family": "serif",
#    "font.sans-serif":['Helvetica'],
    "font.size": "10",
    "figure.facecolor":"white",
    "text.usetex":True,
    "pgf.texsystem" : "pdflatex",
    'savefig.dpi':300,
    'figure.dpi':100,
#    'figure.figsize':[6,5]
#    "text.latex.preamble":r"\usepackage{amsmath}",
})

mpl.style.use('seaborn-ticks')

def converter(instr):
     # when datafram is saved to csv, numpy array becomes string
    return np.fromstring(instr[1:-1],sep =' ') 


def parse_tpa(path):
    path = os.path.join(path,'tpa.csv')
    pd_path = pd.read_csv(path,converters={'Transition dipole moment':converter,\
                                           'dipole moment':converter},index_col = 'Excited state')
    return pd_path

def normalized(sigma):
    x_axis = np.around(np.real(sigma),decimals=0)  
    # Calculating mean and standard deviation
    mean = statistics.mean(x_axis)
    sd = statistics.stdev(x_axis)
    #plt.plot(x_axis, norm.pdf(x_axis, mean, sd),label=r"rhodamin,$\sigma= $"+ str(round(sd,0)))  
    return mean,sd

def find_locmax(sigma,height=1000,distance=200):
    peaks, _ = find_peaks(sigma,height=height,distance=distance)
    peaks = np.multiply(peaks,0.05)
    print("The peaks position(in fs:)",peaks)
       
# Define parameters
root = os.getcwd()
polarization = "parallel"
T_e = np.arange(0,100.,0.05)
pd_iso = parse_tpa(os.path.join(root,'9a','iso'))
pd_sup = parse_tpa(os.path.join(root,'9a','ME'))
#pd_iso = parse_tpa(os.path.join(root,'c4h4n2-d2h-h2o','iso'))
#pd_sup = parse_tpa(os.path.join('c4h4n2-d2h-h2o','20es'))
domi = True
s_iso = 5
s_sup = 4
delta = 1.24
tot_states =10
# =============================================================================
# get the dominant intermediate states
domi_iso = print_contributions(pd_iso,considered_state=s_iso, tot_states=tot_states)
domi_sup = print_contributions(pd_sup,considered_state=s_sup, tot_states=tot_states)
# =============================================================================
# =============================================================================
# Have a look at intermediate states contribution
#exc_ladder(pd_iso,pd_sup,s_iso=s_iso,s_sup=s_sup,tot_states=tot_states,delta=delta,dominant=domi)
# =============================================================================
# check the etpa cross section
#sigma_iso = etpa_time(T_e,pd_iso,considered_state=s_iso,tot_states=tot_states)
#find_locmax(sigma_iso,height=700,distance=100)
print("the maximum etpa cross section for iso: \n",max(sigma_iso))
#sigma_sup = etpa_time(T_e,pd_sup,considered_state=s_sup,tot_states=tot_states)
#find_locmax(sigma_sup,height=400,distance=500)
print("the maximum etpa cross section for sup: \n",max(sigma_sup))
## =============================================================================
### plot it
#print(normalized(sigma_iso))
#print(normalized(sigma_sup))
fig = plt.figure()
gs = gridspec.GridSpec(2,3)  
ax1 = fig.add_subplot(gs[0,:2])
plt.yscale("log")  
plt.plot(T_e, sigma_iso,label='chromophore')
 #ax1.set_ylim(bottom=1)
plt.legend(loc='lower right')
ax2 = fig.add_subplot(gs[1,:2])#,sharey=ax1)
plt.yscale("log")  
plt.plot(T_e, sigma_sup,color='C2',label='FDET (frozen B)')
plt.legend(loc='lower right')
ax3 = fig.add_subplot(gs[:,2])
exc_ladder(pd_iso,pd_sup,s_iso=s_iso,s_sup=s_sup,tot_states=tot_states,delta=delta,title=" ",dominant=domi)
plt.subplots_adjust(wspace=0.5)
fig.text(0.35, 0.04, r"$T_e$ (fs)",fontsize=11,ha='center', va='center')
fig.text(0.05, 0.5, "ETPA cross section (a.u.)", fontsize=11,ha='center', va='center', rotation='vertical') 
plt.show()
