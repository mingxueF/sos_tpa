#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 10:19:41 2021

@author: mingxue
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from tools import print_contributions
#from cycler import cycler
import os
plt.rcParams.update({
    "font.family": "serif",
#    "font.sans-serif":['Helvetica'],
    "font.size": "7",
    "figure.facecolor":"white",
    "text.usetex":True,
    "pgf.texsystem" : "pdflatex",
    'savefig.dpi':300,
    'figure.dpi':100,
#    'axes.prop_cycle':cycler(color='brgcmyk'),
#    'figure.figsize':[3,7]
#    "text.latex.preamble":r"\usepackage{amsmath}",
})
#mpl.style.use('seaborn-whitegrid')

def exc_ladder(pd_iso,pd_sup,tot_states,s_iso,s_sup,title=" ",yval="Excitation energy",dominant=False):
#    fig, ax = plt.subplots()
    i = 0
    if dominant:
        domi_iso,sigma_iso_is,sigma_iso = print_contributions(pd_iso,s_iso, tot_states)
        domi_sup,sigma_sup_is,sigma_sup = print_contributions(pd_sup,s_sup, tot_states)
        z = [sigma_iso_is,sigma_sup_is]
        print(z)
        y = [pd_iso.loc[domi_iso,yval], pd_sup.loc[domi_sup,yval]]
        x = ['iso', 'sup']
        print(y)
    else:
        y = [pd_iso.loc["1":str(tot_states),yval],pd_sup.loc["1":str(tot_states),yval]]
        x = ['iso', 'sup']
        domi_iso,sigma_iso_is,sigma_iso = print_contributions(pd_iso,s_iso, tot_states)
        domi_sup,sigma_sup_is,sigma_sup = print_contributions(pd_sup,s_sup, tot_states)
        z = [sigma_iso,sigma_sup]
    for xe, ye in zip(x, y):
        plt.scatter([xe] * len(ye), ye,c=z[i],cmap="magma_r",s=2000,marker='_',linewidth=1,label=xe)
        #plt.colorbar()
        i = i + 1
        #plt.legend(loc="upper center",markerscale=0.5)
    plt.xticks(label=['iso', 'sup'])
    plt.grid(None)
    plt.grid(axis='y',linewidth = 0.5)
    plt.grid(which='major',linestyle = '--', axis='y',linewidth = 0.5)
    plt.title(title)
    plt.show()

def converter(instr):
     # when datafram is saved to csv, numpy array becomes string
    return np.fromstring(instr[1:-1],sep =' ') 

def parse_tpa(root,system,iso=""):
    path = os.path.join(root,system,iso,'tpa.csv')
    pd_path = pd.read_csv(path,converters={'Transition dipole moment':converter,\
                                           'dipole moment':converter},index_col = 'Excited state')
    return pd_path

path = '/Users/mingxue/lab/tpa/sos_tpa'  
pd_iso_1e = parse_tpa(path,'1e','iso')
pd_sup_1e = parse_tpa(path,'1e','sup')
pd_iso_1a = parse_tpa(path,'1a','iso')
pd_sup_1a = parse_tpa(path,'1a','sup')
pd_iso_9a = parse_tpa(path,'9a','iso')
pd_sup_9a = parse_tpa(path,'9a','sup')
pd_iso_c2v = parse_tpa(path,'c4h4n2-c2v','20es')
pd_sup_c2v = parse_tpa(path,'c4h4n2-c2v-h2o')
pd_iso_d2h = parse_tpa(path,'c4h4n2-d2h','20es')
pd_sup_d2h = parse_tpa(path,'c4h4n2-d2h-h2o','20es')
pd_iso_c2h4_1h2o = parse_tpa(path,'c2h4','acqz/A_MP2')
pd_sup_c2h4_1h2o = parse_tpa(path,'c2h4','acqz/AB_MP2')
# =============================================================================
# get the dominant intermediate states
#domi_iso = print_contributions(pd_iso_1e,considered_state=8, tot_states=10)
#domi_sup = print_contributions(pd_sup_1e,considered_state=8, tot_states=10)
#print(pd_iso_1e.loc[domi_iso,"Excitation energy"])
fig = plt.figure()
#exc_ladder(pd_iso_1e,pd_sup_1e,tot_states=10,title="9a",domi_iso=domi_iso,domi_sup=domi_sup)
domi = True
ax1 = fig.add_subplot(231)
exc_ladder(pd_iso_9a,pd_sup_9a,s_iso=5,s_sup=4,tot_states=10,title="9a",dominant=domi)
ax2 = fig.add_subplot(232)
exc_ladder(pd_iso_1a,pd_sup_1a,s_iso=8,s_sup=8,tot_states=10,title="1a",dominant=domi)
ax3 = fig.add_subplot(233)
exc_ladder(pd_iso_c2v,pd_sup_c2v,s_iso=11,s_sup=12,tot_states=20,title="c2v",dominant=domi)
ax4 = fig.add_subplot(234)
exc_ladder(pd_iso_d2h,pd_sup_d2h,s_iso=10,s_sup=11,tot_states=20,title="d2h",dominant=domi)
ax5 = fig.add_subplot(235)
exc_ladder(pd_iso_c2h4_1h2o,pd_sup_c2h4_1h2o,s_iso=6,s_sup=6,tot_states=10,title="c2h4-1h2o",dominant=domi)
#ax6 = fig.add_subplot(246)
#exc_ladder(pd_iso_c2h4_3h2o,pd_sup_c2h4_3h2o,s_iso=11,s_sup=12,tot_states=20,title="c2h4-3h2o",dominant=domi)
#ax7 = fig.add_subplot(247)
#exc_ladder(pd_iso_c2h4_5h2o,pd_sup_c2h4_5h2o,s_iso=11,s_sup=12,tot_states=20,title="c2h4-5h2o",dominant=domi)
ax8 = fig.add_subplot(236)
exc_ladder(pd_iso_1e,pd_sup_1e,s_iso=8,s_sup=8,tot_states=10,title="1e",dominant=domi)
