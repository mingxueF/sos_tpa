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
from tools import print_contributions,etpa_time
#from cycler import cycler
from scipy.optimize import curve_fit
import os
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
#    'axes.prop_cycle':cycler(color='brgcmyk'),
#    'figure.figsize':[3,7]
#    "text.latex.preamble":r"\usepackage{amsmath}",
})
mpl.style.use('seaborn-ticks')

def exc_ladder(pd_iso,pd_sup,tot_states,s_iso,s_sup,delta=1,title=" ",yval="Excitation energy",dominant=False):
#    fig, ax = plt.subplots()
#    plt.figure(figsize=(3,4.5))
    i = 0
    x = [1,2]
    if dominant:
        domi_iso,sigma_iso_is,sigma_iso = print_contributions(pd_iso,s_iso, tot_states)
        domi_sup,sigma_sup_is,sigma_sup = print_contributions(pd_sup,s_sup, tot_states)
        z = [sigma_iso_is,sigma_sup_is]
        y = [pd_iso.loc[domi_iso,yval], pd_sup.loc[domi_sup,yval]]
        print(title+'\n',y)
    else:
        y = [pd_iso.loc["1":str(tot_states),yval],pd_sup.loc["1":str(tot_states),yval]]
        domi_iso,sigma_iso_is,sigma_iso = print_contributions(pd_iso,s_iso, tot_states)
        domi_sup,sigma_sup_is,sigma_sup = print_contributions(pd_sup,s_sup, tot_states)
        z = [sigma_iso,sigma_sup]
    for xe, ye in zip(x, y):
        plt.scatter([xe] * len(ye), ye,c=z[i],cmap="copper_r",s=2000,marker='_',linewidth=1.5)
        #plt.colorbar()
        if i == 0:
            for j,s in enumerate(domi_iso):
                plt.annotate("S"+s,xy=(xe,ye[j]),ha='center',size=10)
        elif i == 1:
            for j,s in enumerate(domi_sup):
                plt.annotate("S"+s,xy=(xe,ye[j]),ha='center',size=10) 
        i = i + 1       
#    plt.legend(loc="upper center",markerscale=0.5)
    #plt.yticks(np.arange(min(ye),max(ye),0.2))
    y_min = min(min(y[0]),min(y[1]))
    y_max = max(max(y[0]),max(y[1]))
    center = (y_max+y_min)/2
    plt.text(1.5,center, r"$\frac{\Delta_{iso}}{\Delta_{sup}}=$"+str(delta),ha='center', va='center',size=14)
    iso_label = r"$|f=$"+str(s_iso)+r"$\rangle_{iso}$"
    sup_label = r"$|f=$"+str(s_sup)+r"$\rangle_{sup}$"
    plt.xticks([1,2],[iso_label, sup_label])
    #plt.ylabel("Excitation energy (eV)")
    plt.grid(None)
    plt.grid(axis='y',linewidth = 0.5)
    plt.grid(which='both',linestyle = '--', axis='y',linewidth = 0.5)
    plt.title(title)
    #plt.tight_layout()
    plt.show()
    
def converter(instr):
     # when datafram is saved to csv, numpy array becomes string
    return np.fromstring(instr[1:-1],sep =' ') 

def parse_tpa(root,system,iso=""):
    path = os.path.join(root,system,iso,'tpa.csv')
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

def func(x,a,b):
    return a*x+b
#path = '/Users/mingxue/lab/tpa/sos_tpa'  
#pd_iso_1e = parse_tpa(path,'1e','iso')
#pd_sup_1e = parse_tpa(path,'1e','sup')
#pd_iso_1a = parse_tpa(path,'1a','iso')
#pd_sup_1a = parse_tpa(path,'1a','sup')
#pd_iso_9a = parse_tpa(path,'9a','iso')
#pd_sup_9a = parse_tpa(path,'9a','sup')
#pd_iso_c2v = parse_tpa(path,'c4h4n2-c2v-h2o','iso')
#pd_sup_c2v = parse_tpa(path,'c4h4n2-c2v-h2o')
#pd_iso_d2h = parse_tpa(path,'c4h4n2-d2h-h2o','iso')
#pd_sup_d2h = parse_tpa(path,'c4h4n2-d2h-h2o','20es')
#pd_iso_c2h4_1h2o = parse_tpa(path,'c2h4','acdz/A_MP2')
#pd_sup_c2h4_1h2o = parse_tpa(path,'c2h4','acdz/AB_MP2')
#pd_iso_c2h4_3h2o = parse_tpa(path,'c2h4-3h2o','acdz/A')
#pd_sup_c2h4_3h2o = parse_tpa(path,'c2h4-3h2o','acdz/AB')
# =============================================================================
# get the dominant intermediate states
#print_contributions(pd_sup_1a,considered_state=8, tot_states=10)
#domi = True
#exc_ladder(pd_iso_9a,pd_sup_9a,s_iso=5,s_sup=4,tot_states=10,title="9a",dominant=domi)
#fig = plt.figure()
#ax1 = fig.add_subplot(231)
#exc_ladder(pd_iso_9a,pd_sup_9a,s_iso=5,s_sup=4,tot_states=10,delta=6.00,title="9a",dominant=domi)
#ax2 = fig.add_subplot(232)
#exc_ladder(pd_iso_1a,pd_sup_1a,s_iso=8,s_sup=8,tot_states=10,delta=1.19,title="1a",dominant=domi)
#ax3 = fig.add_subplot(233)
#exc_ladder(pd_iso_c2v,pd_sup_c2v,s_iso=11,s_sup=12,tot_states=20,delta=3.11,title="c2v",dominant=domi)
#ax4 = fig.add_subplot(234)
#exc_ladder(pd_iso_d2h,pd_sup_d2h,s_iso=10,s_sup=11,tot_states=20,delta=2.74,title="d2h",dominant=domi)
#ax5 = fig.add_subplot(235)
#exc_ladder(pd_iso_c2h4_1h2o,pd_sup_c2h4_1h2o,s_iso=6,s_sup=6,tot_states=10,delta=1.01,title="c2h4-1h2o",dominant=domi)
##ax6 = fig.add_subplot(246)
##exc_ladder(pd_iso_c2h4_3h2o,pd_sup_c2h4_3h2o,s_iso=11,s_sup=12,tot_states=20,title="c2h4-3h2o",dominant=domi)
##ax7 = fig.add_subplot(247)
##exc_ladder(pd_iso_c2h4_5h2o,pd_sup_c2h4_5h2o,s_iso=11,s_sup=12,tot_states=20,title="c2h4-5h2o",dominant=domi)
#ax8 = fig.add_subplot(236)
#exc_ladder(pd_iso_1e,pd_sup_1e,s_iso=7,s_sup=6,tot_states=10,delta=0.9,title="1e",dominant=domi)
#fig.text(0.06, 0.5, "Excitation energy (eV)",fontsize =12, ha='center', va='center', rotation='vertical') 
# =============================================================================
# =============================================================================
# a plot of the ratio
x_axis = [1.0,2.1,2.9,1.2,1.1,6.0,1.1]
y_axis = [0.96,2.45,3.35,1.16,1.05,6.32,1.11]
mols = ['ethylene-$H_2O$','pyrazine$(D_{2h})$-$H_2O$','pyrazine$(C_{2v})$-$H_2O$',
        '7-hydroxyquinoline-$H_2O$ ','7-hydroxyquinoline-$NH_3$','benzaldehyde-$2H_2O$',
        'dimethylaminopyridinium cation-$4H_2O$']
popt, pcov = curve_fit(func,x_axis,y_axis)
print(popt)
fit_k = popt[0]
fit_b = popt[-1]
xvals = np.arange(0.5,6,0.1)
yvals = func(xvals,fit_k,fit_b)
plt.scatter(x_axis,y_axis,color='black')
plt.plot(xvals,yvals,color='C2',linestyle='--',linewidth=1.2)
plt.grid(which='both',linestyle = '--',linewidth = 1)
plt.annotate('ethylene-$H_2O$',xy=(1.0,0.96),xytext=(0.7,0.5),arrowprops=dict( arrowstyle="->",connectionstyle='arc3,rad=0.3' ))
plt.annotate('7-hydroxyquinoline-$H_2O$',xy=(1.2,1.16),xytext=(0.7,1.6),arrowprops=dict( arrowstyle="->",connectionstyle='arc3,rad=0.3' ))
plt.annotate('7-hydroxyquinoline-$NH_3$',xy=(1.1,1.05),xytext=(1.3,0.8),arrowprops=dict( arrowstyle="->",connectionstyle='arc3,rad=-0.3' ))
plt.annotate('dimethylaminopyridinium cation-$4H_2O$',xy=(1.1,1.11),xytext=(1.4,1.1),arrowprops=dict( arrowstyle="->",connectionstyle='arc3,rad=0.3'))
plt.annotate('pyrazine$(D_{2h})$-$H_2O$',xy=(2.1,2.45),xytext=(2.2,2.45))#,arrowprops=dict( arrowstyle="->" ))
plt.annotate('pyrazine$(C_{2v})$-$H_2O$',xy=(2.9,3.35),xytext=(3,3.35))#,arrowprops=dict( arrowstyle="->" ))
plt.annotate('benzaldehyde-$2H_2O$',xy=(6.0,6.32),xytext=(4.8,6.32))
plt.xlabel(r"$\frac{\Delta}{\Delta'}$",fontsize=14)
plt.ylabel(r"$\frac{P_{complexd}}{P_{free}}$",fontsize=14)
# =============================================================================
# latex table
#T_e = np.arange(0.,50.,0.05)
#tot_20 = ['c4h4n2-c2v/20es','c4h4n2-c2v-h2o','c4h4n2-d2h/20es','c4h4n2-d2h-h2o/20es']
#mols = ['1e/iso','1e/sup','1a/iso','1a/sup','9a/iso','9a/sup','c4h4n2-c2v/20es',
#       'c4h4n2-c2v-h2o','c4h4n2-d2h/20es','c4h4n2-d2h-h2o/20es','c2h4/acdz/A_MP2',
#       'c2h4/acdz/AB_MP2']
#states = [7,6,8,8,5,4,11,12,10,11,6,6]
#summary = {'molecules':[],'max ETPA':[],'mean':[],'sd':[]}
#for state,mol in zip(states,mols):
#    pd_tpa = parse_tpa(path,mol)
#    if mol in tot_20:
#        tot_states = 20
#        sigma = etpa_time(T_e,pd_tpa,considered_state=state,tot_states=tot_states)
#    else:
#        tot_states = 10
#        sigma = etpa_time(T_e,pd_tpa,considered_state=state,tot_states=tot_states)
#    summary['molecules'].append(mol)
#    summary['max ETPA'].append(max(sigma))
#    mean,sd = normalized(sigma) 
#    summary['mean'].append(mean)
#    summary['sd'].append(sd)
#pd_sum = pd.DataFrame(summary)
#pd_sum.to_csv('sum',sep='&')
    
