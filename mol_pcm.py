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
from tools import read_tpa,read_tpa_pcm,parse_tpa,exc_ladder,print_contributions, etpa_time,get_weight
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.signal import find_peaks
import matplotlib.gridspec as gridspec
import statistics

plt.rcParams.update({
    "font.family": "serif",
#    "font.sans-serif":['Helvetica'],
    "font.size": "11",
    "figure.facecolor":"white",
    "text.usetex":True,
    "pgf.texsystem" : "pdflatex",
    'savefig.dpi':300,
    'figure.dpi':100,
#    'figure.figsize':[6,5]
#    "text.latex.preamble":r"\usepackage{amsmath}",
})

mpl.style.use('seaborn-ticks')


def normalized(sigma):
    x_axis = np.around(np.real(sigma),decimals=0)  
    # Calculating mean and standard deviation
    mean = statistics.mean(x_axis)
    sd = statistics.stdev(x_axis)
    #plt.plot(x_axis, norm.pdf(x_axis, mean, sd),label=r"rhodamin,$\sigma= $"+ str(round(sd,0)))  
    return mean,sd

def find_locmax(sigma,height=1000,distance=200):
    peaks, prop = find_peaks(sigma,height=height,distance=distance)
    peaks = np.multiply(peaks,0.05) # it gives back the indices, for real value
    print("The peaks position(in fs:) {} \nwith the height{} ".format(peaks, prop["peak_heights"]))
    
def find_locmin(sigma,order):
    from scipy.signal import argrelextrema
    mins =argrelextrema(sigma, np.less,order=order)
    mins = np.multiply(mins,0.05) # it gives back the indices, for real value
    print("The local minimum(in fs:)",mins)  

def plot_etpa(T_e, pd_iso,s_iso, tot_states,sigma_iso, label='pyrimidine', domi=True):
    ### plot it
    #print(normalized(sigma_iso))
    #print(normalized(sigma_sup))
    fig,ax = plt.subplots(1,1) 
    plt.yscale("log")  
    plt.plot(T_e, sigma_iso,label=label,color='C2')
    plt.grid(which='major',linestyle = '--',axis='both',linewidth = 1)
    ax.set_ylim(bottom=1)
    plt.legend(loc='lower right')
    ax.tick_params(axis='both', which='major', labelsize=12)
    fig.text(0.5, 0.03, r"$T_e$ (fs)",fontsize=11,ha='center', va='center')
    fig.text(0.05, 0.5, "ETPA cross section (a.u.)", fontsize=11,ha='center', va='center', rotation='vertical') 
    plt.show()
    
def plot_etpa_compare(T_e, pd_iso,pd_sup, s_iso, s_sup, tot_states,sigma_iso, sigma_sup,delta,domi=True):
    ### plot it
    #print(normalized(sigma_iso))
    #print(normalized(sigma_sup))
    fig = plt.figure()
    gs = gridspec.GridSpec(2,3)  
    ax1 = fig.add_subplot(gs[0,:2])
    plt.yscale("log")  
    plt.plot(T_e, sigma_iso,label='pyrimidine')
    plt.grid(which='major',linestyle = '--',axis='both',linewidth = 1)
    ax1.set_ylim(bottom=1)
    plt.legend(loc='lower right')
    ax2 = fig.add_subplot(gs[1,:2],sharey=ax1)
    plt.yscale("log")  
#    plt.plot(T_e, sigma_sup,color='C2',label='pyrimidine in water')
    plt.plot(T_e, sigma_sup,color='C2',label=' embedded pyrimidine')
    plt.grid(which='major',linestyle = '--',axis='both',linewidth = 1)
    plt.legend(loc='lower right')
    ax3 = fig.add_subplot(gs[:,2])
    exc_ladder(pd_iso,pd_sup,s_iso=s_iso,s_sup=s_sup,tot_states=tot_states,title=" ",dominant=domi,delta=delta)
    plt.subplots_adjust(wspace=0.5)
    fig.text(0.35, 0.03, r"$T_e$ (fs)",fontsize=11,ha='center', va='center')
    fig.text(0.05, 0.5, "ETPA cross section (a.u.)", fontsize=11,ha='center', va='center', rotation='vertical') 
    plt.show()
# Define parameters
#root = os.getcwd()
#path = "/home/fu/projects/scn/adc2"
#path = "/home/fu/projects/TPA/sos_tpa/d2h_eth/keep_two_eth_v"
path = "/home/fu/projects/TPA/sos_tpa/d2h_h2o/keep_two_h2o_v"
polarization = "parallel"
T_e = np.arange(0,100.,0.05)
#T_e = [5.2]
#mol = "c4h4n2-d2h-h2o/acdz"
#mol= "d2h_eth/clusters"
#mol = "eth_h2o/H-N-H"
mol = "pp"
#mol = "9a/acdz"
#mol = "9a_4h2o"
#mol = "9a_Dichloro"
#mol = 'adc3'
pd_iso = read_tpa(os.path.join(path,"AB_MP2"))
#pd_sup = read_tpa(os.path.join(root,mol,'supp_dens'))#,correction=False)
#pd_iso = read_tpa(os.path.join(root,mol,'suppostion_nonPrepol'))
#pd_iso = read_tpa(os.path.join(root,mol,'frozen'))
#pd_iso = read_tpa(os.path.join(root,mol,'pcm/iso'))#,correction=True)
#pd_sup = read_tpa_pcm(os.path.join(root,mol,'pcm/sup'),correction=False)
s_iso = 1
s_sup = 10
#s_iso = 10
#s_sup = 11
tot_states =8
# =============================================================================
# get the dominant intermediate states
#domi_iso = print_contributions(pd_iso,considered_state=s_iso, tot_states=tot_states)
#domi_sup = print_contributions(pd_sup,considered_state=s_sup, tot_states=tot_states)
#get_weight(pd_sup,considered_state1=17,considered_state2=17,target_state=11)
#get_weight(pd_sup,considered_state1=16,considered_state2=16,target_state=11)
#get_weight(pd_sup,considered_state1=6,considered_state2=6,target_state=4)
#get_weight(pd_sup,considered_state1=4,considered_state2=6,target_state=4)
# =============================================================================
# Have a look at intermediate states contribution
#yticks= np.arange(4.,10.,1)
#exc_ladder(pd_iso,pd_sup,s_iso=s_iso,s_sup=s_sup,tot_states=tot_states,dominant=True)
# =============================================================================
#sigma_iso = etpa_time(T_e,pd_iso,considered_state=s_iso,tot_states=tot_states)#,exclude_state=iso_exclude)
#find_locmax(sigma_iso,height=max(sigma_iso)-2000, distance=50)    
###find_locmin(sigma_iso,order=50)
##print("the maximum etpa cross section for iso: \n",max(sigma_iso))
#print("the minimum etpa cross section for iso: \n",min(sigma_iso[1:]))
##print(50*"*")
#sigma_sup = etpa_time(T_e,pd_sup,considered_state=s_sup,tot_states=tot_states)#,exclude_state=sup_exclude)
#find_locmax(sigma_sup,height=max(sigma_sup)-2000,distance=50)
##print("the maximum etpa cross section for sup: \n",max(sigma_sup))
#print("the minimum etpa cross section for sup: \n",min(sigma_sup[1:]))
###find_locmin(sigma_sup,order=50)
##plot_etpa(T_e, pd_iso, s_iso,tot_states, sigma_iso,label="S10 of the free pyrimidine",domi=True)
#plot_etpa_compare(T_e, pd_iso,pd_sup, s_iso, s_sup, tot_states, sigma_iso, sigma_sup,delta=False, domi=True)