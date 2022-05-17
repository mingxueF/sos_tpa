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

def normalized(sigma):
    x_axis = np.around(np.real(sigma),decimals=0)  
    # Calculating mean and standard deviation
    mean = statistics.mean(x_axis)
    sd = statistics.stdev(x_axis)
    #plt.plot(x_axis, norm.pdf(x_axis, mean, sd),label=r"rhodamin,$\sigma= $"+ str(round(sd,0)))  
    return mean,sd



#def plot_norm(system,T_e,pd,polarization):
#    """
#    system: string
#    """
#    A,B = system.split('_')
#    sigma_A_p = etpa_time(T_e,pd,considered_state=6,tot_states=10,polarization)
#    #normalized_A = normalized("sig)
#    plt.plot(np.arange(0.,max(x_axis),1), norm.pdf(np.arange(0.,max(x_axis),1), mean, sd),label=system +"-"+ polarization+","+r"$\sigma= $"+ str(round(sd,0)))
#    plt.xlabel("ETPA cross section (a.u.)")
#    plt.ylabel("density of probability")
#    plt.legend()

root = os.getcwd()
filename=[f for f in os.listdir() if os.path.isdir(os.path.join(root,f))]   
filename.remove('__pycache__')
basis = ["acdz","actz","acqz","dacdz","dactz","dacqz"]
considered_state = {'A_MP2':[6,6,6,5,5,5],'emb_ME':[6,6,5,5,5,5],'emb_SE':[5,5,5,5,5,5],'AB_MP2':[6,6,6,6,6,6]}
intermediate_state2 = {'A_MP2':[3,3,3,4,3,3],'emb_ME':[3,3,3,4,3,3],'emb_SE':[4,3,3,4,3,3],'AB_MP2':[4,4,3,5,4,3]}
intermediate_state1 = {'A_MP2':[1,1,1,1,1,1],'emb_ME':[1,1,1,1,1,1],'emb_SE':[1,1,1,1,1,1],'AB_MP2':[2,2,1,2,1,1]}
method = "A_MP2"

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
T_e = 15 #(140 fs, 1fs= 41 a.u.)
path_rho = os.path.join(root,'rhodamin','acdz')
path_benzene = os.path.join(root,'benzene')
path_c10h8 = os.path.join(root,'c10h8')
path_c10h7cl = os.path.join(root,'c10h7cl')
path_c10h7f = os.path.join(root,'c10h7f')
path_c2h4_adc1 = os.path.join(root,"c2h4","dacqz","adc1","10es","A","tpa.csv")
tot_states = 10
considered_state = 2
basis = 'acqz'
#pd_c2h4_5h2o_A = parse_pdTPA(root,"c2h4_5h2o",'A')
#pd_c2h4_5h2o_AB = parse_pdTPA(root,"c2h4_5h2o",'ME')
#pd_c2h4_3h2o_A = parse_pdTPA(root,"c2h4_3h2o",'A')
#pd_c2h4_3h2o_AB = parse_pdTPA(root,"c2h4_3h2o",'ME')
#pd_c2h4_h2o_A = parse_pdTPA(root,"c2h4",'A_MP2',basis)
#pd_c2h4_h2o_AB = parse_pdTPA(root,"c2h4",'AB_MP2',basis)
#pd_c2h4_h2o_ME = parse_pdTPA(root,"c2h4",'emb_ME',basis)
#pd_tpa_rho = read_tpa_adc1(path_rho)
#pd_tpa_benzene = read_tpa(path_benzene)
#pd_tpa_c10h8= read_tpa(path_c10h8)
#pd_tpa_c10h7cl = read_tpa(path_c10h7cl)
#pd_tpa_c10h7f = read_tpa(path_c10h7f)
# =============================================================================
pd_tpa_9a_iso = parse_pdTPA(root,"9a",'iso')
pd_tpa_9a_sup = parse_pdTPA(root,"9a",'sup')
pd_tpa_1a_iso = parse_pdTPA(root,"1a",'iso')
pd_tpa_1a_sup = parse_pdTPA(root,"1a",'sup')
# =============================================================================
# for printing stuff
print_contributions(pd_tpa_9a_sup,considered_state=4, tot_states=10)
#print_contributions_entangled(pd_tpa, 5740, considered_state, tot_states) 
#plot_contributions(pd_tpa_9a_iso, "acdz", considered_state=5, tot_states=10)#, ETPA = False)
#plot_contributions_entangled(pd_c2h4_h2o_A, "dacdz",T_e, considered_state, tot_states)#, ETPA = False)
#fig = plt.figure()  
#ax1 = fig.add_subplot(211)
#plot_contributions(pd_c2h4_h2o_AB, "dacdz", considered_state, tot_states)
#
#ax2 = fig.add_subplot(212,sharey=ax1)#,sharey=ax1)
#plot_contributions(pd_c2h4_h2o_ME, "dacdz", 5, tot_states)

#ax3 = fig.add_subplot(313,sharey=ax1)#,sharey=ax1)
#plot_contributions_entangled(pd_c2h4_h2o_A, "dacdz",15, considered_state, tot_states)

#fig.text(0.5, 0.04, r"Excited state (Three-state model)", fontsize =12,ha='center', va='center')
#fig.text(0.06, 0.5, "TPA cross section (a.u.)",fontsize =12, ha='center', va='center', rotation='vertical') 
# =============================================================================
# For entablged photon cases
# =============================================================================
#T_e = 140 #(140 fs, 1fs= 41 a.u.)
#sigma_entangled = np.load("/Users/mingxue/lab/tpa/sos_tpa/Perpendicular/ent_S1_SOS_rhodamin_0-4000-0.1.npy")[:20000]
polarization = "parallel"
T_e = np.arange(0.,100.,0.05)
# =============================================================================
#sigma_c2h4_h2o_A_p = etpa_time(T_e,pd_c2h4_h2o_A,considered_state=6,tot_states=10,polarization = "vertical")
#sigma_c2h4_h2o_A_fine = etpa_time(T_e,pd_c2h4_h2o_A,considered_state=6,tot_states=10,polarization = "parallel")
#sigma_c2h4_h2o_AB = etpa_time(T_e,pd_c2h4_h2o_AB,considered_state=6,tot_states=10,polarization = "parallel")
#sigma_c2h4_h2o_AB_p = etpa_time(T_e,pd_c2h4_h2o_AB,considered_state=6,tot_states=10,polarization = "vertical")
#sigma_c2h4_3h2o_A = etpa_time(T_e,pd_c2h4_3h2o_A,considered_state=6,tot_states=10,polarization = "parallel")
#sigma_c2h4_3h2o_A_p = etpa_time(T_e,pd_c2h4_3h2o_A,considered_state=6,tot_states=10,polarization = "vertical")
#sigma_c2h4_3h2o_AB = etpa_time(T_e,pd_c2h4_3h2o_AB,considered_state=6,tot_states=10,polarization = "parallel")
#sigma_c2h4_3h2o_AB_p = etpa_time(T_e,pd_c2h4_3h2o_AB,considered_state=6,tot_states=10,polarization = "vertical")
#sigma_c2h4_5h2o_A = etpa_time(T_e,pd_c2h4_5h2o_A,considered_state=6,tot_states=10,polarization = "parallel")
#sigma_c2h4_5h2o_A_p = etpa_time(T_e,pd_c2h4_5h2o_A,considered_state=6,tot_states=10,polarization = "vertical")
#sigma_c2h4_5h2o_AB = etpa_time(T_e,pd_c2h4_5h2o_AB,considered_state=6,tot_states=10,polarization = "parallel")
#sigma_c2h4_5h2o_AB_p = etpa_time(T_e,pd_c2h4_5h2o_AB,considered_state=6,tot_states=10,polarization = "vertical")
#sigma_rho = etpa_time(T_e,pd_tpa_rho,considered_state=1,tot_states=9,polarization = "parallel")
# =============================================================================
#sigma_benzene = etpa_time(T_e,pd_tpa_benzene,considered_state=3,tot_states=10,polarization = "parallel")
#sigma_c10h8_es4 = etpa_time(T_e,pd_tpa_c10h8,considered_state=4,tot_states=10,polarization = "parallel")
#sigma_c10h7cl_es4 = etpa_time(T_e,pd_tpa_c10h7cl,considered_state=4,tot_states=10,polarization = "parallel")
#sigma_c10h8_es8 = etpa_time(T_e,pd_tpa_c10h8,considered_state=8,tot_states=10,polarization = "parallel")
#sigma_c10h7cl_es8 = etpa_time(T_e,pd_tpa_c10h7cl,considered_state=8,tot_states=10,polarization = "parallel")
#sigma_c10h7f_es8 = etpa_time(T_e,pd_tpa_c10h7f,considered_state=8,tot_states=10,polarization = "parallel")
# =============================================================================
#sigma_c4h4n2_c2v_es10 = etpa_time(T_e,pd_tpa_c4h4n2_c2v_10,considered_state=10,tot_states=10,polarization = "parallel")
#sigma_c4h4n2_d2h_es10 = etpa_time(T_e,pd_tpa_c4h4n2_d2h_10,considered_state=10,tot_states=10,polarization = "parallel")
#sigma_c4h4n2_d2h_h2o_es10 = etpa_time(T_e,pd_tpa_c4h4n2_d2h_h2o,considered_state=10,tot_states=10,polarization = "parallel")
#sigma_c4h4n2_d2h_h2o_es20 = etpa_time(T_e,pd_tpa_c4h4n2_d2h_h2o_es20,considered_state=11,tot_states=20,polarization = "parallel")
#sigma_c4h4n2_d2h_h2o_MEes20 = etpa_time(T_e,pd_tpa_c4h4n2_d2h_h2o_MEes20,considered_state=10,tot_states=20,polarization = "parallel")
#sigma_c4h4n2_c2v_es20_80fs = etpa_time(T_e,pd_tpa_c4h4n2_c2v_20,considered_state=11,tot_states=20,polarization = "parallel")
#sigma_c4h4n2_c2v_es20 = etpa_time(T_e,pd_tpa_c4h4n2_c2v_20,considered_state=11,tot_states=20,polarization = "parallel")
#sigma_c4h4n2_c2v_h2o_es20 = etpa_time(T_e,pd_tpa_c4h4n2_c2v_h2o_es20,considered_state=12,tot_states=20,polarization = "parallel")
#sigma_c4h4n2_d2h_es20 = etpa_time(T_e,pd_tpa_c4h4n2_d2h_20,considered_state=10,tot_states=20,polarization = "parallel")
# =============================================================================
#sigma_9a_iso = etpa_time(T_e,pd_tpa_9a_iso,considered_state=5,tot_states=10,polarization = "parallel")
#sigma_9a_sup = etpa_time(T_e,pd_tpa_9a_sup,considered_state=4,tot_states=10,polarization = "parallel")
#sigma_1a_iso = etpa_time(T_e,pd_tpa_1a_iso,considered_state=7,tot_states=10,polarization = "parallel")
#sigma_1a_sup = etpa_time(T_e,pd_tpa_1a_sup,considered_state=7,tot_states=10,polarization = "parallel")
# =============================================================================
# Creta a heatmap for all intermedaite states contribution
# =============================================================================
#heatmap_IS(pd_tpa_1a_sup,tot_states=10)
# =============================================================================
# =============================================================================
#plt.plot(T_e, sigma_benzene,label="benzene") 
#plt.legend()
#plt.yscale("log")
#entangled_mean = np.mean(sigma_c2h4_h2o_AB)
#print('the mean value of ETPA:',entangled_mean)
##plt.axhline(y=entangled_mean,linestyle=':',linewidth=1.5)
##plt.title(r"ETPA cross sections for $rhodamin$ S1 vary with $T_e$")
#plt.ylabel("ETPA microscopic cross section (a.u.)")
#plt.xlabel(r"$T_e$ (fs)")
#plt.show()
# =============================================================================
#fig = plt.figure()  
#ax1 = fig.add_subplot(211)
#plt.yscale("log")  
#plt.plot(T_e, sigma_9a_iso,label='iso-s5')
# #ax1.set_ylim(bottom=1)
#plt.legend()
#ax2 = fig.add_subplot(212,sharey=ax1)#,sharey=ax1)
#plt.yscale("log")  
#plt.plot(T_e, sigma_9a_sup,color='C1',label='sup-s4')
#plt.legend()
##ax3 = fig.add_subplot(212,sharey=ax1)#,sharey=ax1)
##plt.yscale("log")  
##plt.plot(T_e, sigma_benzene,color='C2',label='c6H6-S3')
##plt.legend()
# #Set common labels
#fig.text(0.5, 0.04, r"$T_e$ (fs)", fontsize =12,ha='center', va='center')
#fig.text(0.06, 0.5, "ETPA cross section (a.u.)",fontsize =12, ha='center', va='center', rotation='vertical') 



    