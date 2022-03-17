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
#pd_tpa_c4h4n2_c2v_10 = read_tpa(os.path.join(root,'c4h4n2-c2v','10es'))
#pd_tpa_c4h4n2_d2h_10 = read_tpa(os.path.join(root,'c4h4n2-d2h','10es'))
#pd_tpa_c4h4n2_c2v_20 = read_tpa(os.path.join(root,'c4h4n2-c2v','20es'))
#pd_tpa_c4h4n2_d2h_20 = read_tpa(os.path.join(root,'c4h4n2-d2h','20es'))
#pd_tpa_c4h4n2_d2h_h2o_es10 = read_tpa(os.path.join(root,'c4h4n2-d2h-h20'))
#pd_tpa_c4h4n2_d2h_h2o_es20 = read_tpa(os.path.join(root,'c4h4n2-d2h-h20','20es'))
#pd_tpa_c4h4n2_c2v_h2o_es20 = read_tpa(os.path.join(root,'c4h4n2-c2v-h2o'))
#pd_tpa_c4h4n2_d2h_h2o_MEes10 = read_tpa(os.path.join(root,'c4h4n2-d2h-h20','10es','ME'))
#pd_tpa_c4h4n2_d2h_h2o_MEes20 = read_tpa(os.path.join(root,'c4h4n2-d2h-h20','20es',"ME"))
# =============================================================================
# for printing stuff
#print_contributions(pd_tpa_c10h7cl,considered_state, tot_states)
#print_contributions_entangled(pd_tpa, 5740, considered_state, tot_states) 
#plot_contributions(pd_tpa_c10h7cl, "acdz", considered_state=2, tot_states=10)#, ETPA = False)
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
T_e = np.arange(0.,150.,0.05)
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
sigma_c4h4n2_c2v_es20 = etpa_time(T_e,pd_tpa_c4h4n2_c2v_20,considered_state=11,tot_states=20,polarization = "parallel")
sigma_c4h4n2_c2v_h2o_es20 = etpa_time(T_e,pd_tpa_c4h4n2_c2v_h2o_es20,considered_state=12,tot_states=20,polarization = "parallel")
#sigma_c4h4n2_d2h_es20 = etpa_time(T_e,pd_tpa_c4h4n2_d2h_20,considered_state=10,tot_states=20,polarization = "parallel")
# =============================================================================
# Creta a heatmap for all intermedaite states contribution
# =============================================================================
#heatmap_IS(pd_tpa_c4h4n2_c2v_h2o_es20,tot_states=20)
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
fig = plt.figure()  
ax1 = fig.add_subplot(211)
plt.yscale("log")  
plt.plot(T_e, sigma_c4h4n2_c2v_es20,label='c2v-s11')
 #ax1.set_ylim(bottom=1)
plt.legend()
ax2 = fig.add_subplot(212,sharey=ax1)#,sharey=ax1)
plt.yscale("log")  
plt.plot(T_e, sigma_c4h4n2_c2v_h2o_es20,color='C1',label='c2v_h2o-s12')
plt.legend()
#ax3 = fig.add_subplot(212,sharey=ax1)#,sharey=ax1)
#plt.yscale("log")  
#plt.plot(T_e, sigma_benzene,color='C2',label='c6H6-S3')
#plt.legend()
# #Set common labels
fig.text(0.5, 0.04, r"$T_e$ (fs)", fontsize =12,ha='center', va='center')
fig.text(0.06, 0.5, "ETPA cross section (a.u.)",fontsize =12, ha='center', va='center', rotation='vertical') 
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

    