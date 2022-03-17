#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:14:17 2019

@author: fu
"""
import os
import pandas as pd
import numpy as np
import CCJob as ccj
from CCParser.QChem import ADC
from itertools import product
import matplotlib.pyplot as plt

def read_tpa(path):
    """
    input: the folder with the qchem output
    read:    ADC(2) TPA calulation
    return:  dataframe with TPA-related info
    """
    tpa = {'Excited state':[],'Excitation energy':[],'Oscillator strength':[],'dipole moment':[],\
           'Total dipole':[],'TPA cross section(sos)':[],'TPA cross section(dI)':[],"Transition dipole moment":[]}
    try:
        out = ccj.find_output(path)
    except AssertionError:
        pass
    else:
        with open(out,'r') as f:
            rl = f.readlines()
            es = ADC()
            # I need emb MP2 dipole for ground state
            gs_linenum = [i for i,line in enumerate(rl) if "MP(2) Summary" in line][-1]
            gs_dip = es.dipole_moment(gs_linenum+4,rl)
            tpa['Excited state'].append('0')
            tpa['dipole moment'].append(gs_dip)
            tpa['Total dipole'].append(float(rl[gs_linenum+5].split()[-1]))
            tpa['Excitation energy'].append(0)
            tpa['Oscillator strength'].append("-")
            tpa['TPA cross section(sos)'].append("-")
            tpa['TPA cross section(dI)'].append("-")
            tpa['Transition dipole moment'].append("-")
            for i,line in enumerate(rl):        
                matches = ["Excited state","[converged]"]
                if all(match in line for match in matches):                #if "Excited state" in line:
                     ex_num = line.split()[2] # ['Excited', 'state', '1'] 
                     ex_energy = es.exc_energies(i+5,rl)
                     ex_osi = es.osc_strength(i+7,rl)
                     ex_cro_sos = float(rl[i+10].split()[-1])/30
                     if "Two-photon absorption cross-section" in rl[i+15]:                    
                         ex_cro_dI = float(rl[i+15].split()[-1])/30
                         ex_cro_dI = round(ex_cro_dI,2)
                         ex_dip = es.dipole_moment(i+21,rl)
                         tot_dipole = float(rl[i+22].split()[-1])
                     else:
                         ex_cro_dI = "-"
                         ex_dip = es.dipole_moment(i+16,rl)
                         tot_dipole = float(rl[i+17].split()[-1])
                     #ex_amp = es.amplitudes(i+28,rl).to_dataframe()
                     ex_transDip = es.transition_dipole(i+8,rl)
                     tpa["Excited state"].append(ex_num)
                     tpa['Excitation energy'].append(ex_energy)
                     tpa['Oscillator strength'].append(round(ex_osi,3))
                     tpa['TPA cross section(sos)'].append(round(ex_cro_sos,2))
                     tpa['TPA cross section(dI)'].append(ex_cro_dI)
                     tpa['Transition dipole moment'].append(ex_transDip)
                     tpa['dipole moment'].append(ex_dip)
                     tpa['Total dipole'].append(tot_dipole)
                     #tpa['Transition amplitudes'].append(ex_amp.set_index('weight'))  
                if "Transition from excited state" in line:
                    ex_num_start = line.split()[4]
                    ex_num_end = line.split()[8]
                    ex_energy = es.exc_energies(i+3,rl)
                    ex_osi = es.osc_strength(i+4,rl)
                    ex_transDip = es.transition_dipole(i+5,rl)
                    tpa["Excited state"].append(ex_num_start+'->'+ex_num_end)
                    tpa['Excitation energy'].append(ex_energy)
                    tpa['Oscillator strength'].append(round(ex_osi,3))
                    tpa['TPA cross section(sos)'].append("-")
                    tpa['TPA cross section(dI)'].append("-")
                    tpa['Transition dipole moment'].append(ex_transDip)
                    tpa['dipole moment'].append("-")
                    tpa['Total dipole'].append("-")
        pd_tpa = pd.DataFrame(tpa).set_index('Excited state')
        pd_tpa.to_csv(os.path.join(path,'tpa.csv'))
        print("Saved tpa.csv in File!")
        return pd_tpa
    
def read_tpa_adc1(path):
    """
    input: the folder with the qchem output
    read:    ADC(1) TPA calulation
    return:  dataframe with TPA-related info
    """
    tpa = {'Excited state':[],'Excitation energy':[],'Oscillator strength':[],'dipole moment':[],\
           'Total dipole':[],'TPA cross section(sos)':[],"Transition dipole moment":[]}
    try:
        out = ccj.find_output(path)
    except AssertionError:
        pass
    else:
        with open(out,'r') as f:
            rl = f.readlines()
            es = ADC()
            # I need emb HF dipole for ground state
            gs_linenum = [i for i,line in enumerate(rl) if "HF Summary" in line][-1]
            gs_dip = es.dipole_moment(gs_linenum+3,rl)
            tpa['Excited state'].append('0')
            tpa['dipole moment'].append(gs_dip)
            tpa['Total dipole'].append(float(rl[gs_linenum+4].split()[-1]))
            tpa['Excitation energy'].append(0)
            tpa['Oscillator strength'].append("-")
            tpa['TPA cross section(sos)'].append("-")
            tpa['Transition dipole moment'].append("-")
            for i,line in enumerate(rl):        
                matches = ["Excited state","[converged]"]
                if all(match in line for match in matches):
                     print(line) 
                     ex_num = line.split()[2] # ['Excited', 'state', '1'] 
                     ex_energy = es.exc_energies(i+5,rl)
                     ex_osi = es.osc_strength(i+7,rl)
                     ex_cro_sos = float(rl[i+10].split()[-1])/30
                     ex_dip = es.dipole_moment(i+16,rl)
                     #ex_amp = es.amplitudes(i+28,rl).to_dataframe()
                     ex_transDip = es.transition_dipole(i+8,rl)
                     tpa["Excited state"].append(ex_num)
                     tpa['Excitation energy'].append(ex_energy)
                     tpa['Oscillator strength'].append(round(ex_osi,3))
                     tpa['TPA cross section(sos)'].append(round(ex_cro_sos,2))
                     tpa['Transition dipole moment'].append(ex_transDip)
                     tpa['dipole moment'].append(ex_dip)
                     tpa['Total dipole'].append(float(rl[i+17].split()[-1]))
                     #tpa['Transition amplitudes'].append(ex_amp.set_index('weight'))  
                if "Transition from excited state" in line:
                    try:
                        ex_num_start = line.split()[4]
                        ex_num_end = line.split()[8]
                        ex_energy = es.exc_energies(i+3,rl)
                        ex_osi = es.osc_strength(i+4,rl)
                        ex_transDip = es.transition_dipole(i+5,rl)
                        tpa["Excited state"].append(ex_num_start+'->'+ex_num_end)
                        tpa['Excitation energy'].append(ex_energy)
                        tpa['Oscillator strength'].append(round(ex_osi,3))
                        tpa['TPA cross section(sos)'].append("-")
                        tpa['Transition dipole moment'].append(ex_transDip)
                        tpa['dipole moment'].append("-")
                        tpa['Total dipole'].append("-")
                    except IndexError:
                        pass
        pd_tpa = pd.DataFrame(tpa).set_index('Excited state')
        pd_tpa.to_csv(os.path.join(path,'tpa.csv'))
        print("Saved tpa.csv in File!")
        return pd_tpa


def get_M(Cart_turple,pd_tpa,M, considered_state, intermediate_state):
    """
    Cart_turple: turple type eg. (0,0)
    calculate the TPA transition moment M matrix (3,3)
    one intermediate state(IS) case ---devided by E unit:eV
    """
    E_1 = pd_tpa.loc[str(considered_state),'Excitation energy']/2
    E_k = pd_tpa.loc[str(intermediate_state),'Excitation energy']
    E_diff = E_1 - E_k
    a, b = Cart_turple
    k = intermediate_state
    if k < considered_state:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] = (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/E_diff*27.211
    elif k > considered_state:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] = (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/E_diff*27.211
    else:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] = (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *(pd_tpa.loc[str(considered_state),'dipole moment'][b]-pd_tpa.loc['0','dipole moment'][b])\
                 + (pd_tpa.loc[str(considered_state),'dipole moment'][a]-pd_tpa.loc['0','dipole moment'][a])\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/E_diff*27.21   
    return M

def get_M_entangled(Cart_turple,pd_tpa,M, T_e, considered_state, intermediate_state):
    """
    Fro the entangled photons considering one inermediate state.
    Cart_turple: turple type eg. (0,0)
    calculate the TPA transition moment M matrix (3,3)
    one intermediate state(IS) case ---devided by E unit:eV
    """
    T_e = T_e * 41
    E_1 = pd_tpa.loc[str(considered_state),'Excitation energy']/2
    E_k = pd_tpa.loc[str(intermediate_state),'Excitation energy']
    E_diff = E_k - E_1
    phase = T_e*E_diff/27.211/2
    a, b = Cart_turple
    k = intermediate_state
    if k < considered_state:
            M[a][b] = (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])*2*np.sin(phase)/E_diff*27.211
    elif k > considered_state:
            M[a][b] = (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])*2*np.sin(phase)/E_diff*27.211
    else:
            M[a][b] = (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *(pd_tpa.loc[str(considered_state),'dipole moment'][b]-pd_tpa.loc['0','dipole moment'][b])\
                 + (pd_tpa.loc[str(considered_state),'dipole moment'][a]-pd_tpa.loc['0','dipole moment'][a])\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])*2*np.sin(phase)/E_diff*27.21   
    return M

def get_M_sos(Cart_turple,pd_tpa,M,considered_state, tot_state):
    """
    The SOS expression using fluctuation operator...
    Cart_turple: turple type eg. (0,0)         
    calculate the TPA transition moment M matrix
    sos intermediate state(IS) case --devided by E unit:eV=1/27 a.u. 
    """
    a, b = Cart_turple
    E_1 = pd_tpa.loc[str(considered_state),'Excitation energy']/2
    tot_state += 1
    for k in range(1,tot_state):
        E_k = pd_tpa.loc[str(k),'Excitation energy']
        if k < considered_state:
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.211
        elif k > considered_state:
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.211
        else:
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *(pd_tpa.loc[str(considered_state),'dipole moment'][b]-pd_tpa.loc['0','dipole moment'][b])\
                 + (pd_tpa.loc[str(considered_state),'dipole moment'][a]-pd_tpa.loc['0','dipole moment'][a])\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.21    
    #print(M)
    return M

def get_M_sos_entangled(Cart_turple,pd_tpa,M,T_e,considered_state, tot_state):
    """
    The SOS expression using fluctuation operator for entangled photons.
    Cart_turple: turple type eg. (0,0)         
    calculate the TPA transition moment M matrix
    sos intermediate state(IS) case --devided by E unit:eV=1/27 a.u. 1fs = 41a.u. 
    """
    T_e = T_e * 41
    a, b = Cart_turple
    E_1 = pd_tpa.loc[str(considered_state),'Excitation energy']/2
    tot_state += 1
    for k in range(1,tot_state):
        E_k = pd_tpa.loc[str(k),'Excitation energy']
        phase = T_e*(E_k-E_1)/27.211/2
        if k < considered_state:
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])*(2*np.sin(phase)*np.exp(1j*phase))/(E_k-E_1)*27.211
        elif k > considered_state:
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])*(2*np.sin(phase)*np.exp(1j*phase))/(E_k-E_1)*27.211
        else:
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *(pd_tpa.loc[str(considered_state),'dipole moment'][b]-pd_tpa.loc['0','dipole moment'][b])\
                 + (pd_tpa.loc[str(considered_state),'dipole moment'][a]-pd_tpa.loc['0','dipole moment'][a])\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])*(2*np.sin(phase)*np.exp(1j*phase))/(E_k-E_1)*27.21    
    #print(M)
    return M

def get_M_sos_randomOmega(Cart_turple,E_1,pd_tpa,M,considered_state, tot_state):
    """
    The SOS expression using fluctuation operator for random photon frequency...
    Cart_turple: turple type eg. (0,0)         
    calculate the TPA transition moment M matrix
    sos intermediate state(IS) case --devided by E unit:eV=1/27 a.u. 
    """
    a, b = Cart_turple
    tot_state += 1
    for k in range(1,tot_state):
        if k < considered_state:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.211
        elif k > considered_state:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.211
        else:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *(pd_tpa.loc[str(considered_state),'dipole moment'][b]-pd_tpa.loc['0','dipole moment'][b])\
                 + (pd_tpa.loc[str(considered_state),'dipole moment'][a]-pd_tpa.loc['0','dipole moment'][a])\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.21    
    #print(M)
    return M

def get_M_sos_zero(Cart_turple,pd_tpa,M,considered_state, tot_state):
    """
    This is another SOS expression where k includes the ground state.
    Cart_turple: turple type eg. (0,0)         
    calculate the TPA transition moment M matrix
    sos intermediate state(IS) case --devided by E unit:eV=1/27 a.u. 
    """
    a, b = Cart_turple
    E_1 = pd_tpa.loc[str(considered_state),'Excitation energy']/2
    tot_state += 1
    for k in range(1,tot_state):
        if k < considered_state:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(k)+'->'+str(considered_state)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.211
        elif k > considered_state:
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            M[a][b] += (pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state) + '->'+str(k)+':','Transition dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.211
        else:    
            E_k = pd_tpa.loc[str(k),'Excitation energy']
            E_0 = 0
            M[a][b] += (pd_tpa.loc['0','dipole moment'][a]*pd_tpa.loc[str(k),'Transition dipole moment'][b]\
                 +pd_tpa.loc['0','dipole moment'][b]*pd_tpa.loc[str(k),'Transition dipole moment'][a])/(E_0-E_1)*27.21\
                 +(pd_tpa.loc[str(k),'Transition dipole moment'][a]\
                 *pd_tpa.loc[str(considered_state),'dipole moment'][b]\
                 + pd_tpa.loc[str(considered_state),'dipole moment'][a]\
                 *pd_tpa.loc[str(k),'Transition dipole moment'][b])/(E_k-E_1)*27.21    
    #print(M)
    return M

def get_S(M):
    """
    input: M transition moment
    output: S tensor(9,3,3)
    """
    S = np.zeros((9,3,3))
    S[0,:,:] = M[0][0] * M
    S[1,:,:] = M[0][1] * M
    S[2,:,:] = M[0][2] * M
    S[3,:,:] = M[1][0] * M
    S[4,:,:] = M[1][1] * M
    S[5,:,:] = M[1][2] * M
    S[6,:,:] = M[2][0] * M
    S[7,:,:] = M[2][1] * M
    S[8,:,:] = M[2][2] * M
    return S

def get_sigma(M):
    """
    sigma = 1/30(F*sigma_f + G*sigma_g + H*sigma_h), F=G=H=2 for two parallel incident lights
    sigma_f = Suuvv,sigma_g = Suvuv,sigma_H = Suvvu
    unit: a.u.
    """
    sigma_f = 0
    sigma_g = 0 
    sigma_h = 0
    for i in range(3):
        for j in range(3):
            sigma_f += M[i][i]*M[j][j]
            #print(sigma_f)
            sigma_g += M[i][j]*M[i][j]
            #print(sigma_g)
            sigma_h += M[i][j]*M[j][i]
    #print("sigma_F:",sigma_f)
    #print("sigma_G:",sigma_g)
    #print("sigma_H:",sigma_h)
    sigma = (1/30)*(2*sigma_f + 2*sigma_g + 2*sigma_h)
    return sigma

def get_sigma_entangled(M,polarization='parallel'):
    """
    for entangled cases where phases are invloved.
    sigma = 1/30(F*sigma_f + G*sigma_g + H*sigma_h)
    sigma_f = Suuvv,sigma_g = Suvuv,sigma_H = Suvvu
    unit: a.u.
    """
    sigma_f = 0
    sigma_g = 0 
    sigma_h = 0
    for i in range(3):
        for j in range(3):
            sigma_f += M[i][i]*np.conjugate(M[j][j])
            #print('diagnol:',sigma_f)
            sigma_g += M[i][j]*np.conjugate(M[i][j])
            #print('ij*ij*',sigma_g)
            sigma_h += M[i][j]*np.conjugate(M[j][i])
            #print('ij*ji*',sigma_h)
    if polarization =='vertical':
        sigma = (1/30)*(-1*sigma_f + 4*sigma_g + -1*sigma_h)
    elif polarization =='parallel':
        sigma = (1/30)*(2*sigma_f + 2*sigma_g + 2*sigma_h)
    return sigma

def get_sigma_IScontri(pd_tpa,considered_state,intermediate_state):
    """
    input: dataframe
    output: the averaged TPA corss section using Three-state-model
    """
    cart = product(range(3), repeat = 2)
    M = np.zeros((3,3))
    for item in cart:    
        M_os = get_M(item,pd_tpa, M, considered_state, intermediate_state)
    sigma_os = get_sigma(M_os)
    return sigma_os

def plot_contributions(pd_tpa, basis, considered_state, tot_states):
    """
    input: dataframe, eg. tpa.csv
    basis: string, maybe not necessary in the graph
    ouput: plot the bar graph of the contributions from intemediate states of one
           considered state.
    """
    M = np.zeros((3,3))
    cart = product(range(3), repeat = 2)
    for item in cart:
        M_sos = get_M_sos(item, pd_tpa, M, considered_state,tot_states)
    sigma_sos = get_sigma(M_sos)
    tot_states += 1
    for j in range(1,tot_states):
        cart = product(range(3), repeat = 2)
        for item in cart:
             M_os = get_M(item,pd_tpa,M,considered_state,j)
        sigma = get_sigma(M_os)
        plt.bar(j, sigma, width=0.5, color='tab:blue')
        plt.xticks(np.arange(1,tot_states,1))
        #plt.title("Final TPA cross section for ES "+ str(considered_state) + ": " + str(round(sigma_sos,2)))
        #plt.ylabel("TPA cross section contribution (a.u.)")
        #plt.xlabel("Excited state (Three-state model)")
        #plt.tight_layout()
        
def plot_contributions_entangled(pd_tpa, basis,T_e, considered_state, tot_states):
    """
    input: dataframe, eg. tpa.csv
    basis: string, maybe not necessary in the graph
    ouput: plot the bar graph of the contributions from intemediate states of one
           considered state.
    """
    M = np.zeros((3,3),dtype=complex)
    cart = product(range(3), repeat = 2)
    for item in cart:
        M_sos = get_M_sos_entangled(item, pd_tpa, M,T_e, considered_state,tot_states)
    sigma_sos = get_sigma_entangled(M_sos)
    tot_states += 1
    for j in range(1,tot_states):
        cart = product(range(3), repeat = 2)
        for item in cart:
             M_os = get_M_entangled(item,pd_tpa,M,T_e,considered_state,j)
        sigma = get_sigma_entangled(M_os)
        plt.bar(j, sigma, width=0.5, color='tab:blue')
        plt.xticks(np.arange(1,tot_states,1))
        #plt.title("Microscopic ETPA cross section for ES "+ str(considered_state) + ": " + str(round(sigma_sos,2)))
        #plt.ylabel("TPA cross section contribution (a.u.)")
        #plt.xlabel(r"Excited state (Three-state model)  $T_e (fs)$ = "+str(round(T_e,2)))
        #plt.tight_layout()

def plot_contributions_e(pd_tpa, basis, considered_state, tot_states, ETPA = False):
    """
    input: dataframe, eg. tpa.csv
    basis: string, maybe not necessary in the graph
    ouput: plot the bar graph of the contributions from intemediate states of one 
           considered state.
    """
    tot_states += 1
    if ETPA == False:
        M = np.zeros((3,3))
        cart = product(range(3), repeat = 2)
        for item in cart:
            M_sos = get_M_sos(item, pd_tpa, M, considered_state,tot_states)
        sigma_sos = get_sigma(M_sos)        
        for j in range(1,tot_states):
            cart = product(range(3), repeat = 2)
            for item in cart:
                 M_os = get_M(item,pd_tpa,M,considered_state,j)
            sigma = get_sigma(M_os)
            plt.bar(j, sigma, width=0.5, color='tab:blue')
            print(j)
            plt.xticks(np.arange(1,11,1))
            #plt.title(basis + "---ES"+ str(considered_state) + ": " + str(round(sigma_sos,2)))
            plt.title("Final TPA cross section for ES "+ str(considered_state) + ": " + str(round(sigma_sos,2)))
            plt.ylabel("TPA cross section contribution (a.u.)")
            plt.xlabel("Excited state (Three-state model)")
            #plt.tight_layout()
    if ETPA == True:
        T_e = 5740
        M = np.zeros((3,3),dtype=complex)
        cart = product(range(3), repeat = 2)
        for item in cart:
            M_sos = get_M_sos_entangled(item, pd_tpa, M, T_e,considered_state,tot_states)
        sigma_sos = get_sigma_entangled(M_sos)        
        for j in range(1,tot_states):
            cart = product(range(3), repeat = 2)
            for item in cart:
                 M_os = get_M_entangled(item,pd_tpa,M,T_e,considered_state,j)
            sigma = np.real(get_sigma_entangled(M_os))
            plt.bar(j, sigma, width=0.5, color='tab:blue')
            plt.xticks(np.arange(1,11,1))
            #plt.title(basis + "---ES"+ str(considered_state) + ": " + str(round(sigma_sos,2)))
            plt.title("Final ETPA cross section for ES "+ str(considered_state) + ": " + str(round(np.real(sigma_sos),2)))
            plt.ylabel("TPA cross section contribution (a.u.)")
            plt.xlabel("Excited state (Three-state model)")
            #plt.tight_layout()        
        


def print_contributions(pd_tpa,considered_state, tot_states):
    """
    input: dataframe, eg. tpa.csv
    h_bar = 6.6*10**(-16) # unit: ev.s, here h_bar = 1--atomic units
    output: print the sos expression and all cotributions from each intemediate state
    """
    lines = 50*'-'
    M = np.zeros((3,3))
    # =============================================================================
    # 1.generate 9 combinations from [0,1,2] for M_a,b
    # =============================================================================
    cart = product(range(3), repeat = 2) # (0,0),(0,1),(0,2)
    # =============================================================================
    # 2.now generate M matrix
    # =============================================================================
    for item in cart:
        M_sos = get_M_sos(item, pd_tpa, M, considered_state,tot_states)
    print("M martrix for sos expression:\n",M_sos)
    sigma_sos = get_sigma(M_sos)
    ## =============================================================================
    # 3. generate the transition strength tensor S
    # =============================================================================
    #S = np.zeros((9,3,3))
    #S = get_S(M)
    ## =============================================================================
    print("The SOS TPA cross section for the " + str(considered_state) + "th state is:")
    print(sigma_sos)
    print(lines)
    # =============================================================================
    # Determine which intermediate state contributes the most
    # =============================================================================
    #sigma_sum = 0
    tot_states += 1
    #fig, axs = plt.subplots(2, 3, sharex=True, sharey=False)
    for j in range(1,tot_states):
        cart = product(range(3), repeat = 2)
        for item in cart:
             M_os = get_M(item,pd_tpa,M,considered_state,j)
        sigma = get_sigma(M_os)
        print("If only consider the " + str(j) + "th intermediate stata (Three-states model), the final cross section:")
        print(sigma)
        print("\n")
    #    sigma_sum += sigma
    #    print("Sum of previous states:  " + str(sigma_sum))
    #    print("\n")
        
def print_contributions_entangled(pd_tpa,T_e,considered_state, tot_states):
    """
    input: dataframe, eg. tpa.csv
    h_bar = 6.6*10**(-16) # unit: ev.s, here h_bar = 1--atomic units
    output: print the sos expression and all cotributions from each intemediate state
    """
    lines = 50*'-'
    M = np.zeros((3,3),dtype=complex)
    # =============================================================================
    # 1.generate 9 combinations from [0,1,2] for M_a,b
    # =============================================================================
    cart = product(range(3), repeat = 2) # (0,0),(0,1),(0,2)
    # =============================================================================
    # 2.now generate M matrix
    # =============================================================================
    for item in cart:
        M_sos = get_M_sos_entangled(item, pd_tpa, M, T_e,considered_state,tot_states)
    sigma_sos = get_sigma_entangled(M_sos)
    ## =============================================================================
    # 3. generate the transition strength tensor S
    # =============================================================================
    #S = np.zeros((9,3,3))
    #S = get_S(M)
    ## =============================================================================
    print("The SOS ETPA cross section for the " + str(considered_state) + "th state is:")
    print(sigma_sos)
    print(lines)
    # =============================================================================
    # Determine which intermediate state contributes the most
    # =============================================================================
    #sigma_sum = 0
    tot_states += 1
    #fig, axs = plt.subplots(2, 3, sharex=True, sharey=False)
    for j in range(1,tot_states):
        cart = product(range(3), repeat = 2)
        for item in cart:
             M_os = get_M_entangled(item,pd_tpa,M,T_e,considered_state,j)
        sigma = get_sigma_entangled(M_os)
        print("If only consider the " + str(j) + "th intermediate stata (Three-states model), the final cross section:")
        print(sigma)
        print("\n")
    #    sigma_sum += sigma
    #    print("Sum of previous states:  " + str(sigma_sum))
    #    print("\n")
    
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_xlabel(r"$|k\rangle$")
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)
    ax.set_ylabel(r"$|f\rangle$")
    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")
    ax.grid(visible=None)

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def heatmap_IS(pd_tpa,tot_states):
    #fig, ax = plt.subplots()
    M = np.zeros((3,3))
    three_state_model = []
    for i in range(1,tot_states+1):
        considered_state = i
        #print("final state:",i)
        for j in range(1,tot_states+1):
            cart = product(range(3), repeat = 2)
            for item in cart:
                 M_os = get_M(item,pd_tpa,M,considered_state,j)
                 sigma = get_sigma(M_os)
            three_state_model.append(sigma)
    np_three_state = np.array(three_state_model)
    data = np_three_state.reshape(tot_states,tot_states)
    x = ["{}".format(i) for i in range(1, tot_states+1)]
    #print(x)
    y = ["{}".format(i) for i in range(1,tot_states+1)]
    #print(y)
    im, _ = heatmap(data,y,x,vmin=0,cmap="magma_r",cbarlabel="intermediate TPA contribution (a.u.)")
    #annotate_heatmap(im, valfmt="{x:d}", size=7, threshold=20,textcolors=("red", "white"))