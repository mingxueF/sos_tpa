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

def converter(instr):
     # when datafram is saved to csv, numpy array becomes string
    return np.fromstring(instr[1:-1],sep =' ') 


def parse_tpa(path):
    path = os.path.join(path,'tpa.csv')
    pd_path = pd.read_csv(path,converters={'Transition dipole moment':converter,\
                                           'dipole moment':converter},index_col = 'Excited state')
    return pd_path

def lines2array(line):
    """
    input string format from Qchem, eg ['\n','1.1 2.2 3.3\n']
    return: a numpy array with the float format.
    """
    ar = []
    for i in line:
        if len(i) > 2:
            i = i.replace('\n','')
            ar.append([float(x) for x in i.split()])
    ar = np.asarray(ar)
    return ar

def read_modes(path):
    """
    read a vibration analysis output from Qchem calculations.
    return a (3N,3N) matrix of eigenvections of H_vib L = W L
     with frequency w (3N,) 
    """
    out = ccj.find_output(path)
    with open(out,"r")as f:
        done = False
        start = []
        end = []
        rl = f.readlines()
        for i, line in enumerate(rl):
            if "$molecule" in line:
                start.append(i)
                s = start[0]
                done = True
            if done and ("$end" in line):
                end.append(i)
                e = end[0]
            if "Eigenvectors of Proj. Mass-Weighted Hessian Matrix" in line:
                vec_num = i
            if "Vibrational Frequencies in atomic units" in line:
                w_num = i
        N = e - s -2     # get the shape 
        vec_lines = rl[vec_num+1:w_num]
        w_lines = rl[w_num+1:int(w_num+3*N/6+1)]
    vec = lines2array(vec_lines)
    a,b = vec.shape
    col = int(a/(3*N))
    vec_align = np.zeros((3*N,3*N))
    for i in range(1,col+1):  # to rearrange the array
        vec_align[:,6*(i-1):6*i,] = vec[3*N*(i-1):3*N*i,:]
    w = lines2array(w_lines)
    w = np.ndarray.flatten(w)  
    return w,vec_align
        

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
                        tpa['TPA cross section(dI)'].append("-")
                        tpa['Transition dipole moment'].append(ex_transDip)
                        tpa['dipole moment'].append("-")
                        tpa['Total dipole'].append("-")
                    except IndexError:
                        pass
        pd_tpa = pd.DataFrame(tpa).set_index('Excited state')
        pd_tpa.to_csv(os.path.join(path,'tpa.csv'))
        print("Saved tpa.csv in File!")
        return pd_tpa
def read_tpa_pcm(path,correction=True):
    """
    input: the folder with the qchem output
    read:    ADC(2)-PCM TPA calulation
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
                if all(match in line for match in matches) and correction == True:                #if "Excited state" in line:
                     ex_num = line.split()[2] # ['Excited', 'state', '1'] 
                     ex_energy = es.exc_energies(i+12,rl) # ADC(2)/ptSS-PCM(PTD)
                     ex_osi = es.osc_strength(i+24,rl)
                     ex_cro_sos = float(rl[i+27].split()[-1])/30
                     if "Two-photon absorption cross-section" in rl[i+32]:                    
                         ex_cro_dI = float(rl[i+32].split()[-1])/30
                         ex_cro_dI = round(ex_cro_dI,2)
                         ex_dip = es.dipole_moment(i+38,rl)
                         tot_dipole = float(rl[i+39].split()[-1])
                     else:
                         ex_cro_dI = "-"
                         ex_dip = es.dipole_moment(i+33,rl)
                         tot_dipole = float(rl[i+34].split()[-1])
                     #ex_amp = es.amplitudes(i+28,rl).to_dataframe()
                     ex_transDip = es.transition_dipole(i+25,rl)
                     tpa["Excited state"].append(ex_num)
                     tpa['Excitation energy'].append(ex_energy)
                     tpa['Oscillator strength'].append(round(ex_osi,3))
                     tpa['TPA cross section(sos)'].append(round(ex_cro_sos,2))
                     tpa['TPA cross section(dI)'].append(ex_cro_dI)
                     tpa['Transition dipole moment'].append(ex_transDip)
                     tpa['dipole moment'].append(ex_dip)
                     tpa['Total dipole'].append(tot_dipole)
                     #tpa['Transition amplitudes'].append(ex_amp.set_index('weight'))  
                if all(match in line for match in matches) and correction == False:                #if "Excited state" in line:
                     ex_num = line.split()[2] # ['Excited', 'state', '1'] 
                     ex_energy = es.exc_energies(i+8,rl) # Excitation energy (PCM 0th order):
                     ex_osi = es.osc_strength(i+24,rl)
                     ex_cro_sos = float(rl[i+27].split()[-1])/30
                     if "Two-photon absorption cross-section" in rl[i+32]:                    
                         ex_cro_dI = float(rl[i+32].split()[-1])/30
                         ex_cro_dI = round(ex_cro_dI,2)
                         ex_dip = es.dipole_moment(i+38,rl)
                         tot_dipole = float(rl[i+39].split()[-1])
                     else:
                         ex_cro_dI = "-"
                         ex_dip = es.dipole_moment(i+33,rl)
                         tot_dipole = float(rl[i+34].split()[-1])
                     #ex_amp = es.amplitudes(i+28,rl).to_dataframe()
                     ex_transDip = es.transition_dipole(i+25,rl)
                     tpa["Excited state"].append(ex_num)
                     tpa['Excitation energy'].append(ex_energy)
                     tpa['Oscillator strength'].append(round(ex_osi,3))
                     tpa['TPA cross section(sos)'].append(round(ex_cro_sos,2))
                     tpa['TPA cross section(dI)'].append(ex_cro_dI)
                     tpa['Transition dipole moment'].append(ex_transDip)
                     tpa['dipole moment'].append(ex_dip)
                     tpa['Total dipole'].append(tot_dipole)
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
                        tpa['TPA cross section(dI)'].append("-")
                        tpa['Transition dipole moment'].append(ex_transDip)
                        tpa['dipole moment'].append("-")
                        tpa['Total dipole'].append("-")
                    except IndexError:
                        pass
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
    M_ab = <0|u^a|k><k|u^b|f>/(w_k-w_f/2) + a->b
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
    The SOS expression using fluctuation operator. 
    M_ab =sum_k{ <0|u^a|k><k|u^b|f>/(w_k-w_f/2) + a->b}
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

def get_M_sos_entangled(Cart_turple,pd_tpa,M,T_e,considered_state, tot_state,exclude_state=None):
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
    if not exclude_state:
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
    if exclude_state:
        for k in range(1,tot_state):
            if k not in exclude_state:
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
        plt.title("Final TPA cross section for ES "+ str(considered_state) + ": " + str(round(sigma_sos,2)))
        plt.ylabel("TPA cross section contribution (a.u.)")
        plt.xlabel("Excited state (Three-state model)")
#        plt.tight_layout()
        
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


def print_contributions(pd_tpa,considered_state, tot_states):
    """
    input: dataframe, eg. tpa.csv
    h_bar = 6.6*10**(-16) # unit: ev.s, here h_bar = 1--atomic units
    output: print the sos expression and all cotributions from each intemediate state
            and return a list of dominant intermediate states
    """
    lines = 50*'-'
    domin_is = []
    sigma_r = []
    sigma_is = []
    transDipole = []
    detuning = []
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
        E_1 = pd_tpa.loc[str(considered_state),'Excitation energy']/2
        E_k = pd_tpa.loc[str(j),'Excitation energy']
        E_diff = (E_1 - E_k)
        cart = product(range(3), repeat = 2)
        for item in cart:
             M_os = get_M(item,pd_tpa,M,considered_state,j)
        sigma = get_sigma(M_os)
        print("If only consider the " + str(j) + "th intermediate stata (Three-states model), the final cross section:")
        print(sigma)
        sigma_r.append(sigma)
#        print("\n")
        if sigma >= sigma_sos*0.05:
            transDipole.append(round(sigma*E_diff/27.21,2))
            domin_is.append(str(j))
            sigma_is.append(round(sigma,2))
            detuning.append(round(E_diff,3))
    print("The dominant intermediate states are: \n",domin_is)
    print("With the crorss section (au): \n",sigma_is)
#    print("The transition dipole averaged (au):\n",transDipole)
#    print("The detuning energy (eV):\n",detuning)
    print(lines)
    return domin_is,sigma_is,sigma_r
    #    sigma_sum += sigma
    #    print("Sum of previous states:  " + str(sigma_sum))
    #    print("\n")
        

def trans_dipoleD(pd_tpa,considered_state,target_state):
    """
    ouput: a 3*3 matrix with D_a,b = <0|r^a|k><k|r^b|f>+ a->b
    """
    M = np.zeros((3,3))
    cart = product(range(3), repeat = 2) # (0,0),(0,1),(0,2)
    # =============================================================================
    # 2.now generate M matrix
    # =============================================================================
#    E_1 = pd_tpa.loc[str(target_state),'Excitation energy']/2
#    E_k = pd_tpa.loc[str(considered_state),'Excitation energy']
#    E_diff = (E_1 - E_k)/27.21
    for item in cart:
        M_os = get_M(item,pd_tpa,M,target_state,considered_state)
#    D = np.multiply(M_os,E_diff)
    D = M_os
    return D
    
def get_weight(pd_tpa,considered_state1,considered_state2,target_state):
    D_s1 = trans_dipoleD(pd_tpa,considered_state1,target_state)
#    print("The 3*3 matrix of D_k for state "+ str(considered_state1),D_s1)
    D_s2 = trans_dipoleD(pd_tpa,considered_state2,target_state)
#    print("The 3*3 matrix of D_k for state "+ str(considered_state2),D_s2)
    sigma_f = 0
    sigma_g = 0 
    sigma_h = 0
    for i in range(3):
        for j in range(3):
            sigma_f += D_s1[i][i]*D_s2[j][j]
            #print(sigma_f)
            sigma_g += D_s1[i][j]*D_s2[i][j]
            #print(sigma_g)
            sigma_h += D_s1[i][j]*D_s2[j][i]
    #print("sigma_F:",sigma_f)
    #print("sigma_G:",sigma_g)
    #print("sigma_H:",sigma_h)
    weight = (1/30)*(2*sigma_f + 2*sigma_g + 2*sigma_h)
    print("The weight between state "+ str(considered_state1)\
          + " <-> State " + str(considered_state2)+" is: \n",weight)
  

def etpa_time(T_e,pd_tpa,considered_state,tot_states,polarization = "parallel",exclude_state=None):
    """
    input: list of entanglement time in the unit of fs
            exclude_state: a list of states do not take in account
    ouput: list- etpa value
    """
    M_matrix = {}
    sigma_entangled = []
    for t in T_e:    
        #print(t)
        cart = product(range(3), repeat = 2)
        M = np.zeros((3,3),dtype = complex)
        for item in cart:
            get_M_sos_entangled(item,pd_tpa,M,t,considered_state, tot_states,exclude_state)
            #get_M_entangled(item,pd_tpa,M,t,considered_state, 1)
        M_matrix[t] = M
        sig = get_sigma_entangled(M,polarization)
#        sig = sig/(t)
        #print(sig)
        sigma_entangled.append(sig)
    sigma_entangled = np.around(np.real(sigma_entangled),decimals=2)
    print('finished the last one:',sigma_entangled[-1])
    return sigma_entangled

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
    ax.set_xticks(np.arange(data.shape[1]))#, labels=col_labels)
    ax.set_xticklabels(col_labels)
    ax.set_xlabel(r"$|k\rangle$")
    ax.set_yticks(np.arange(data.shape[0]))#, labels=row_labels)
    ax.set_yticklabels(row_labels)
    ax.set_ylabel(r"$|f\rangle$")
    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")
    ax.grid(visible=None)

    # Turn spines off and create white grid.
    for key,spine in ax.spines.items():
        spine.set_visible(False)
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
    plt.show()
    #annotate_heatmap(im, valfmt="{x:d}", size=7, threshold=20,textcolors=("red", "white"))

def exc_ladder(pd_iso,pd_sup,tot_states,s_iso,s_sup,title=" ",yval="Excitation energy",yticks = None,dominant=False):
#    fig, ax = plt.subplots() 
#    plt.figure(figsize=(3,4.5))
    i= 0
    x = [1,2]
    domi_iso,sigma_iso_is,sigma_iso = print_contributions(pd_iso,s_iso, tot_states)
    domi_sup,sigma_sup_is,sigma_sup = print_contributions(pd_sup,s_sup, tot_states)
    if dominant:
        z = [sigma_iso_is,sigma_sup_is]
        y = [pd_iso.loc[domi_iso,yval], pd_sup.loc[domi_sup,yval]]
        print(y)
    else:
        y = [pd_iso.loc["1":str(tot_states),yval],pd_sup.loc["1":str(tot_states),yval]]
        z = [sigma_iso,sigma_sup]
    for xe, ye in zip(x, y):
        plt.scatter([xe] * len(ye), ye,c=z[i],cmap="copper_r",s=2000,marker='_',linewidth=1.5)
        if dominant:
            if i == 0 :
                for j,s in enumerate(domi_iso):
                    plt.annotate("S"+s,xy=(xe+0.1,ye[j]),ha='center',size=12,weight="bold")
            if i == 1 :
                for j,s in enumerate(domi_sup):
                    plt.annotate("S"+s,xy=(xe-0.1,ye[j]),ha='center',size=12,weight="bold") 
        if not dominant:
            if i == 0:
                for j in range(1,tot_states+1):
                    plt.annotate(str(j),xy=(xe+0.2,ye[j-1]),ha='center',size=10,weight="bold")
            if i == 1 :
                for j in range(1,tot_states+1):
                    plt.annotate(str(j),xy=(xe-0.2,ye[j-1]),ha='center',size=10,weight="bold")
        i = i + 1 
    if dominant:
        y_min = min(min(y[0]),min(y[1]))
        y_max = max(max(y[0]),max(y[1]))
        center = (y_max+y_min)/2
        delta = determine_delta(pd_iso,pd_sup,domi_iso,sigma_iso_is,domi_sup,sigma_sup_is,s_iso,s_sup)
        plt.text(1.5,center, r"$\frac{\Delta}{\Delta'}=$"+str(delta),ha='center', va='center',size=14)
        cbar = plt.colorbar(ticks=[min(sigma_sup_is),max(sigma_sup_is)])
        cbar.ax.set_yticklabels(['Low','High'])
    iso_label = r"$|f=$"+str(s_iso)+r"$\rangle_{free}$"
    sup_label = r"$|f'=$"+str(s_sup)+r"$\rangle_{complexed}$"
    plt.xticks([1,2],[iso_label, sup_label],fontsize=11)
    try:
        if yticks.any():
            plt.yticks(yticks)
    except AttributeError:
        pass
    plt.ylabel("Excitation energy (eV)")
#    plt.grid(None)
#    plt.grid(axis='y',linewidth = 0.5)
    plt.grid(which='both',linestyle = '--', axis='y',linewidth = 1)
    plt.title(title)
    plt.show()    
    # the color bar with dual labels
#    cbar = plt.colorbar(pad=0.2)
#    cbar.ax.yaxis.set_ticks_position('left')
#    pos = cbar.ax.get_position()
#    cbar.ax.set_aspect('auto')
#    ax2 = cbar.ax.twinx()
#    ax2.set_ylim(min(sigma_iso),max(sigma_iso))
##    pos.x0 +=0.05
#    cbar.ax.set_position(pos)
#    ax2.set_position(pos)

#    plt.tight_layout()
#    plt.show()    
    
def determine_delta(pd_iso,pd_sup,domi_iso,sigma_iso_is,domi_sup,sigma_sup_is,s_iso,s_sup):
    sort_iso = list(set(sigma_iso_is))
    sort_iso.sort(key=float)
    if len(sort_iso) == 1:
        s_a = domi_iso[-1]
        delta_iso = abs(pd_iso.loc[s_a,"Excitation energy"] - 0.5*pd_iso.loc[str(s_iso),"Excitation energy"])
        print("Only one dominant state from isolated (w_{s_a} - 1/2W_{s_iso}) \
              \n {delta_iso:.3f} : state number S{s_a} ".format(delta_iso =delta_iso,s_a =s_a,s_iso=s_iso))      
    if len(sort_iso) > 1:
        a = sigma_iso_is.index(sort_iso[-1]) # the larges cross section index
        b = sigma_iso_is.index(sort_iso[-2]) # the second largest 
        s_a = domi_iso[a]
        s_b = domi_iso[b]
        delta_iso = abs(pd_iso.loc[s_a,"Excitation energy"] - pd_iso.loc[s_b,"Excitation energy"])
        print("The energy gap between two dominant state from isolated \n {:.3f} : state number S{}, S{} ".format(delta_iso,s_a,s_b))
    sort_sup = list(set(sigma_sup_is))
    sort_sup.sort(key=float)
    if len(sort_sup) == 1:
        s_a = domi_sup[-1]
        delta_sup = abs(pd_sup.loc[s_a,"Excitation energy"] - 0.5*pd_sup.loc[str(s_sup),"Excitation energy"])
        print("Only one dominant state from isolated (w_{s_a} - 1/2W_{s_sup}) \
              \n {delta_sup:.3f} : state number S{s_a} ".format(delta_sup=delta_sup,s_a=s_a,s_sup=s_sup)) 
    if len(sort_sup) > 1:
        a = sigma_sup_is.index(sort_sup[-1]) # the larges cross section index
        b = sigma_sup_is.index(sort_sup[-2]) # the second largest 
        s_a = domi_sup[a]
        s_b = domi_sup[b]
        delta_sup = abs(pd_sup.loc[s_a,"Excitation energy"] - pd_sup.loc[s_b,"Excitation energy"])
        print("The energy gap between two dominant state from complexed \n {:.3f}: state number S{}, S{} ".format(delta_sup,s_a,s_b))
    delta = round(delta_iso/delta_sup,2)
    return delta