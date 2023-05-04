#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:42:18 2022

@author: fu
"""

import numpy as np
import tools
from itertools import product
import os
root = os.getcwd()
#folder = os.path.join(root, "chalcone")#,"dis_bigger")
folder = os.path.join(root, "vibration_coupling")#,"dis_bigger")
folder_a = os.path.join(folder,"dis_0.08","dis_a")
folder_s = os.path.join(folder,"dis_0.08","dis_s")
#folder_a = os.path.join(folder,'dis_0.08',"dis_a")
#folder_s = os.path.join(folder,'dis_0.08',"dis_s")
np_eigen = np.loadtxt(os.path.join(folder,"eigen.csv"),delimiter=",")
lamda = np_eigen[0,:] # vibration frequency in a.u.
L = np_eigen[1:,:]
#lamda, L = tools.read_modes(folder)
tot_states = 7
considered_state = 2
displacement =  0.08
shape = len(lamda) # parse the number 3N
coor_a = {}
coor_s = {}
M_a = np.zeros((9,shape)) 
M_s = np.zeros((9,shape))
for i in range(1,shape+1):
    coor_a[i] = tools.parse_tpa(os.path.join(folder_a,"coord_"+str(i)))
    coor_s[i] = tools.parse_tpa(os.path.join(folder_s,"coord_"+str(i)))    
 #get S tpa transion moment in displacement a + and s subtracting    
    M = np.zeros((3,3))
    cart = product(range(3), repeat = 2)
    for item in cart:
       M_dummy = tools.get_M_sos(item, coor_a[i], M, considered_state,tot_states)
    M_a[:,i-1] = M_dummy.reshape((-1))  
    cart = product(range(3), repeat = 2)
    M = np.zeros((3,3))
    for item in cart:
        M_dummy2= tools.get_M_sos(item, coor_s[i], M, considered_state,tot_states)
    M_s[:,i-1] = M_dummy2.reshape((-1))
  # using finite differential (M_a-M_s)/2d respect to cartian displacement
partial_C = (M_a-M_s)/(2*displacement)
# transform to normal coordinate basis
partial_Q = np.einsum('ji,lj,ij->lj',L,partial_C,L)
#partial_Q = np.dot(partial_C, L)
#partial_Q = np.dot(partial_Q, np.transpose(L))
#partial_Q = np.einsum('jkl,ij->ikl',partial_Q, np.transpose(L)) # this gives back to cartian basis
 # sum of all normal nodes sum_v sigma_v*hbar/2w_v
sigma_vib_sum = 0
for q in range(len(lamda)):
    if lamda[q] != 0: 
        sigma_v = tools.get_sigma(partial_Q[:,q].reshape(3,3))
        sigma_v /= (2*lamda[q])
        print("vibration frequency:",lamda[q])
        print("sigma_v for the mode "+str(q)+" is :",sigma_v)
        sigma_vib_sum += sigma_v
print("vibronic contribution for the electronic state "+str(considered_state)+"-->",sigma_vib_sum)