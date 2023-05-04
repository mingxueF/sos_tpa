#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:22:01 2022

@author: fu
"""

import numpy as np
import tools
from itertools import product
import os
root = os.getcwd()
folder = os.path.join(root, "vibration_coupling")
np_eigen = np.loadtxt(os.path.join(folder,"eigen.csv"),delimiter=",")
out_a = os.path.join(folder,"dis_a_all")
out_s = os.path.join(folder,"dis_s_all")
lamda = np_eigen[0,:]
L = np_eigen[1:,:]
tpa_a = tools.parse_tpa(out_a)
tpa_s = tools.parse_tpa(out_s)
# get S tpa transion moment in displacement a + and s subtracting
tot_states = 10
#for state in range(1,tot_states):
considered_state = 2
M = np.zeros((3,3))
cart = product(range(3), repeat = 2)
for item in cart:
   M_sos_a = tools.get_M_sos(item, tpa_a, M, considered_state,tot_states)
cart = product(range(3), repeat = 2)
M = np.zeros((3,3))
for item in cart:
   M_sos_s = tools.get_M_sos(item, tpa_s, M, considered_state,tot_states)
# using finite differential (M_a-M_s)/2d respect to cartian displacment
partial_s = (M_sos_a-M_sos_s)/0.018
# reshape to 3N,
partial = np.tile(partial_s,(18,1))
partial = np.reshape(partial,(18,3,3))
# transform to normal mode coordinate basis
partial_Q = np.einsum('ij,jkl->ikl',np.transpose(L),partial)
#partial_c = np.einsum('ij,jkl->ikl',np.transpose(L),partial_Q)
#partial_shape = np.reshape(partial_Q,(6,3,3))
# get the sigma 
#sigma_v = tools.get_sigma(partial_shape[1,:,:])
# sum of all normal nodes sum_v sigma_v*hbar/2w_v
sigma_vib_sum = 0
for n in range(len(lamda)):
    if lamda[n] != 0: 
        sigma_v = tools.get_sigma(partial_Q[n,:,:])
        print(sigma_v)
        sigma_vib_sum += sigma_v/(2*lamda[n])
print("vibronic contribution for state "+str(considered_state),sigma_vib_sum)