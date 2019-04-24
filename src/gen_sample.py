

#!/usr/bin/python

import os
import sys
import re
import json
import random
import pandas as pd
import numpy as np
import learnna_json as lna_json
from pdblib.num import *
from common.matvec import lsqfit
from copy import deepcopy
from helix_base import *


step_dic = {'CC':'CG/CG','GG':'GC/GC','CA':'CG/AT','TG':'TA/GC',
            'AA':'AT/AT','TT':'TA/TA','AG':'AT/GC','CT':'CG/TA',
            'GC':'GC/CG','TC':'TA/CG','GA':'GC/AT','CG':'CG/GC',
            'AT':'AT/TA','TA':'TA/AT','GT':'GC/TA','AC':'AT/CG',
           }

def gen_sample(helix,mode="MC",T=1.):

    if mode not in ["MC","random"]:
        print("ERROR(): Does not support mode: %s"%mode)
        print("mode can be 'MC' or 'random'")
        sys.exit()

    seq = deepcopy(helix.seq)
    seqlen = deepcopy(helix.seqlen)

    # randomly select an bp step to swap
    randseed = random.randint(0,seqlen-2)
    #randseed = 1

    cur_step = deepcopy(helix.seq_step[randseed])
    cur_step_idx = deepcopy(helix.step_idx[randseed])
    cur_mdid_idx = deepcopy(helix.mdid_idx[randseed])
    cur_energy = deepcopy(helix.energys[randseed])

    # randomly select an bp step from library given the same sequence
    step_idxes = cluster_step_frame.loc[cluster_step_frame['step_name'] == cur_step].index
    new_step_idx = random.choice(step_idxes)
    new_mdid_idx = cluster_step_frame['step_idx'].ix[new_step_idx] 
    new_energy = -8.314*(25+273.15)*np.log(cluster_step_frame['pop'].ix[new_step_idx])*0.000239

    if mode == "MC":
        if new_energy > cur_energy:
            P_acc = np.exp((cur_energy - new_energy)/T)
            coin = random.uniform(0,1)
            if P_acc <= coin:
                return helix

    # old origin and M_axis
    o_01 = deepcopy(helix.origins[randseed])
    m_01 = deepcopy(helix.M_axes[randseed])
    o_02 = deepcopy(helix.origins[randseed+1])
    m_02 = deepcopy(helix.M_axes[randseed+1])
    
    # new origin and M_axis
    o_11 = np.array([0.0,0.0,0.0])
    m_11 = np.identity(3)
    o_12 = np.array([float(num) for num in cluster_step_frame.ix[new_step_idx]['origin'].split(',')])
    m_12 = np.array([float(num) for num in cluster_step_frame.ix[new_step_idx]['M_axis'].split(',')]).reshape((3,3))

    # rotate and anchor the bp step
    rot1, rmsd1 = lsqfit(m_11,m_01)
    d1 = o_11 - o_01
    o_12_new = np.dot(o_12 - d1,rot1)
    m_12_new = np.dot(m_12,rot1)
    helix.origins[randseed+1] = o_12_new
    helix.M_axes[randseed+1] = m_12_new
    helix.step_idx[randseed] = new_step_idx
    helix.mdid_idx[randseed] = new_mdid_idx
    helix.energys[randseed] = new_energy


    # rotate the entire upper helix
    rot2, rmsd2 = lsqfit(m_02,m_12_new)
    d2 = o_02 - o_12_new

    #if randseed+2 < seqlen:
    #    helix.origins[randseed+2:] = np.dot(deepcopy(helix.origins[randseed+2:]) - d2,rot2)
    #    helix.M_axes[randseed+2:] = np.dot(deepcopy(helix.M_axes[randseed+2:]),rot2)
    helix.origins[randseed+2:] = np.dot(deepcopy(helix.origins[randseed+2:]) - d2,rot2)
    helix.M_axes[randseed+2:] = np.dot(deepcopy(helix.M_axes[randseed+2:]),rot2)

    print(randseed)
    print(d1)
    print(d2)
    print(helix.step_idx)
    print(helix.origins)

    return helix 


