

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
from supimp import *


step_dic = {'CC':'CG/CG','GG':'GC/GC','CA':'CG/AT','TG':'TA/GC',
            'AA':'AT/AT','TT':'TA/TA','AG':'AT/GC','CT':'CG/TA',
            'GC':'GC/CG','TC':'TA/CG','GA':'GC/AT','CG':'CG/GC',
            'AT':'AT/TA','TA':'TA/AT','GT':'GC/TA','AC':'AT/CG',
           }

def gen_sample_full(helix,mol,mode="MC",T=1.):

    if mode not in ["MC","random"]:
        print("ERROR(): Does not support mode: %s"%mode)
        print("mode can be 'MC' or 'random'")
        sys.exit()

    seq = deepcopy(helix.seq)
    seqlen = deepcopy(helix.seqlen)

    # randomly select an bp step to swap
    randseed = random.randint(0,seqlen-2)

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
                helix.isSwap = False
                return helix, mol

    helix.step_idx[randseed] = new_step_idx
    helix.mdid_idx[randseed] = new_mdid_idx
    helix.energys[randseed] = new_energy
    helix.isSwap = True

    # old bp step
    bp_old = Mol()
    bp_old.segs = [Segment(),Segment()]
    bp_old.segs[0].reses = deepcopy(mol.segs[0].reses[randseed:randseed+2])
    bp_old.segs[1].reses = deepcopy(mol.segs[1].reses[seqlen-randseed-2:seqlen-randseed])

    # new bp step
    bp_new = deepcopy(cluster_step_pdb.mds[new_step_idx])

    # supimpose new bp step to old bp step
    supimpose(bp_new,bp_old,[0,1],[0,1],[0,1],[0,1])

    # upper helix
    upper = Mol()
    upper.segs = [Segment(),Segment()]
    upper.segs[0].reses = mol.segs[0].reses[randseed+1:]
    upper.segs[1].reses = mol.segs[1].reses[0:seqlen-randseed-1]

    # supimpose upper helix to new bp step
    supimpose(upper,bp_new,[0,1],[0,len(upper.segs[1].reses)-1],[0,1],[1,0])
    
    # write new coordinate
    
    newmol = Mol()
    newmol.segs = [Segment(),Segment()]
    newmol.segs[0].reses = mol.segs[0].reses[0:randseed] + bp_new.segs[0].reses + mol.segs[0].reses[randseed+2:seqlen]
    newmol.segs[1].reses = mol.segs[1].reses[0:seqlen-randseed-2] + bp_new.segs[1].reses + mol.segs[1].reses[seqlen-randseed:seqlen]
    
    return helix, newmol 
