
"""
helix_base, version 1.0, 04-19-2017
Written by Honglue Shi
Copyright: Honglue Shi and Al-hashimi's group @Duke University

"""
#############################################
#                 Imports                   #
#############################################

# Standard library import
import os
import re
import pandas as pd
import numpy as np
from pdblib.num import *
from sys import stdout,exit
from operator import add
from common.matvec import lsqfit
from common.base import divide

#############################################
#            Global Variables               #
#############################################

# uniform backbone atom library

step_dic = {'CC':'CG/CG','GG':'GC/GC','CA':'CG/AT','TG':'TA/GC',
            'AA':'AT/AT','TT':'TA/TA','AG':'AT/GC','CT':'CG/TA',
            'GC':'GC/CG','TC':'TA/CG','GA':'GC/AT','CG':'CG/GC',
            'AT':'AT/TA','TA':'TA/AT','GT':'GC/TA','AC':'AT/CG',
           }

cluster_step_frame = pd.read_csv("../resources/DNA_step_cluster_frame.csv")
cluster_step_pdb  = Pdb("../resources/DNA_step_cluster.pdb")


#############################################
#                Classes                    #
#############################################

# Internal use classes

#====== Continue ============================================
class Continue(Exception):
    pass

# Public use classes

#====== Atom =================================================
class Helix:

    seq = ""
    seqlen = 0
    name = ""
    seq_step = []
    step_idx = np.array([])
    mdid_idx = np.array([])
    energys = np.array([])
    origins = np.array([])
    M_axes = np.array([])
    isSwap = False

    #---------------------------------------------------------
    def __init__(self,seq,name):
        
        self.seq = seq
        self.seqlen = len(seq)
        self.name = name
        self.seq_step = [step_dic[step] for step in [seq[i:i+2] for i in range(len(seq)-1)]]
        self.step_idx = np.array([-1]*(len(seq)-1))
        self.mdid_idx = np.array([-1]*(len(seq)-1))
        self.energys = np.array([np.inf]*(len(seq)-1))
        self.isSwap = False

    #---------------------------------------------------------
    def show(self):
        
        output = []
        output.append('Starting printing helix class:\n')
        output.append(('seq: ' + '\n' + '%s' + '\n')%self.seq)
        output.append(('seqlen: ' + '\n' + '%s' + '\n')%self.seqlen)
        output.append(('seq_step: ' + '\n' + '%s' + '\n')%self.seq_step)
        output.append(('step_idx: ' + '\n' + '%s'+'\n')%self.step_idx)
        output.append(('mdid_idx: ' + '\n' + '%s'+'\n')%self.mdid_idx)
        output.append(('energys: ' + '\n' + '%s'+'\n')%self.energys)
        output.append(('isSwap: ' + '\n' + '%s'+'\n')%self.isSwap)
        #output.append(('\n' + 'origins: ' + '\n'))
        #for origin in self.origins:
        #    output.append(('%s'+'\n')%origin)
        #output.append(('\n' + 'M_axes: ' + '\n'))
        #for M_axis in self.M_axes:
        #    output.append(('%s'+'\n')%M_axis)
        output.append('Ending printing helix class\n')
        output.append('\n')
        return stdout.writelines(output)
    
    def write(self,fname):

        omol = Mol()
        omol.segs = [Segment(),Segment()]

        for i, idx in enumerate(self.step_idx[::2]):

            o_01 = self.origins[2*i]
            m_01 = self.M_axes[2*i]
            o_02 = self.origins[2*i+1]
            m_02 = self.M_axes[2*i+1]

            bpstep = deepcopy(cluster_step_pdb.mds[idx])

            o_11 = np.array([0.0,0.0,0.0])
            m_11 = np.identity(3)
            o_12 = np.array([float(num) for num in cluster_step_frame.ix[idx]['origin'].split(',')])
            m_12 = np.array([float(num) for num in cluster_step_frame.ix[idx]['M_axis'].split(',')]).reshape((3,3))

            rot, rmsd = lsqfit(m_11,m_01)
            d0 = o_11 - o_01

            m = getats(bpstep)
            M = getmat(m)
            M -= d0
            M = np.dot(M,rot)
            putmat(m,M)

            omol.segs[0].reses += bpstep.segs[0].reses
            omol.segs[1].reses += bpstep.segs[1].reses[::-1]
        
        omol.segs[1].reses = omol.segs[1].reses[::-1]
        omol.renumber(1,1)
        omol.write(fname)


#############################################
#             Module methods                #
#############################################


