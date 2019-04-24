

#!/usr/bin/python

import os
import sys
import re
import json
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

def gen_init_seed(seq,name):

    helix = Helix(seq,name)

    os.system("fiber -b -seq=%s init_%s.pdb"%(seq,name))
    os.system("reduce -DB reduce_het_dict.txt -noadj -Quiet init_%s.pdb > init_%s_h.pdb"%(name,name))

    os.system('x3dna-dssr -i=init_%s_h.pdb --more > /dev/null'%name)
    os.system('cp dssr-stems.pdb init_%s_h.pdb'%name)
    os.system('x3dna-dssr --cleanup')

    mol = Mol('./init_%s_h.pdb'%name)

    o1 = np.array([0.0,0.0,0.0])
    A1 = np.array([[0.0,0.0,-1.0],[-1.0,0.0,0.0],[0.0,1.0,0.0]])
    o0 = np.array([0.0,0.0,0.0])
    A0 = np.identity(3)

    rot, rmsd = lsqfit(A1,A0)
    d0 = o1 - o0

    m = getats(mol)
    M = getmat(m)
    M -= d0
    M = np.dot(M,rot)
    putmat(m,M)

    newmol = Mol()
    newmol.segs = [Segment(),Segment()]
    newmol.segs[0].reses = mol.segs[0].reses[0:len(seq)]
    newmol.segs[1].reses = mol.segs[0].reses[len(seq):]
    newmol.write('./init_%s_rot.pdb'%name)
   
    os.system("rm init_%s.pdb init_%s_h.pdb"%(name,name))
    os.system('x3dna-dssr -i=init_%s_rot.pdb -o=init_%s_rot.json --more --json'%(name,name))
    os.system('x3dna-dssr --cleanup')
    
    json_f = './init_%s_rot.json'%name
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
        data = json.load(json_data)
    
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file
    
    O = []
    M = []
    
    bps = najson.json_file['pairs']
    for bp in bps:
        O.append(bp['frame']['origin'])
        M.append([bp['frame']['x_axis'],bp['frame']['y_axis'],bp['frame']['z_axis']])

    helix.origins = np.array(O)
    helix.M_axes  = np.array(M)
    
    os.system("rm ./init_%s_rot.json"%name)

    return helix

