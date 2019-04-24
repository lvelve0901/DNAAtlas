
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


step_dic = {'CC':'CG/CG','GG':'GC/GC','CA':'CG/AT','TG':'TA/GC',
            'AA':'AT/AT','TT':'TA/TA','AG':'AT/GC','CT':'CG/TA',
            'GC':'GC/CG','TC':'TA/CG','GA':'GC/AT','CG':'CG/GC',
            'AT':'AT/TA','TA':'TA/AT','GT':'GC/TA','AC':'AT/CG',
           }

seq = 'GCATCGATTGGC'
seq_step = [seq[i:i+2] for i in range(len(seq)-1)]
seq_step = [step_dic[step] for step in seq_step]

print(seq_step)


os.system("fiber -b -seq=%s init_coor.pdb"%seq)
os.system("reduce -DB reduce_het_dict.txt -noadj -Quiet init_coor.pdb > init_coor_h.pdb")

os.system('x3dna-dssr -i=init_coor_h.pdb --more > /dev/null')
os.system('cp dssr-stems.pdb init_coor_h.pdb')
os.system('x3dna-dssr --cleanup')

mol = Mol('./init_coor_h.pdb')

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
newmol.write('./init_coor_rot.pdb')

os.system('x3dna-dssr -i=init_coor_rot.pdb -o=init_coor_rot.json --more --json')
os.system('x3dna-dssr --cleanup')

json_f = './init_coor_rot.json'
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

O = np.array(O)
M = np.array(M)


lib = Pdb('../resources/DNA_step_cluster.pdb')
lib_frame = pd.read_csv("../resources/DNA_step_cluster_frame.csv")
bpstep = deepcopy(lib.mds[13])


# original frame from 2 to 4 

o_01 = O[7]
o_02 = O[8]
o_03 = O[9]

m_01 = M[7]
m_02 = M[8]
m_03 = M[9]

# CA step

o_11 = np.array([0.0,0.0,0.0])
m_11 = np.identity(3)

o_12 = np.array([float(num) for num in lib_frame.ix[13]['origin'].split(',')])
m_12 = np.array([float(num) for num in lib_frame.ix[13]['M_axis'].split(',')]).reshape((3,3))


# rotate CA step to the reference helix

rot, rmsd = lsqfit(m_11,m_01)
d0 = o_11 - o_01

o_12_new = np.dot(o_12 - d0,rot)
m_12_new = np.dot(m_12,rot)

print("old o_02")
print(o_02)
print("new o_02")
print(o_12_new)

m = getats(bpstep)
M = getmat(m)
M -= d0
M = np.dot(M,rot)
putmat(m,M)

# rotate the upper helix part relative to the CA step

rot, rmsd = lsqfit(m_02,m_12_new)
d0 = o_02 - o_12_new

O_upper = np.dot(O[9:12] - d0,rot)
M_upper = np.dot(M[9:12],rot)

print("old o_03")
print(O[9])
print("new o_03")
print(O_upper[0])

newmol_upper = Mol()
newmol_upper.segs = [Segment(),Segment()]
newmol_upper.segs[0].reses = newmol.segs[0].reses[9:12]
newmol_upper.segs[1].reses = newmol.segs[1].reses[0:3]

m = getats(newmol_upper)
M = getmat(m)
M -= d0
M = np.dot(M,rot)
putmat(m,M)


# final output the structure

finalmol = Mol()
finalmol.segs = [Segment(),Segment()]
finalmol.segs[0].reses = newmol.segs[0].reses[0:7] + bpstep.segs[0].reses + newmol_upper.segs[0].reses
finalmol.segs[1].reses = newmol_upper.segs[1].reses + bpstep.segs[1].reses + newmol.segs[1].reses[5:12]
finalmol.renumber(1,1)

finalmol.write('test_final.pdb')

