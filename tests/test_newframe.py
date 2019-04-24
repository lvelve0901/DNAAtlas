
import json # Handle JSON DSSR files
import os
import re
import sys
import pandas as pd
import numpy as np
import learnna_json as lna_json
from pdblib.num import *
from common.matvec import lsqfit
from copy import copy, deepcopy


json_f = './4i9v_0_12.json'
najson = lna_json.NA_JSON()  #initialize class objects
with open(json_f) as json_data:  #read each json file
    data = json.load(json_data)

najson.set_json(data)  #pass json file to class pbject
najson.read_idx()  #set index from own json file


stem_len = 12
mol = Mol('./4i9v_0_12.pdb')
oo = np.array([0,0,0])
MM = np.identity(3)
bps =  najson.json_file['pairs']
idx = 4
frame0 = bps[idx]['frame']
frame1 = bps[idx+1]['frame']

print(" === forward === ")

origin0 = np.array(frame0['origin'])
origin1 = np.array(frame1['origin'])
M_axis0 = np.array([np.array(frame0['x_axis']),np.array(frame0['y_axis']),np.array(frame0['z_axis'])])
M_axis1 = np.array([np.array(frame1['x_axis']),np.array(frame1['y_axis']),np.array(frame1['z_axis'])])
d = origin0 - oo
dd = origin1 - origin0
rot, rmsd = lsqfit(M_axis0,MM)
M_axis00 = np.dot(M_axis0,rot)
M_axis11 = np.dot(M_axis1,rot)
ddd = np.dot(dd,rot)

print("origin0")
print(origin0)
print("origin1")
print(origin1)
print("d")
print(d)
print("rot")
print(rot)
print("rmsd")
print(rmsd)
print("M_axis0")
print(M_axis0)
print("M_axis1")
print(M_axis1)
print("dd")
print(dd)
print("M_axis00")
print(M_axis00)
print("M_axis11")
print(M_axis11)
print("ddd")
print(ddd)

newmol = Mol()
newmol.segs = [Segment(),Segment()]
mol.segs[0].reses[idx].show()
mol.segs[0].reses[idx+1].show()
mol.segs[1].reses[stem_len-2-idx].show()
mol.segs[1].reses[stem_len-2+1-idx].show()
newmol.segs[0].reses = [deepcopy(mol.segs[0].reses[idx]),deepcopy(mol.segs[0].reses[idx+1])]
newmol.segs[1].reses = [deepcopy(mol.segs[1].reses[stem_len-2-idx]),deepcopy(mol.segs[1].reses[stem_len-2+1-idx])]
#newmol.renumber(1,1)

m = getats(newmol)
M = getmat(m)
M -= d
M = np.dot(M,rot)
putmat(m,M)
newmol.mdid = idx+1
newmol.write('./4i9v_0_12_new.pdb')

print(" === reverse === ")

origin0_ = np.array(frame1['origin'])
origin1_ = np.array(frame0['origin'])
M_axis0_ = np.array([np.array(frame1['x_axis']),np.array(frame1['y_axis'])*-1,np.array(frame1['z_axis'])*-1])
M_axis1_ = np.array([np.array(frame0['x_axis']),np.array(frame0['y_axis'])*-1,np.array(frame0['z_axis'])*-1])
d_ = origin0_ - oo
dd_ = origin1_ - origin0_
rot_, rmsd_ = lsqfit(M_axis0_,MM)
M_axis00_ = np.dot(M_axis0_,rot_)
M_axis11_ = np.dot(M_axis1_,rot_)
ddd_ = np.dot(dd_,rot_)

print("origin0_")
print(origin0_)
print("origin1_")
print(origin1_)
print("d_")
print(d_)
print("M_axis0_")
print(M_axis0_)
print("M_axis1_")
print(M_axis1_)
print("dd_")
print(dd_)
print("rot_")
print(rot_)
print("rmsd_")
print(rmsd_)
print("M_axis00_")
print(M_axis00_)
print("M_axis11_")
print(M_axis11_)
print("ddd_")
print(ddd_)

newmol_ = Mol()
newmol_.segs = [Segment(),Segment()]
mol.segs[1].reses[stem_len-2-idx].show()
mol.segs[1].reses[stem_len-2+1-idx].show()
mol.segs[0].reses[idx].show()
mol.segs[0].reses[idx+1].show()
newmol_.segs[0].reses = [deepcopy(mol.segs[1].reses[stem_len-2-idx]),deepcopy(mol.segs[1].reses[stem_len-2+1-idx])]
newmol_.segs[1].reses = [deepcopy(mol.segs[0].reses[idx]),deepcopy(mol.segs[0].reses[idx+1])]
#newmol_.renumber(1,1)

m_ = getats(newmol_)
M_ = getmat(m_)
M_ -= d_
M_ = np.dot(M_,rot_)
putmat(m_,M_)
newmol_.mdid = idx+1
newmol_.write('./4i9v_0_12_new_.pdb')

print(newmol.segs[0].reses[0])
print(newmol_.segs[1].reses[0])

