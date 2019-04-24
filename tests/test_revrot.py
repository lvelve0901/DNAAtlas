
import json # Handle JSON DSSR files
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import learnna_json as lna_json
from common.matvec import lsqfit
from pdblib.num import *


json_f1 = './A2dna.json'
najson1 = lna_json.NA_JSON()  #initialize class objects
with open(json_f1) as json_data:  #read each json file
    data1 = json.load(json_data)

najson1.set_json(data1)  #pass json file to class pbject
najson1.read_idx()  #set index from own json file

json_f2 = './revA2dna.json'
najson2 = lna_json.NA_JSON()  #initialize class objects
with open(json_f2) as json_data:  #read each json file
    data2 = json.load(json_data)

najson2.set_json(data2)  #pass json file to class pbject
najson2.read_idx()  #set index from own json file

bps1 =  najson1.json_file['pairs']
bps2 =  najson2.json_file['pairs']

print bps1[0]['nt1'] + "&" + bps1[0]['nt2']
print bps1[0]['frame']
print bps1[1]['nt1'] + "&" + bps1[1]['nt2']
print bps1[1]['frame']
print bps2[10]['nt1'] + "&" + bps2[10]['nt2']
print bps2[10]['frame']
print bps2[11]['nt1'] + "&" + bps2[11]['nt2']
print bps2[11]['frame']

o_ref1 = np.array(bps1[0]['frame']['origin'])
o_ref2 = np.array(bps1[1]['frame']['origin'])
M_ref1 = np.array([bps1[0]['frame']['x_axis'],bps1[0]['frame']['y_axis'],bps1[0]['frame']['z_axis']])
M_ref2 = np.array([bps1[1]['frame']['x_axis'],bps1[1]['frame']['y_axis'],bps1[1]['frame']['z_axis']])


o_ref1_ = np.array(bps2[10]['frame']['origin'])
o_ref2_ = np.array(bps2[11]['frame']['origin'])
M_ref1_ = np.array([bps2[10]['frame']['x_axis'],bps2[10]['frame']['y_axis'],bps2[10]['frame']['z_axis']])
M_ref2_ = np.array([bps2[11]['frame']['x_axis'],bps2[11]['frame']['y_axis'],bps2[11]['frame']['z_axis']])

print o_ref1
print o_ref2
print M_ref1
print M_ref2

print o_ref1_
print o_ref2_
print M_ref1_
print M_ref2_

o0 = np.array([0,0,0])
A0 = np.identity(3)

rot, rmsd = lsqfit(M_ref1,A0)
d = o_ref2 - o_ref1
B1 = np.dot(M_ref1,rot)
B2 = np.dot(M_ref2,rot)
dd = np.dot(d,rot)

print dd
print np.sqrt(np.sum((dd)**2))
print B2

rot, rmsd = lsqfit(M_ref1_,A0)
d = o_ref2_ - o_ref1_
B1 = np.dot(M_ref1_,rot)
B2 = np.dot(M_ref2_,rot)
dd = np.dot(d,rot)

print dd
print np.sqrt(np.sum((dd)**2))
print B2

