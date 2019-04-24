
import json # Handle JSON DSSR files
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import learnna_json as lna_json
from pdblib.num import *


#json_f1 = './A2dna.json'
#najson1 = lna_json.NA_JSON()  #initialize class objects
#with open(json_f1) as json_data:  #read each json file
#    data1 = json.load(json_data)
#
#najson1.set_json(data1)  #pass json file to class pbject
#najson1.read_idx()  #set index from own json file
#
#json_f2 = './revA2dna.json'
#najson2 = lna_json.NA_JSON()  #initialize class objects
#with open(json_f2) as json_data:  #read each json file
#    data2 = json.load(json_data)
#
#najson2.set_json(data2)  #pass json file to class pbject
#najson2.read_idx()  #set index from own json file
#
#bps1 =  najson1.json_file['pairs']
#bps2 =  najson2.json_file['pairs']
#
#print bps1[0]['frame']
#print bps1[1]['frame']
#print bps2[10]['frame']
#print bps2[11]['frame']




#json_f3 = './3kz8_0_4_new.json'
#najson3 = lna_json.NA_JSON()  #initialize class objects
#with open(json_f3) as json_data:  #read each json file
#    data3 = json.load(json_data)
#
#najson3.set_json(data3)  #pass json file to class pbject
#najson3.read_idx()  #set index from own json file
#
#bps3 =  najson3.json_file['pairs']
#print bps3[0]['nt1']
#print bps3[0]['frame']
#print bps3[1]['nt1']
#print bps3[1]['frame']
#
#json_f4 = './3kz8_0_4_new_.json'
#najson4 = lna_json.NA_JSON()  #initialize class objects
#with open(json_f4) as json_data:  #read each json file
#    data4 = json.load(json_data)
#
#najson4.set_json(data4)  #pass json file to class pbject
#najson4.read_idx()  #set index from own json file
#
#bps4 =  najson4.json_file['pairs']
#print bps4[0]['nt1']
#print bps4[0]['frame']
#print bps4[1]['nt1']
#print bps4[1]['frame']


#json_f5 = './4i9v_0_12.json'
#najson5 = lna_json.NA_JSON()  #initialize class objects
#with open(json_f5) as json_data:  #read each json file
#    data5 = json.load(json_data)
#
#najson5.set_json(data5)  #pass json file to class pbject
#najson5.read_idx()  #set index from own json file
#
#bps5 =  najson5.json_file['pairs']
##print bps5[0]['nt1'] + "&" + bps5[0]['nt2']
##print bps5[0]['frame']
#print bps5[4]['nt1'] + "&" + bps5[4]['nt2']
#print bps5[4]['bp']
#print bps5[4]['frame']
#print bps5[5]['nt1'] + "&" + bps5[5]['nt2']
#print bps5[5]['bp']
#print bps5[5]['frame']
#print bps5[6]['nt1'] + "&" + bps5[6]['nt2']
#print bps5[6]['bp']
#print bps5[6]['frame']
##print bps5[3]['nt1'] + "&" + bps5[3]['nt2']
##print bps5[3]['frame']
#stem5_pair = najson5.json_file['stems'][0]['pairs']
#for pair in stem5_pair:
#    print((pair['nt1'],pair['nt2']))

json_f5 = './test_step.json'
najson5 = lna_json.NA_JSON()  #initialize class objects
with open(json_f5) as json_data:  #read each json file
    data5 = json.load(json_data)

najson5.set_json(data5)  #pass json file to class pbject
najson5.read_idx()  #set index from own json file

bps5 =  najson5.json_file['pairs']
print bps5[0]['nt1'] + "&" + bps5[0]['nt2']
print bps5[0]['frame']
print bps5[1]['nt1'] + "&" + bps5[1]['nt2']
print bps5[1]['frame']

