
#!/usr/bin/python

import re
import os
import sys
import json
import pandas as pd
import numpy as np
import learnna_json as lna_json
from commontool import read, readchar
from common.matvec import lsqfit
from pdblib.num import *


# load DSSR json file
def load_json(filename):
    
    json_f = filename
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
        data = json.load(json_data)

    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file
    
    return najson.json_file  

# generate a dictionary of nts index
def gen_ntsidx(nts):

    nts_idx = {}  #nts index dictionary with nt_id as key
    for i, nt in enumerate(nts):
        ntid = nt['nt_id']
        if ntid in nts_idx.keys():
            print("ERROR(): Potential same labeling of nucleotides")
            print("%s"%ntid)
            sys.exit()
        nts_idx[ntid] = i

    return nts_idx

# generate a dictionary of bps index
def gen_bpsidx(bps):

    bps_idx = {}
    for i, bp in enumerate(bps):
        bpid  = bp['nt1'] + "&" + bp['nt2']
        bpid_ = bp['nt2'] + "&" + bp['nt1']
        if bpid in bps_idx.keys():
            print("ERROR(): Potential same labeling of base pairs")
            print("%s"%bpid)
            sys.exit()
        if bpid_ in bps_idx.keys():
            print("ERROR(): Potential same labeling of base pairs")
            print("%s"%bpid_)
            sys.exit()
        bps_idx[bpid] = i

    return bps_idx

# find the bp_id and bp_info
def find_bpid(nt1_id,nt2_id,bps,bps_idx):

    try:
        isSwap = False
        bp_id = nt1_id + "&" + nt2_id
        bp = bps[bps_idx[bp_id]]
    except KeyError, e:
        print("WARNING(): Cannot find bp_id: %s"%bp_id)
        print("Try swapping nt1 and nt2")
        try:
            isSwap = True
            bp_id = nt2_id + "&" + nt1_id
            bp = bps[bps_idx[bp_id]]
            bp_id = nt1_id + "&" + nt2_id
            bp_id_ = nt2_id + "&" + nt1_id
        except KeyError, e:
            print("ERROR(): Still cannot find bp_id: %s"%bp_id)
            sys.exit()
    
    origin  = bp['frame']['origin']
    x_axis  = bp['frame']['x_axis']
    if isSwap == False:
        y_axis  = bp['frame']['y_axis']
        z_axis  = bp['frame']['z_axis']
        bp_name = "".join([bp['bp'][0],bp['bp'][-1]])
    else:
        y_axis  = bp['frame']['y_axis'] * -1
        z_axis  = bp['frame']['z_axis'] * -1
        bp_name = "".join([bp['bp'][0],bp['bp'][-1]])[::-1]

    # normalization
    x_axis  = x_axis  / np.linalg.norm(x_axis)
    y_axis  = y_axis  / np.linalg.norm(y_axis)
    z_axis  = z_axis  / np.linalg.norm(z_axis)


    bp_info = [bp_id,bp_name,origin,x_axis,y_axis,z_axis]

    return bp_id, bp_info

# assign the stem helix form
def gen_sform(sforms):

    sfrom = '.'
    if "A" in sforms and "B" not in sforms and "Z" not in sforms: sform = "A"
    elif "A" not in sforms and "B" in sforms and "Z" not in sforms: sform = "B"
    elif "A" not in sforms and "B" not in sforms and "Z" in sforms: sform = "Z"
    elif "A" in sforms and "B" in sforms and "Z" not in sforms: sform = "AB"
    elif "A" in sforms and "B" not in sforms and "Z" in sforms: sform = "AZ"
    elif "A" not in sforms and "B" in sforms and "Z" in sforms: sform = "BZ"
    elif "A" not in sforms and "B" not in sforms and "Z" not in sforms: sform = "."
    else:
        print("ERROR(): Unassigned mixed helix form in the stem: %s"%(sforms))
        sys.exit()
    
    return sform

# rotate the reference frame of base pair 1
def rot_frame(o1,o2,A1,A2):
    
    o0 = np.array([0,0,0])
    A0 = np.identity(3)

    rot, rmsd = lsqfit(A1,A0)
    d = o2 - o1
    B1 = np.dot(A1,rot)
    B2 = np.dot(A2,rot)
    dd = np.dot(d,rot)

    return dd, B2
    
def main():

    step_idx = 0
    jsondir = '/mnt/hs189/NAfinder_3.0/Crystal/Json'
    pdblist = []

    df = pd.read_csv('/mnt/hs189/NAfinder_3.0/doc/NDB/PDB_NDB/Final_library_DNA_Representative.csv')
    pdblist = df.loc[(df.hasLigand == False) & (df.reso < 2.0)]['pdbid'].tolist()
    pdblist.sort()

    #for file in read('/mnt/hs189/NAfinder_3.0/Crystal/doc/all_assemb_id.txt'):
    #    pdblist.append(file[0])
    #for file in read('/mnt/hs189/NAfinder_3.0/Crystal/doc/all_bundle_file.txt'):
    #    pdblist.append(file[0])
    #pdblist.sort()  #sort the pdblist

    #pdblist = ['1bna','3kz8']

    output = [['pdbid','stem_idx','sform','step_idx','step_name','step_id','step_form','origin','M_axis']]

    for i, pdb in enumerate(pdblist):

        pdbid = pdb.split('-')[0]
        print("--- Working on pdbid %s [%d/%d] ---"%(pdbid,i,len(pdblist)))
        
        jsonf = os.path.join(jsondir,pdb+'.json') 
        json_file = load_json(jsonf)
   
        if 'stems' not in json_file.keys():
            continue

        nts = json_file['nts']
        nts_idx = gen_ntsidx(nts)
        bps = json_file['pairs']
        bps_idx = gen_bpsidx(bps)
        stems = json_file['stems']
    
        for j, stem in enumerate(stems):
            
            stem_idx = j
            stem_pair = stem['pairs']
            stem_forms = stem['helix_form']
            sform = gen_sform(stem_forms)

            # ignore all non-Bform DNA
            if sform != 'B':
                continue
    
            bp_pairs = []
            
            for pair in stem_pair:
                bp_pairs.append((pair['nt1'],pair['nt2']))
            
            for idx in range(0,len(bp_pairs)-1):
                bp1_id, bp1_info = find_bpid(bp_pairs[idx][0],  bp_pairs[idx][1],bps,bps_idx  )
                bp2_id, bp2_info = find_bpid(bp_pairs[idx+1][0],bp_pairs[idx+1][1],bps,bps_idx)
                step_id = bp1_id + '&&' + bp2_id
                step_name = bp1_info[1] + '/' + bp2_info[1]
                
                # ignore modified base and GT mismatch
                if step_name.isupper() == False:
                    continue
                elif ('U' in step_name) == True:
                    continue
                elif any(s in step_name for s in ['GT','Gt','gT','gt','TG','Tg','tG','tg']):
                    continue
                
                 

                

                step_form = stem_forms[idx]

                o_ref1 = np.array(bp1_info[2])
                o_ref2 = np.array(bp2_info[2])
                M_ref1 = np.array(bp1_info[3:6])
                M_ref2 = np.array(bp2_info[3:6])

                origin, M_axis = rot_frame(o_ref1,o_ref2,M_ref1,M_ref2)
                origin = ",".join(format(x, ".3f") for x in origin)
                M_axis = ",".join(format(x, ".3f") for x in M_axis.flatten())

                step_idx += 1

                #print([pdbid,stem_idx,sform,step_idx,step_name,step_id,step_form,origin,M_axis])
                output.append([pdbid,stem_idx,sform,step_idx,step_name,step_id,step_form,origin,M_axis])

    output_df = pd.DataFrame(output[1:],columns=output[0])
    output_df.to_csv("DNA_step_library_frame.csv",index=False)

if __name__ == "__main__":
    main()


