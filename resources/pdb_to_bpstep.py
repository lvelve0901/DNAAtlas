
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

# find the bp_id, bp_id_, bp_info and bp_info_ (both forward and reverse)
def find_bpid(nt1_id,nt2_id,bps,bps_idx):

    try:
        isSwap = False
        bp_id = nt1_id + "&" + nt2_id
        bp_id_ = nt2_id + "&" + nt1_id
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
    
    # forward bp direction
    origin  = bp['frame']['origin']
    x_axis  = np.array(bp['frame']['x_axis'])
    if isSwap == False:
        y_axis  = np.array(bp['frame']['y_axis'])
        z_axis  = np.array(bp['frame']['z_axis'])
        bp_name = "".join([bp['bp'][0],bp['bp'][-1]])
    else:
        y_axis  = np.array(bp['frame']['y_axis']) * -1
        z_axis  = np.array(bp['frame']['z_axis']) * -1
        bp_name = "".join([bp['bp'][0],bp['bp'][-1]])[::-1]

    # reverse bp direction
    origin_ = origin
    x_axis_ = x_axis
    y_axis_ = y_axis * -1 
    z_axis_ = z_axis * -1
    bp_name_ = bp_name[::-1]

    # normalization
    x_axis  = x_axis  / np.linalg.norm(x_axis)
    y_axis  = y_axis  / np.linalg.norm(y_axis)
    z_axis  = z_axis  / np.linalg.norm(z_axis)
    x_axis_ = x_axis_ / np.linalg.norm(x_axis_)
    y_axis_ = y_axis_ / np.linalg.norm(y_axis_)
    z_axis_ = z_axis_ / np.linalg.norm(z_axis_)

    bp_info = [bp_id,bp_name,origin,x_axis,y_axis,z_axis]
    bp_info_ = [bp_id_,bp_name_,origin_,x_axis_,y_axis_,z_axis_]

    return bp_id, bp_id_, bp_info, bp_info_

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
    d0 = o1 - o0
    o2_new = np.dot(o2 - d0,rot)
    
    A1_new = np.dot(A1,rot)
    A2_new = np.dot(A2,rot)

    return o2_new, A2_new, d0, rot

# input two bp_info and perform origin and axis Matrix rotation on bp1 then bp2
def rot_matrix(bp1_info,bp2_info):

    o_ref1 = np.array(bp1_info[2])
    o_ref2 = np.array(bp2_info[2])
    M_ref1 = np.array(bp1_info[3:6])
    M_ref2 = np.array(bp2_info[3:6])
    
    origin, M_axis, d0, rot = rot_frame(o_ref1,o_ref2,M_ref1,M_ref2)
    origin = ",".join(format(x, ".3f") for x in origin)
    M_axis = ",".join(format(x, ".3f") for x in M_axis.flatten())
  
    return origin, M_axis


def main():

    step_idx = 0
    jsondir = '/mnt/hs189/NAfinder_3.0/Crystal/Json'
    fragdir = '/mnt/hs189/DNAAtlas/resources/refine_frag/refine_stem'
    pdblist = []

    df = pd.read_csv('/mnt/hs189/NAfinder_3.0/doc/NDB/PDB_NDB/Final_library_DNA_Representative.csv')
    pdblist = df.loc[(df.hasLigand == False) & (df.reso < 2.0)]['pdbid'].tolist()
    pdblist.sort()

    output = [['pdbid','stem_idx','stem_len','sform','step_idx','step_name','step_id','step_form','origin','M_axis']]

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
            stem_len = len(stem_pair)
            stem_file = os.path.join(fragdir,pdbid + '_' + str(j) + '_' + str(stem_len) + '.pdb')

            stem_forms = stem['helix_form']
            sform = gen_sform(stem_forms)

            # ignore all non-Bform DNA
            if sform != 'B':
                continue
    
            bp_pairs = []
            
            for pair in stem_pair:
                bp_pairs.append((pair['nt1'],pair['nt2']))
            
            # skip terminal base pairs
            for idx in range(1,len(bp_pairs)-2):


                bp1_id, bp1_id_, bp1_info, bp1_info_ = find_bpid(bp_pairs[idx][0],  bp_pairs[idx][1],  bps,bps_idx)
                bp2_id, bp2_id_, bp2_info, bp2_info_ = find_bpid(bp_pairs[idx+1][0],bp_pairs[idx+1][1],bps,bps_idx)
                step_id  = bp1_id  + '&&' + bp2_id
                step_id_ = bp2_id_ + '&&' + bp1_id_
                step_name  = bp1_info[1]  + '/' + bp2_info[1]
                step_name_ = bp2_info_[1] + '/' + bp1_info_[1]

                # ignore modified base and GT mismatch
                if step_name.isupper() == False:
                    continue
                elif ('U' in step_name) == True:
                    continue
                elif any(s in step_name for s in ['GT','Gt','gT','gt','TG','Tg','tG','tg']):
                    continue

                step_form = stem_forms[idx]

                origin,  M_axis,  = rot_matrix(bp1_info, bp2_info )
                origin_, M_axis_, = rot_matrix(bp2_info_,bp1_info_)

                step_idx += 1
                output.append([pdbid,stem_idx,stem_len,sform,step_idx,step_name,step_id,step_form,origin,M_axis])

                step_idx += 1
                output.append([pdbid,stem_idx,stem_len,sform,step_idx,step_name_,step_id_,step_form,origin_,M_axis_])

    output_df = pd.DataFrame(output[1:],columns=output[0])
    output_df.to_csv("DNA_step_library_frame.csv",index=False)

if __name__ == "__main__":
    main()


