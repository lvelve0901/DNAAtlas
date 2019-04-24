#!/usr/bin/python

from pdblib.base import *
from os import listdir

dn = '../examples/A2dna_0.6_movie'
fns = listdir(dn)
fns.sort(key=lambda x: int(x.split('_')[-1][:-4]))

pdb = Pdb()
for i,fn in enumerate(fns):
    mol = Mol(dn+'/'+fn)
    mol.mdid = i+1
    pdb.mds.append(mol)

pdb.write('A2_movie.pdb')

