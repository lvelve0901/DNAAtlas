
import numpy as np
from common.matvec import lsqfit
from pdblib.num import *

o0 = np.array([1.,1.,1.])
x0 = np.array([1.,0.,0.])
y0 = np.array([0.,1.,0.])
z0 = np.array([0.,0.,1.])

o1 = np.array([-0.058, 6.885, -0.611])
x1 = np.array([-0.928, -0.041, -0.37])
y1 = np.array([-0.372, 0.073, 0.926])
z1 = np.array([-0.011, 0.997, -0.083])

o2 = np.array([0.606, 10.298, -1.2])
x2 = np.array([-0.952, -0.0, 0.305])
y2 = np.array([0.305, 0.059, 0.951])
z2 = np.array([-0.018, 0.998, -0.057])


A0 = np.array([x0,y0,z0])
A1 = np.array([x1,y1,z1])
A2 = np.array([x2,y2,z2])



rot, rmsd = lsqfit(A1,A0)
d = o1 - o0
print np.dot(A1,rot)
print np.dot(A2,rot)
o2_new = o2 - d
print o2_new
print np.dot(o2_new,rot)
dd = o2 - o1
print dd
print np.dot(dd,rot)

mol = Mol('3kz8_0_4.pdb')
newmol = Mol()
newmol.segs = [Segment(),Segment()]
newmol.segs[0].reses = mol.segs[0].reses[2:4]
newmol.segs[1].reses = mol.segs[1].reses[0:2]

#bp1 = [newmol.segs[0].reses[0],newmol.segs[1].reses[1]]
#bp2 = [newmol.segs[0].reses[1],newmol.segs[1].reses[0]]

#bp1_atoms = [atom for res in bp1 for atom in res.atoms]
#bp2_atoms = [atom for res in bp2 for atom in res.atoms]

#X = getmat(bp1_atoms)
#Y = getmat(bp2_atoms)

#print X[0]

m = getats(newmol)
M = getmat(m)

M -= d
M = np.dot(M,rot)

putmat(m,M)

#X = np.dot(X,rot)
#Y = np.dot(Y,rot)

#X -= d1
#Y -= d2

#putmat(bp1_atoms,X)
#putmat(bp2_atoms,Y)

newmol.write('3kz8_0_4_new.pdb')

