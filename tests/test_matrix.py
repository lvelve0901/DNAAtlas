
import numpy as np
from common.matvec import lsqfit


M0 = np.identity(3)
M1 = np.array([[-0.03,0.004,1.],[-1.,0.003,-0.03],[-0.003,-1.,0.004]])

print(M1)
print(np.linalg.norm(M1,axis=1))
print(M1/np.linalg.norm(M1,axis=1))
M1 = M1/np.linalg.norm(M1,axis=1)

rot, rmsd = lsqfit(M1,M0)
print(rot)
print(rmsd)

M2 = np.array([0.867, -0.015, -0.498])
print(np.linalg.norm(M2))

