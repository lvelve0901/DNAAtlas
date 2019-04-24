from pdblib.num import *
from common.base import range2list
import pdb

atns = ("C1'","C2'","C5'","O3'","O4'","N1","N3","N9","C4","C6")

#===============================================================================
def supimpose(m1, m0, segl1, resl1, segl0, resl0):
    '''superimpose m1 on m0
    '''
    
    res_list1 = [m1.segs[segid].reses[resid] for segid, resid in zip(segl1,resl1)]
    res_list0 = [m0.segs[segid].reses[resid] for segid, resid in zip(segl0,resl0)]

    ats_m1 = []
    for res in res_list1:
        for at in atns:
            at1 = res.getat(at)
            if at1:
                ats_m1.append(res.getat(at))
    ats_m0 = []
    for res in res_list0:
        for at in atns:
            at0 = res.getat(at)
            if at0:
                ats_m0.append(res.getat(at))
    if len(ats_m1) != len(ats_m0):
        #pdb.set_trace()
        print('supimpose: The overlaid regions have different size!')
        print('  source: ' + ' '.join(map(lambda x: '%d-%s'%(x.resi,x.name),
              ats_m1)))
        print('  target: ' + ' '.join(map(lambda x: '%d-%s'%(x.resi,x.name),
              ats_m0)))
        exit(1)
    # align
    if (isinstance(m1, list) and isinstance(m1[0], Residue)):
        m1 = reduce(add, map(lambda x: x.atoms, m1), [])
    M = getmat(m1)
    A = getmat(ats_m1)
    B = getmat(ats_m0)
    mc_A = mean(A, 0)
    A -= mc_A
    mc_B = mean(B, 0)
    B -= mc_B
    rot,rmsd = matvec.lsqfit(A, B)
    if rot is not None:
        M -= mc_A
        M = dot(M,rot) + mc_B
        putmat(m1, M)
        return rmsd
    else:
        raise Exception, 'The alignment failed!'

#===============================================================================
def res_del(mol, rg):
    ri = range2list(rg)
    for seg in mol.segs:
        seg.reses = [res for res in seg.reses if res.resi not in ri]

#===============================================================================
def res_ins(sg, reses, resi):
    idx = sg.getindex(resi)
    sg.reses[idx+1:idx+1] = reses

#===============================================================================
def getreslist(sg, rg):
    rg = range2list(rg)
    return [res for res in sg.reses if res.resi in rg]

#===============================================================================
def glue(mol1, mol2):
    '''Glue mol2 to mol1
    '''
    nsg1 = len(mol1.segs)
    nsg2 = len(mol2.segs)
    if nsg1==1 and nsg2==1:
        mol1.segs[0].reses = mol1.segs[0].reses + mol2.segs[0].reses
    elif nsg1==1 and nsg2==2:
        mol1.segs[0].reses = mol2.segs[0].reses + mol1.segs[0].reses \
                           + mol2.segs[1].reses
    elif nsg1==2 and nsg2==1:
        mol1.segs[0].reses = mol1.segs[0].reses + mol2.segs[0].reses \
                           + mol1.segs[1].reses
        mol1.segs[1:] = []
    elif nsg1==2 and nsg2==2:
        mol1.segs[0].reses = mol1.segs[0].reses + mol2.sges[0].reses
        mol1.segs[1].reses = mol2.segs[1].reses + mol1.sges[1].reses
    else:
        print('Error: mol1 has %d nt, mol2 has %d nt!'%(nsg1,nsg2))
        exit(1)
    #mol1.renumber(1, 1)
