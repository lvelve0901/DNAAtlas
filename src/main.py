
#!/usr/bin/python


from gen_init_seed import *
#from gen_sample import gen_sample
from gen_sample_full import *
from copy import deepcopy

def main():
    seq = "GCATCGATTGGC"
    name = "A2dna"
    nstep = 100000
    T = 100.

    os.system("mkdir -p ../examples/%s_%s"%(name,T))

    helix = gen_init_seed(seq,name)
    mol = Mol('init_%s_rot.pdb'%name)


    for i in range(1,nstep+1):
        helix, mol = gen_sample_full(helix,mol,mode="MC",T=T)
       
        if i%100 == 0:
            print("Working on [%6d/%6d] step (%6.4f %%)"%(i,nstep,i*100./nstep))
            mol.renumber(1,1)
            mol.write('../examples/%s_%s/%s_output_%d.pdb'%(name,T,name,i/100))

if __name__ == "__main__":
    main()


