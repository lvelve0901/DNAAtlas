
#!/usr/bin/python


from gen_init_seed import *
#from gen_sample import gen_sample
from gen_sample_mc import *
from copy import deepcopy
from pprint import pprint

def main():
    seq = "GCATCGATTGGC"
    name = "A2dna"
    #helix = gen_init_seed(seq)
    #helix.write('test_output_0.pdb')
    #helix.show()
    
    #for i in range(1,10+1):
    #    helix = gen_sample(helix,mode="MC")
    #    helix.write('test_output_%d.pdb'%i)
    
    #helix.show()
    #helix.write('test_output.pdb')

    helix = gen_init_seed(seq,name)
    helix.show()

    nstep = 10000


    Ts = [10000.,1000.,100.,10.,1.,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,1e-3,1e-4,1e-5,1e-6]

    for T in Ts:
        swaps = []
        for i in range(1,nstep+1):
            helix = gen_sample_mc(helix,T=T)
            if helix.isSwap == True:
                swaps.append(1)
            else:
                swaps.append(0)

        pprint("Temp: %13.7f MC Pacc: %5d/%5d (%7.3f %%)"%(T,swaps.count(1),nstep,swaps.count(1)*100./nstep))

if __name__ == "__main__":
    main()



