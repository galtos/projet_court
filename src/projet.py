from pprint import pprint
import urllib.request
import os
import sys
import pbxplore as pbx

"""
Creer classe avec soit :
- nom structure
- nombre chaines
- nombre de frame
Dataframe :
- liste de liste avec assignation de l'alphabet
"""
input_topology,_ = [sys.argv[1],sys.argv[1]] #gro
input_trajectory,_ = [sys.argv[2],sys.argv[2]] #xtc
print(input_topology, input_trajectory)

class StructureAssign():
    def __init__(self, input_topology, input_trajectory):
        self.sa_md = []
        
        for chain_name, chain in pbx.chains_from_trajectory(input_trajectory, input_topology):
	        dihedrals = chain.get_phi_psi_angles()
	        pb_seq = pbx.assign(dihedrals)
	        self.sa_md.append(pb_seq)
	        #print('* {}'.format(chain_name))
	        #print('  {}'.format(pb_seq))

class Statistique():
    def __init__(self, sa_md):
        self.sa_md = sa_md

    def shannon_entropy(self):
        print(len(self.sa_md),len(self.sa_md[0]))
        #for i in range(len(self.sa_md[])):
        
        
        

#### MAIN ####

dt = StructureAssign(input_topology, input_trajectory)
print(dt.sa_md[1][1])
stat = Statistique(dt.sa_md)
stat.shannon_entropy()


