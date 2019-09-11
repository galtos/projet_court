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
input_topology, _ = [sys.argv[1],sys.argv[1]] #gro
input_trajectory, _ = [sys.argv[2],sys.argv[2]] #xtc
print(input_topology, input_trajectory)

class Structure_assign():
    def __init__(self, input_topology, input_trajectory):
	    for chain_name, chain in pbx.chains_from_trajectory(input_trajectory, input_topology):
	        dihedrals = chain.get_phi_psi_angles()
	        pb_seq = pbx.assign(dihedrals)
	        print('* {}'.format(chain_name))
	        print('  {}'.format(pb_seq))

#### MAIN ####

dt = Structure_assign(input_topology, input_trajectory)


