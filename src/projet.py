from pprint import pprint
import urllib.request
import os
import sys
import numpy as np
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
        
        
    def matrice_p_a(self):
        """
        compte le nombre de lettre du SA dans la MD 
        """
        nSA = {"a":0,"b":0,"c":0,"d":0,"e":0,"f":0,"g":0,"h":0,"i":0,"j":0,"k":0,"l":0,"m":0,"n":0,"o":0,"p":0,"Z":0}
        
        nLetters = len(self.sa_md[0])
        nFrames = len(self.sa_md)
        mat_probabilities = []

        for i in range(len(self.sa_md)):
            mat_probabilities.append([])
            for j in range(len(self.sa_md[0])):
                nSA[self.sa_md[i][j]] +=1
            for k in nSA:
                if nSA[k] == 0:
                    
                    mat_probabilities[i].append(0)
                else:
                    mat_probabilities[i].append(nSA[k]/nLetters)   

                nSA[k] = 0
            


        print(nSA)  
        
        
          
    def mutual_information(self):
        print(len(self.sa_md),len(self.sa_md[0]))
        #for i in range(len(self.sa_md[])):
        
        
        

#### MAIN ####

dt = StructureAssign(input_topology, input_trajectory)
print(dt.sa_md[1][1])
stat = Statistique(dt.sa_md)
stat.mutual_information()
stat.matrice_p_a()

