from pprint import pprint
import urllib.request
import os
import sys
import numpy as np
import pandas as pd
import pbxplore as pbx

AS_SIZE = 17
AS_LETTERS = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","Z"]

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

        #compte le nombre de lettre du SA dans la MD 

        nSA = {"a":0,"b":0,"c":0,"d":0,"e":0,"f":0,"g":0,"h":0,"i":0,"j":0,"k":0,"l":0,"m":0,"n":0,"o":0,"p":0,"Z":0}
      
        nLetters = len(self.sa_md[0])
        nFrames = len(self.sa_md)

        mat_probabilities = pd.DataFrame(np.zeros((AS_SIZE,len(self.sa_md[0]))), index = AS_SIZE)
        print(mat_probabilities)

        for i in range(len(self.sa_md[0])):
            for j in range(len(self.sa_md)):
                nSA[self.sa_md[j][i]] +=1

            for k in nSA:
                if nSA[k] == 0:
                    
                    mat_probabilities[k][i] = 0
                else:
                    mat_probabilities[k][i] = nSA[k]/len(self.sa_md) 

                nSA[k] = 0

        print( mat_probabilities)
        

         
"""
    def matrice_p_ab(self, A, B):

        mat_count_occurence = pd.DataFrame(np.zeros((17,17)), columns = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","Z"],\
        index = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","Z"]) # FIXME : taille matrice en dur

        for i in range(len(A)):
            mat_count_occurence[A[i]][B[i]] += 1
        print(mat_count_occurence) 
        ####probabilité jointe ####
        
        print(mat_count_occurence/(len(self.sa_md[0])/2))
        
    def mutual_information(self):
        print(len(self.sa_md),len(self.sa_md[0]))
        
        MI = pd.DataFrame(np.zeros(len(self.sa_md[0]),len(self.sa_md[0])))
        mat_pa = self.matrice_p_a()
        
        for i in range(len(self.sa_md)):
            for j in range(len(self.sa_md)):
                mat_p_ab = matrice_p_ab(self.sa_md[0:len(self.sa_md)][i],\
                                       self.sa_md[0:len(self.sa_md)][j])
                for k in SA_LETTERS:
                    for l in range(SA_LETTERS):    
                        I = mat_p_ab[k][l]*np.log2(mat_p_ab[k][l]/(mat_p_a[k][l]*mat_p_a[k][l]))
        
        print(MI)
        """
#### MAIN ####

dt = StructureAssign(input_topology, input_trajectory)
print(dt.sa_md[1][1])
stat = Statistique(dt.sa_md)

stat.matrice_p_a()
#stat.mutual_information()




