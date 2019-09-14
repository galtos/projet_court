from pprint import pprint
import urllib.request
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pymol
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
input_pdb = sys.argv[3]

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

        nSA = pd.DataFrame(np.zeros(AS_SIZE), index = AS_LETTERS)
      
        nLetters = len(self.sa_md[0])
        nFrames = len(self.sa_md)

        mat_probabilities = pd.DataFrame(np.zeros(AS_SIZE), index = AS_LETTERS)

        for i in range(len(self.sa_md)):
            for j in range(len(self.sa_md[0])):
                mat_probabilities[0][self.sa_md[i][j]] += 1
                
        return(mat_probabilities/(len(self.sa_md)*len(self.sa_md[0])))
        
     

    def matrice_p_ab(self):

        mat_count_occurence = pd.DataFrame(np.zeros((AS_SIZE, AS_SIZE)), columns = AS_LETTERS, index = AS_LETTERS) 

        for i in range(len(self.sa_md)-1):
            for j in range(len(self.sa_md[0])):
                mat_count_occurence[self.sa_md[i][j]][self.sa_md[i+1][j]] += 1
        
         
        ####probabilité jointe ####
        print(mat_count_occurence/(len(self.sa_md)*len(self.sa_md[0])))
        return(mat_count_occurence/(len(self.sa_md)*len(self.sa_md[0])))

    def mutual_information(self):
        print(len(self.sa_md),len(self.sa_md[0]))
        
        MI = pd.DataFrame(np.zeros((len(self.sa_md[0]),len(self.sa_md[0]))))
        
        mat_p_a = self.matrice_p_a()
        mat_p_ab = self.matrice_p_ab()
        
        for i in range(len(self.sa_md[0])):
            for j in range(0, len(self.sa_md[0])): ## pour remplir MI
                I = 0.
                for k in range(len(self.sa_md)):#Pour iterer sur les frame FIXME 
                    if mat_p_ab[self.sa_md[k][i]][self.sa_md[k][j]] == 0 or mat_p_a[0][self.sa_md[k][i]] == 0 or mat_p_a[0][self.sa_md[k][j]] == 0:  
                        I += 0
                    else:               
                        I += mat_p_ab[self.sa_md[k][i]][self.sa_md[k][j]]*\
                    np.log2(mat_p_ab[self.sa_md[k][i]][self.sa_md[k][j]]/\
                    (mat_p_a[0][self.sa_md[k][i]]*mat_p_a[0][self.sa_md[k][j]]))
       
                MI[i][j] = I
        return(MI)
        
class Visualisation():
    def __init__(self, matrice_MI):
        self.matrice_MI = matrice_MI
    def visualize_matrice_mi(self):

        sns.heatmap(self.matrice_MI)
        plt.show()
    def show_pymol_network(self):
        pymol.finish_launching(["pymol","-q"]) 
        for i in range(len(self.matrice_MI)):
            for j in range(len(self.matrice_MI)):
                if self.matrice_MI[i][j] > SEUIL:
                    cmd.distance("resi"+string(i),"resi"+string(j))
#### MAIN ####

dt = StructureAssign(input_topology, input_trajectory)
print(dt.sa_md[1][1])
stat = Statistique(dt.sa_md)

#stat.matrice_p_a()
#stat.matrice_p_ab()
#MI = stat.mutual_information()
#MI = [[1,2,3],[1,1,1],[3,3,3]]
mat_visual = Visualisation(MI)
#mat_visual.visualize_matrice_mi()
mat_visual.show_pymol_network()








