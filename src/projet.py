from pprint import pprint
import urllib.request
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#import pymol
import pbxplore as pbx

SA_SIZE = 17
SA_LETTERS = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","Z"]



class AssignationPbxplore():
    def __init__(self, input_topology, input_trajectory):
        self.md_sa_seq = []
        
        for chain_name, chain in pbx.chains_from_trajectory(input_trajectory, input_topology):
	        dihedrals = chain.get_phi_psi_angles()
	        pb_seq = pbx.assign(dihedrals)
	        self.md_sa_seq.append(pb_seq)


class Statistique():
    def __init__(self, md_sa_seq):
        self.md_sa_seq = md_sa_seq
        
       
    def matrice_p_a(self):

        #compte le nombre de lettre du SA dans la MD 

        nSA = pd.DataFrame(np.zeros(SA_SIZE), index = SA_LETTERS)
      
        nLetters = len(self.md_sa_seq[0])
        nFrames = len(self.md_sa_seq)

        mat_probabilities = pd.DataFrame(np.zeros(SA_SIZE), index = SA_LETTERS)

        for i in range(len(self.md_sa_seq)):
            for j in range(len(self.md_sa_seq[0])):
                mat_probabilities[0][self.md_sa_seq[i][j]] += 1
                
        return(mat_probabilities/(len(self.md_sa_seq)*len(self.md_sa_seq[0])))
        
     

    def matrice_p_ab(self):

        mat_count_occurence = pd.DataFrame(np.zeros((SA_SIZE, SA_SIZE)), columns = SA_LETTERS, index = SA_LETTERS) 

        for i in range(len(self.md_sa_seq)-1):
            for j in range(len(self.md_sa_seq[0])):
                mat_count_occurence[self.md_sa_seq[i][j]][self.md_sa_seq[i+1][j]] += 1
        
         
        ####probabilité jointe ####
        print(mat_count_occurence/(len(self.md_sa_seq)*len(self.md_sa_seq[0])))
        return(mat_count_occurence/(len(self.md_sa_seq)*len(self.md_sa_seq[0])))

    def mutual_information(self):
        print(len(self.md_sa_seq),len(self.md_sa_seq[0]))
        
        MI = pd.DataFrame(np.zeros((len(self.md_sa_seq[0]),len(self.md_sa_seq[0]))))
        
        mat_p_a = self.matrice_p_a()
        mat_p_ab = self.matrice_p_ab()
        
        for i in range(len(self.md_sa_seq[0])):
            for j in range(0, len(self.md_sa_seq[0])): ## pour remplir MI
                I = 0.
                for k in range(len(self.md_sa_seq)):#Pour iterer sur les frame FIXME 
                    if mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]] == 0 or mat_p_a[0][self.md_sa_seq[k][i]] == 0 or mat_p_a[0][self.md_sa_seq[k][j]] == 0:  
                        I += 0
                    else:               
                        I += mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]]*\
                    np.log2(mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]]/\
                    (mat_p_a[0][self.md_sa_seq[k][i]]*mat_p_a[0][self.md_sa_seq[k][j]]))
       
                MI[i][j] = I
        return(MI)
        
class Visualisation():
    def __init__(self, matrice_MI):
        self.matrice_MI = matrice_MI
    def visualize_matrice_mi(self):

        sns.heatmap(self.matrice_MI)
        plt.show()
"""
    def show_pymol_network(self):
        pymol.finish_launching(["pymol","-q"]) 
        for i in range(len(self.matrice_MI)):
            for j in range(len(self.matrice_MI)):
                if self.matrice_MI[i][j] > SEUIL:
                    cmd.distance("resi"+string(i),"resi"+string(j))
"""
#### MAIN ####

input_topology,_ = [sys.argv[1],sys.argv[1]] #gro
input_trajectory,_ = [sys.argv[2],sys.argv[2]] #xtc
#input_pdb = sys.argv[3]

print(input_topology, input_trajectory)

dt = AssignationPbxplore(input_topology, input_trajectory)
print(dt.md_sa_seq[1][1])
stat = Statistique(dt.md_sa_seq)

#stat.matrice_p_a()
#stat.matrice_p_ab()
MI = stat.mutual_information()
#MI = [[1,2,3],[1,1,1],[3,3,3]]
mat_visual = Visualisation(MI)
mat_visual.visualize_matrice_mi()
#mat_visual.show_pymol_network()








