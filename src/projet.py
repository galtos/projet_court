import sys
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import pbxplore as pbx

SA_SIZE = 17
SA_LETTERS = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k",\
              "l", "m", "n", "o", "p", "Z"]

#class Configuration

class AssignationPbxplore():
    """
    This class assign structural letters to a protein structure with pbxplore 
    to each frame of a molecular dynamics. The class takes 3 arguments:
    
    - the topology of the structure
    - the trajectory of the molecular dynamics
    - the argument to know if a fasta output file is needed
    
    the set of structural letters is stored in : md_sa_seq
    """
    def __init__(self, input_topology, input_trajectory, args):
        self.md_sa_seq = []
        if args.f:
            file_output_fasta = open(args.f,"w")
        for chain_name, chain in pbx.chains_from_trajectory(input_trajectory,\
                                                            input_topology):
	        dihedrals = chain.get_phi_psi_angles()
	        pb_seq = pbx.assign(dihedrals)
	        self.md_sa_seq.append(pb_seq)
	        if args.f:
	            file_output_fasta.write(">{}\n".format(chain_name))
	            file_output_fasta.write("{}\n".format(pb_seq))
        if args.f:
            file_output_fasta.close()

class Statistique():
    """
    This class contain the differents methods to perform the statistics 
    analysis of the set of structural letters of the structure during the 
    trajectory.
    This class actually compute the mutual information.
    It takes one argument:
    - the set of structural letters of the(md_sa_seq)
    """
    def __init__(self, md_sa_seq):
        self.md_sa_seq = md_sa_seq
        
       
    def matrice_p_a(self):
        """
        Method calculate the probability for a structural letter to occur in 
        the molecular dynamics. P(a)
        Return a matrix 1*SA_SIZE with the probability for each letters
        to occur.
        """
        mat_probabilities = pd.DataFrame(np.zeros(SA_SIZE), index = SA_LETTERS)

        for i in range(len(self.md_sa_seq)):
            for j in range(len(self.md_sa_seq[0])):
                mat_probabilities[0][self.md_sa_seq[i][j]] += 1
                
        return(mat_probabilities/(len(self.md_sa_seq)*len(self.md_sa_seq[0])))
        
    
    def matrice_p_ab(self):
        """
        Method calculate the probability for a two set of structural letters to 
        occur one after an other during the molecular dynamics. It count the 
        number of time a set of structural letts is seen between two frame of 
        the molecular dynamics.
        ex: in the first frame in the position 1 the letter "a" is seen and in 
        the frame 2 at the position 1 the letter "b" is seen.
        the method count these occurences and divide them with the number of 
        positions in the structure * the number of frames of the molecular
        dynamics.
        Return a matrix of size SA_SIZE*SA_SIZE with the probability for each 
        pair of letters. P(a,b)
        """
        mat_count_occurence = pd.DataFrame(                       \
                                     np.zeros((SA_SIZE, SA_SIZE)),\
                                     columns = SA_LETTERS, index = SA_LETTERS) 

        for i in range(len(self.md_sa_seq)-1):
            for j in range(len(self.md_sa_seq[0])):
                mat_count_occurence[self.md_sa_seq[i][j]]\
                                   [self.md_sa_seq[i+1][j]] += 1
        
        ####probabilité jointe ####
        return(mat_count_occurence/\
              (len(self.md_sa_seq)*len(self.md_sa_seq[0])))
    def mutual_information_i(self, i, j, mat_p_a, mat_p_ab):
        I = 0
        for k in range(len(self.md_sa_seq)): 
                if mat_p_ab[self.md_sa_seq[k][i]]\
                           [self.md_sa_seq[k][j]] == 0 \
                or mat_p_a[0][self.md_sa_seq[k][i]] == 0 \
                or mat_p_a[0][self.md_sa_seq[k][j]] == 0:  
                    I += 0
                else:               
                    I += mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]]*\
                 np.log2(mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]]/\
                (mat_p_a[0][self.md_sa_seq[k][i]]*\
                 mat_p_a[0][self.md_sa_seq[k][j]]))
    
    
        return I
    def mutual_information(self):
        """
        Compute the mutual information for each pair of positions in the 
        structure.
        return the matrix of size (number of position inf the structure)**2 with
        each mutual information calculate.
        """
        MI = pd.DataFrame(np.zeros((len(self.md_sa_seq[0]),\
                                    len(self.md_sa_seq[0]))))
        
        mat_p_a = self.matrice_p_a()
        mat_p_ab = self.matrice_p_ab()
        
        for i in range(len(self.md_sa_seq[0])):
            for j in range(0, len(self.md_sa_seq[0])):
                MI[i][j] = self.mutual_information_i(i, j, mat_p_a, mat_p_ab)
        return(MI)
        
class Visualisation():
    def __init__(self, matrice_MI):
        self.matrice_MI = matrice_MI
    def write_mi_csv(self, args):
        file_csv = open(args.omi,"w")
        for i in range(len(self.matrice_MI)):
            for j in range(len(self.matrice_MI)):
                file_csv.write("{},".format(self.matrice_MI[i][j]))
            file_csv.write("\n")
        file_csv.close()       
    def visualize_matrice_mi(self):

        sns.heatmap(self.matrice_MI)
        plt.show()

    def color_palette(self):
        color =[[ 0.00,  0.00,  1.00],\
                [ 0.24,  0.16,  0.75],\
                [ 0.36,  0.24,  0.63],\
                [ 0.53,  0.34,  0.47],\
                [ 0.77,  0.50,  0.22],\
                [ 1.00,  0.63,  0.00],\
                [ 1.00,  0.50,  0.00],\
                [ 1.00,  0.37,  0.00],\
                [ 1.00,  0.24,  0.00],\
                [ 1.00,  0.10,  0.00],\
                [ 1.00,  0.00,  0.00]]
        return color
        
            
    def show_pymol_network(self, input_pdb):
        input_pymol = open("pymol_tmp.pml","w")
        input_pymol.write("load "+input_pdb+"\n")
        input_pymol.write("set dash_gap, 0\n")
        input_pymol.write("set dash_width, 2\n")
        input_pymol.write("preset.pretty(selection='all')\n")
        
        line_treshold = np.max(np.max(self.matrice_MI))/2
        print(line_treshold)
        
        color = np.max(self.matrice_MI)/np.max(np.max(self.matrice_MI))

        color = 10*(color - color.min()) / (color.max() - color.min())
        palette = self.color_palette()
        print(palette)
        print(color)
        for i in range(len(self.matrice_MI)):
            input_pymol.write("set_color color_{} =   [{} , {} , {}]\n".format(i,\
                                palette[int(color[i])][0],\
                                palette[int(color[i])][1],\
                                palette[int(color[i])][2]))
            input_pymol.write("color color_{}, resi {}\n".format(i,i+1))
        
            for j in range(i, len(self.matrice_MI)):
                if self.matrice_MI[i][j] > line_treshold:

                    input_pymol.write("distance {}_{},".format(i + 1,j + 1)+\
            "name CA and resi {}, ".format(i + 1)+ "name CA and resi {}\n".format(j + 1))
                    input_pymol.write("hide labels, "+"{}_{}\n".format(i + 1,j + 1))

                    

        input_pymol.close()
        os.system("pymol -p pymol_tmp.pml")
 


#### MAIN ####
if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", nargs = "?", help = "output a fasta file containing the generated structural alphabet", const = "pbxplore_sequence.fasta")
    parser.add_argument("-py", nargs = "?", help = "generate the structure network according to the MI matrix on pymol\ntake in argument the pdb file of the structure")
    parser.add_argument("-omi", nargs = "?", help = "generate an output file containing the mutual information matrix ", const = "mi_matrix.csv")
    
    parser.add_argument("input_topology", help = "number of molecular dynamics to analyse")
    parser.add_argument("input_trajectory", help = "number of molecular dynamics to analyse")
    
    args = parser.parse_args()

    dt = AssignationPbxplore(args.input_topology, args.input_trajectory, args)

    stat = Statistique(dt.md_sa_seq)

    #stat.matrice_p_a()
    #stat.matrice_p_ab()
    MI = stat.mutual_information()
    #MI2 = [[1,2,3],[1,1,1],[3,3,3]]
    mat_visual = Visualisation(MI)
    ##mat_visual.visualize_matrice_mi()
    mat_visual.write_mi_csv(args)
    print(MI)
    if args.py:
        mat_visual.show_pymol_network(args.py)








