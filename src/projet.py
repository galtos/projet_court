#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read a molecular dynamics and assign for each frames protein blocs[1]
The program then compute the mutual information matrix and represent the matrix
in the form of a heatmap and/or a structure network with pymol[2].
This program is inspired by GSATools[3].

2019 - Guillaume OLLITRAULT

[1]
Barnoud J, Santuz H, Craveur P, Joseph AP, Jallu V, de Brevern AG, Poulain P,
PBxplore: a tool to analyze local protein structure and deformability with
Protein Blocks.
PeerJ 5:e4013 https://doi.org/10.7717/peerj.4013 (2017).
[2]
The PyMOL Molecular Graphics System, Version 1.8.4.0 Schrödinger, LLC.
[3]
Pandini A, Fornili A, Fraternali F, Kleinjung J
GSATools: analysis of allosteric communication and functional local motions
using a Structural Alphabet
Bioinformatics 29(16):2053-2055, 2013
"""

import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pbxplore as pbx

SA_SIZE = 17
SA_LETTERS = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k",\
              "l", "m", "n", "o", "p", "Z"]

class AssignationPbxplore(object):
    """
    This class assign structural letters to a protein structure with pbxplore
    to each frame of a molecular dynamics. The class takes 3 arguments:

    - the topology of the structure
    - the trajectory of the molecular dynamics
    - the argument to know if a fasta output file is needed

    the set of structural letters is stored in : md_sa_seq
    """
    def __init__(self, args):
        self.md_sa_seq = []
        if args.f:
            file_output_fasta = open(args.f, "w")
        for chain_name, chain in pbx.chains_from_trajectory(\
                                                        args.input_trajectory,\
                                                        args.input_topology):
            dihedrals = chain.get_phi_psi_angles()
            pb_seq = pbx.assign(dihedrals)
            self.md_sa_seq.append(pb_seq)
            if args.f:
                file_output_fasta.write(">{}\n".format(chain_name))
                file_output_fasta.write("{}\n".format(pb_seq))
        if args.f:
            file_output_fasta.close()

class Statistics(object):
    """
    This class contains the different methods to perform the statistics
    analysis of the set of structural letters of the structure during the
    trajectory.
    This class compute the mutual information.
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
        mat_probabilities = pd.DataFrame(np.zeros(SA_SIZE), index=SA_LETTERS)

        for i in range(len(self.md_sa_seq)):
            for j in range(len(self.md_sa_seq[0])):
                mat_probabilities[0][self.md_sa_seq[i][j]] += 1

        return mat_probabilities/(len(self.md_sa_seq)*len(self.md_sa_seq[0]))

    def matrice_p_ab(self):
        """
        Method calculate the probability for a two set of structural letters to
        occur one after an other during the molecular dynamics. It count the
        number of time a set of structural letters is seen between two frame of
        the molecular dynamics.
        ex: in the first frame in the position 1 the letter "a" is seen and in
        the frame 2 at the position 1 the letter "b" is seen.
        the method count these occurences and divide them with the number of
        positions in the structure * the number of frames of the molecular
        dynamics.
        Return a matrix of size SA_SIZE*SA_SIZE with the probability for each
        pair of letters. P(a,b)
        """
        mat_count_occurence = pd.DataFrame(                            \
                                          np.zeros((SA_SIZE, SA_SIZE)),\
                                          columns=SA_LETTERS, index=SA_LETTERS)

        for i in range(len(self.md_sa_seq)-1):
            for j in range(len(self.md_sa_seq[0])):
                mat_count_occurence[self.md_sa_seq[i][j]]\
                                   [self.md_sa_seq[i+1][j]] += 1

        return mat_count_occurence/(len(self.md_sa_seq)*len(self.md_sa_seq[0]))

    def mutual_information_i(self, i, j, mat_p_a, mat_p_ab):
        """
        Compute the mutual information for two sets of positions in the
        structure. MI(X;Y) =  ∑ P(x;y)*log(P(x;y)/(P(x)*P(y)))
        """
        mutual_information = 0
        for k in range(len(self.md_sa_seq)):
            if mat_p_ab[self.md_sa_seq[k][i]]\
                       [self.md_sa_seq[k][j]] == 0 \
            or mat_p_a[0][self.md_sa_seq[k][i]] == 0 \
            or mat_p_a[0][self.md_sa_seq[k][j]] == 0:
                mutual_information += 0
            else:
                mutual_information += \
                    mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]]*\
           np.log2(mat_p_ab[self.md_sa_seq[k][i]][self.md_sa_seq[k][j]]/\
           (mat_p_a[0][self.md_sa_seq[k][i]]*mat_p_a[0][self.md_sa_seq[k][j]]))

        return mutual_information

    def mutual_information(self):
        """
        Compute the mutual information for each pair of positions in the
        structure.
        return the matrix of size (number of position inf the structure)**2
        with each mutual information calculate.
        """
        mutual_information = pd.DataFrame(np.zeros((len(self.md_sa_seq[0]),\
                                    len(self.md_sa_seq[0]))))

        mat_p_a = self.matrice_p_a()
        mat_p_ab = self.matrice_p_ab()

        for i in range(len(self.md_sa_seq[0])):
            for j in range(0, len(self.md_sa_seq[0])):
                mutual_information[i][j] = self.mutual_information_i(\
                                                i, j, mat_p_a, mat_p_ab)
        return mutual_information

class Visualization(object):
    """
    Contains the methods for the visualization of the mutual information matrix
    with heatmap and a network representation with pymol. This class also allow
    the user to write a cvs file with the mutual information matrix.
    """
    def __init__(self, matrice_mi):
        self.matrice_mi = matrice_mi

    def write_mi_csv(self, args):
        """
        Write a csv file containing the mutual inforamtion matrix
        """
        file_csv = open(args.omi, "w")
        for i in range(len(self.matrice_mi)):
            for j in range(len(self.matrice_mi)):
                file_csv.write("{},".format(self.matrice_mi[i][j]))
            file_csv.write("\n")
        file_csv.close()
    def visualize_matrice_mi(self):
        """
        Represent the mutual information matrix in form of a heatmap
        """
        ax_heatmap = plt.axes()
        sns.heatmap(self.matrice_mi, ax=ax_heatmap)
        ax_heatmap.set_title(\
       "Heatmap of the compute mutual information for the different positions")
        plt.show()

    def color_palette(self):
        """
        r,g,b numbers for the color palette to visualize with pymol
        the important residues involved in the network.
        The colors goes from the blue to red.
        """
        color = [[0.00, 0.00, 1.00],\
                 [0.24, 0.16, 0.75],\
                 [0.36, 0.24, 0.63],\
                 [0.53, 0.34, 0.47],\
                 [0.77, 0.50, 0.22],\
                 [1.00, 0.63, 0.00],\
                 [1.00, 0.50, 0.00],\
                 [1.00, 0.37, 0.00],\
                 [1.00, 0.24, 0.00],\
                 [1.00, 0.10, 0.00],\
                 [1.00, 0.00, 0.00]]
        return color

    def show_pymol_network(self, input_pdb, line_threshold):
        """
        launch pymol and read a pml file which color the residues according
        their significance in the network made with the mutual information.
        Put a line between the important residues (>threshold).
        the threshold is the maximum value of the matrix divided by 2 but
        can be given by the user.
        """
        input_pymol = open("pymol_tmp.pml", "w")
        input_pymol.write("load "+input_pdb+"\n")
        input_pymol.write("set dash_gap, 0\n")
        input_pymol.write("set dash_width, 2\n")
        input_pymol.write("preset.pretty(selection='all')\n")

        color = np.max(self.matrice_mi)/np.max(np.max(self.matrice_mi))

        color = 10*(color - color.min()) / (color.max() - color.min())
        palette = self.color_palette()

        for i in range(len(self.matrice_mi)):
            input_pymol.write("set_color color_{} =   [{} , {} , {}]\n"\
                                .format(i,\
                                palette[int(color[i])][0],\
                                palette[int(color[i])][1],\
                                palette[int(color[i])][2]))
            input_pymol.write("color color_{}, resi {}\n".format(i, i+1))

            for j in range(i, len(self.matrice_mi)):
                if self.matrice_mi[i][j] > line_threshold:

                    input_pymol.write("distance {}_{},".format(i + 1, j + 1)+\
                                      "name CA and resi {}, ".format(i + 1) +\
                                      "name CA and resi {}\n".format(j + 1))
                    input_pymol.write("hide labels, "+"{}_{}\n".format(i + 1,\
                                                                       j + 1))

        input_pymol.close()
        os.system("pymol -p pymol_tmp.pml")


def analysis_mutual_information():
    """
    Perform the anylisis by computing the mutual information and dealing with
    the different arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", nargs="?",\
        help="output a fasta file containing the generated structural" + \
             " alphabet",\
        const="pbxplore_sequence.fasta")

    parser.add_argument("-py", nargs="?", \
        help="generate the structure network according to the MI matrix" + \
             " on pymol\ntake in argument the pdb file of the structure\n" + \
             "save a .pml file that can be reused with pymol to see again\n" +\
             "the visualization." + \
             "warning: your residue numbering have to start at 1 in your" + \
             " pdb file.")
    parser.add_argument("-omi", nargs="?",\
        help="generate an output file containing the mutual information" + \
             " matrix ",\
        const="mi_matrix.csv")

    parser.add_argument("-hmap", \
        action="store_true",\
        help="show the heatmap of the mutual information matrix")

    parser.add_argument("input_topology",\
        help="topology file of the molecular dynamics")
    parser.add_argument("input_trajectory",\
        help="trajectory file of the molecular dynamics")

    args = parser.parse_args()

    pbxplore_seq_md = AssignationPbxplore(args)

    stats_mi = Statistics(pbxplore_seq_md.md_sa_seq)

    matrix_mi = stats_mi.mutual_information()

    mat_visual = Visualization(matrix_mi)

    if args.hmap:
        mat_visual.visualize_matrice_mi()

    if args.omi:
        mat_visual.write_mi_csv(args)

    if args.py:
        change_threshold = input(\
                    "You asked for the network visualization with pymol\n" + \
                    "do you want to change the threshold for the network " + \
                    "representation ? (y/n) :")
        while change_threshold != "y" and change_threshold != "n":
            change_threshold = input(\
                    "bad input.\nDo you want to change the threshold for" + \
                    " the network representation ? (y/n) : ")
        if change_threshold == "y":
            threshold = max(matrix_mi.max()) + 1
            while threshold >= max(matrix_mi.max()) \
               or threshold <= min(matrix_mi.max()):
                while True:
                    try:
                        threshold = int(input(\
    "Value of the threshold: (max value = {:.1f} ; min value = {:.1f}):"\
    .format(max(matrix_mi.max()), min(matrix_mi.min()))))
                        break
                    except ValueError:
                        print("This is not a valid number. Try again.")
        else:
            threshold = np.max(np.max(matrix_mi))/2
        mat_visual.show_pymol_network(args.py, threshold)


#### MAIN ####
if __name__ == "__main__":
    analysis_mutual_information()
