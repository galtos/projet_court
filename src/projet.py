import sys

file_pdb = sys.argv[1]
print(sys.argv[1])

with open(file_pdb,"r") as file_pdb:
    res_count = 0
    for line in file_pdb:
        if line[0:6].strip() == "ATOM":
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_num = int(line[22:26])
            if atom_name == "CA":
                res_count += 1
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                print(res_name, res_num, x, y, z)
"""
Extract from a pdb file the xyz coordonées 
"""
#### MAIN ####

