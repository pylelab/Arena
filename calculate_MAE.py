#!/usr/bin/python3

'''
Calculate the MAE between an input PDB compared to a reference PDB.
Only include atoms present in both structures.
How to run: python calculate_RMSD.py -r reference.pdb -i input.pdb
'''

import argparse
parser = argparse.ArgumentParser(description= __doc__)
parser.add_argument("-r",type=str,help="reference PDB")
parser.add_argument("-i",type=str,help="input PDB")
args = parser.parse_args()

# loop over atoms in reference PDB
ref_information = {}
ref_pdb = open(args.r)
for line in ref_pdb.read().splitlines():
	if line[0:6] == "ATOM  ":
		atom_name = line[12:16]
		residue_number = line[22:26]
		x = float(line[30:38])
		y = float(line[38:46])
		z = float(line[46:54])
		ref_information[atom_name, residue_number] = [x, y, z]

ref_pdb.close()

# loop over atoms in input PDB
in_information = {}
in_pdb = open(args.i)
for line in in_pdb.read().splitlines():
	if line[0:6] == "ATOM  ":
		atom_name = line[12:16]
		residue_number = line[22:26]
		x = float(line[30:38])
		y = float(line[38:46])
		z = float(line[46:54])
		in_information[atom_name, residue_number] = [x, y, z]

in_pdb.close()

# calculate MAE
count = 0
error = 0
# loop through keys in input
for key in in_information.keys():
	# extract the x,y,z coordinates
	ref_atom_coord = ref_information.get(key) # key may not be in reference so use 'get'
	in_atom_coord = in_information[key]
	# only calculate RMSD for those atoms present in the reference
	if ref_atom_coord != None:
		# calculate MAE
		count += 1
		error += (abs(ref_atom_coord[0]-in_atom_coord[0])+abs(ref_atom_coord[1]-in_atom_coord[1])+abs(ref_atom_coord[2]-in_atom_coord[2]))/3

mae = error/count

print(mae)
