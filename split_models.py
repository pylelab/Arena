#!/usr/bin/python

#Reference: Zion R Perry, Anna Marie Pyle, Chengxin Zhang (2023) Arena: rapid and accurate reconstruction of full atomic RNA structures from coarse-grained models.


import sys,os

def split_models(filename):

	'''
	Purpose: Split a PDB file into separate models
	Input: PDB file with multiple models
	Output: Folder with a PDB file for each model

	Usage: python split_models.py input.pdb
	'''

	cwd = os.getcwd()
	directory = filename.strip(".pdb") + "_models"
	path = os.path.join(cwd, directory)
	os.mkdir(path)

	# open input file in read-only mode
	multiple_models = open(filename, "r")
	model_number = 1
	# return as a single string and split by keyword (faster than by line)
	txt = multiple_models.read()
	for block in txt.split("\nMODEL"):
		# check if ATOM is in the block
		if "ATOM " in block:
			model_filename = "model_{}.pdb".format(model_number)
			# create file for the first model in write mode
			path_to_file = path + "/" + model_filename
			model_file = open(path_to_file, "w")
			# write block to file
			model_file.write(block)
			model_file.close()
			model_number += 1

	multiple_models.close()

#allow command line inputs
if __name__=="__main__":
	if len(sys.argv) < 2:
		print(split_models.__doc__)
	if len(sys.argv) == 2:	
		split_models(sys.argv[1])
	else:	
		print("Incorrect input. Run 'python split_models.py' for help.")



