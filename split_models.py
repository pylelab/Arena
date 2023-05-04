#!/usr/bin/python3

#Reference: Zion R Perry, Anna Marie Pyle, Chengxin Zhang (2023) Arena: rapid and accurate reconstruction of full atomic RNA structures from coarse-grained models.


import sys,os
import shutil

def split_models(filename, outfile):

    '''
    Usage: python3 split_models.py input.pdb
    Purpose: Split a PDB file into separate models
    Input: PDB file with multiple models
    Output: Folder with a PDB file for each model

    Usage: python3 split_models.py input.pdb output.pdb
    Purpose: Split a PDB file into separate models, run Arena on individual
             models, and combine the Arena output into one multi-model file.
    Output: A PDB file with full-atomic model for all models.

    '''

    cwd = os.getcwd()
    directory = filename.strip(".pdb") + "_models"
    path = os.path.join(cwd, directory)
    if not os.path.isdir(path):
        os.makedirs(path)
    else:
        print("WARNING! "+path+" already exist")

    # open input file in read-only mode
    multiple_models = open(filename, "r")
    model_number = 1
    # return as a single string and split by keyword (faster than by line)
    txt = multiple_models.read()
    filename_list=[]
    for block in txt.split("\nMODEL"):
        # check if ATOM is in the block
        if "ATOM " in block:
            model_filename = "model_{}.pdb".format(model_number)
            # create file for the first model in write mode
            path_to_file = path + "/" + model_filename
            model_file = open(path_to_file, "w")
            # write block to file
            if block.startswith(' ') and block[8] in "1234567890":
                block="MODEL"+block
            model_file.write(block)
            model_file.close()
            model_number += 1
            filename_list.append(path_to_file)

    multiple_models.close()

    if outfile=='':
        return

    Arena=os.path.join(os.path.dirname(os.path.abspath(__file__)),"Arena")
    if os.name=="nt":
        Arena+=".ext"
    if not os.path.isfile(Arena):
        print("ERRPR! Cannot locate the 'Arena' executable at "+Arena)
        return

    txt=''
    for path_to_file in filename_list:
        print("parsing "+path_to_file)
        fp=open(path_to_file)
        for line in fp:
            if line.startswith("MODEL "):
                txt+=line
                break
        fp.close()
        cmd=Arena+' '+path_to_file+' '+path_to_file+".full.pdb"
        os.system(cmd)
        fp=open(path_to_file+".full.pdb")
        for line in fp.read().splitlines():
            if line.startswith("ATOM  ") or line.startswith("TER"):
                txt+=line+'\n'
        fp.close()
        txt+="ENDMDL\n"
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    if os.path.isfile(outfile) and os.path.isdir(path):
        shutil.rmtree(path)
    return

#allow command line inputs
if __name__=="__main__":
    if len(sys.argv) < 2:
        print(split_models.__doc__)
    if len(sys.argv) in [2,3]:
        filename = sys.argv[1]
        outfile  = ''
        if len(sys.argv)>2:
            outfile=sys.argv[2]
        filename_list=split_models(filename, outfile)
    else:    
        print("Incorrect input. Run 'python3 split_models.py' for help.")
