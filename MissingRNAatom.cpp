const char* docstring=""
"MissingRNAatoms input.pdb output.pdb [option]\n"
"    Find and fill missing atoms in input.pdb\n"
"    Output filled model to output.pdb\n"
"option:\n"
"    0 - (default) only check missing atoms\n"
"    1 - remove non-standard atoms\n"
"    2 - remove non-standard atoms and residues with missing atoms\n"
"    3 - remove non-standard atoms, fill atoms for residues with at\n"
"        least three atoms\n"
"    4 - remove non-standard atoms, fill atoms for residues with at\n"
"        least three atoms, remove residues with less than three atoms\n"
"    5 - remove non-standard atoms, fill atoms for residues with at\n"
"        least one atom\n"
;

#include <iostream>
#include "MissingRNAatom.hpp"

int main(int argc,char **argv)
{
    if (argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string infile =argv[1];
    string outfile=argv[2];
    int option    =(argc>3)?atoi(argv[3]):0;
    
    int atomic_detail=2;
    int allowX=0;
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);
    map<string, map<string,vector<float> > >ideal_rna=parse_ideal_rna();
    MissingRNAatom(pdb_entry,ideal_rna,option);
    map<string, map<string,vector<float> > >().swap(ideal_rna);
    standardize_pdb_order(pdb_entry);
    write_pdb_structure(outfile.c_str(),pdb_entry);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 0;
}
