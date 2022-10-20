const char* docstring=""
"MissingRNAatoms input.pdb output.pdb [option]\n"
"    Find and fill missing atoms in input.pdb\n"
"    Output filled model to output.pdb\n"
"option:\n"
"    0 - only check missing atoms\n"
"    1 - remove non-standard atoms\n"
"    2 - remove non-standard atoms and residues with missing atoms\n"
"    3 - remove non-standard atoms, fill atoms for residues with at\n"
"        least three atoms\n"
"    4 - remove non-standard atoms, fill atoms for residues with at\n"
"        least three atoms, remove residues with less than three atoms\n"
"    5 - (default) remove non-standard atoms, fill atoms for residues with at\n"
"        least one atom\n"
;

#include <iostream>
#include <iomanip>
#include "MissingRNAatom.hpp"
#include "BondLengths.hpp"
#include "BaseConformation.hpp"
#include "cssr.hpp"
#include "AtomicClashes.hpp"

int main(int argc,char **argv)
{
    if (argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string infile =argv[1];
    string outfile=argv[2];
    int option    =(argc>3)?atoi(argv[3]):5;
    double tolerance =(argc>4)?atoi(argv[4]):0.1; // default tolerance is 1%
    
    int atomic_detail=2;
    int allowX=0;
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);
    const ModelUnit pep;
    
    vector<vector<size_t> >res_str_vec;
    vector<pair<double,vector<size_t> > > bp_vec;
    if (option>=5)
    {
        cssr(pdb_entry, res_str_vec, bp_vec);
        filter_bp(bp_vec);
    }
   
    map<string, map<string,vector<double> > >ideal_rna=parse_ideal_rna();
    MissingRNAatom(pdb_entry,ideal_rna,bp_vec,option);

    if (option>=5)
    {
        standardize_pdb_order(pdb_entry);

        for (int c=0; c<pep.chains.size(); c++){
            int length = pep.chains[c].residues.size();
            int iterations = round(0.04*length+100);

            for (int t=0; t<iterations; t++){
        	    int moved = fix_bond_lengths(pdb_entry, tolerance);
                
                if (t%10==0)
                {
                    size_t bp;
                    for (bp=0;bp<res_str_vec.size();bp++) res_str_vec[bp].clear();
                    res_str_vec.clear();
                    for (bp=0;bp<bp_vec.size();bp++) bp_vec[bp].second.clear();
                    bp_vec.clear();

                    cssr(pdb_entry, res_str_vec, bp_vec);
                    filter_bp(bp_vec);
                }
                moved+=fixBaseConformation(pdb_entry, ideal_rna, bp_vec);
                moved+=fix_clashes(pdb_entry);
     
                // if no atoms are moved, immediately break out of this for loop
        	    if (moved==0) break;
                cout<<"t="<<t<<" moved="<<moved<<endl;
            }
    	}
    }
    
    map<string, map<string,vector<double> > >().swap(ideal_rna);
    write_pdb_structure(outfile.c_str(),pdb_entry);
    vector<ChainUnit>().swap(pdb_entry.chains);
    vector<vector<size_t> >().swap(res_str_vec);
    vector<pair<double,vector<size_t> > > ().swap(bp_vec);
    return 0;
}
