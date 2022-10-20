#ifndef MissingRNAatom_HPP
#define MissingRNAatom_HPP 1

#include <iostream>
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "IdealRNA.hpp"

/* return true if residue is standard nucleotide
 * return false otherwise */
bool checkMissingRNAatom(ResidueUnit &residue, 
    map<string, map<string,vector<double> > >&ideal_rna,
    vector<string> &missing_atom_vec,
    vector<size_t> &nonstd_atom_vec, const int option=0)
{
    missing_atom_vec.clear();
    nonstd_atom_vec.clear();
    if (ideal_rna.count(residue.resn)==0)
    {
        cout<<residue.resn<<residue.resi<<residue.icode
            <<"\tnon-standard residue"<<endl;
        return false;
    }
    size_t a,i,j;
    for (a=0;a<residue.atoms.size();a++)
        if (ideal_rna[residue.resn].count(residue.atoms[a].name)==0)
            nonstd_atom_vec.push_back(a);
    if (nonstd_atom_vec.size())
    {
        cout<<residue.resn<<residue.resi<<residue.icode<<"\t"
            <<nonstd_atom_vec.size()<<" non-standard atoms\t";
        for (i=0;i<nonstd_atom_vec.size();i++)
            cout<<residue.atoms[nonstd_atom_vec[i]].name<<' ';
        cout<<endl;
        if (option)
        {
            for (i=0;i<nonstd_atom_vec.size();i++)
                residue.atoms[nonstd_atom_vec[i]].name.clear();
            for (i=0;i<residue.atoms.size();i++)
            {
                if (residue.atoms[i].name.size()) continue;
                for (j=i+1;j<residue.atoms.size();j++)
                {
                    if (residue.atoms[j].name.size()==0) continue;
                    residue.atoms[i].name  =residue.atoms[j].name;
                    residue.atoms[i].xyz[0]=residue.atoms[j].xyz[0];
                    residue.atoms[i].xyz[1]=residue.atoms[j].xyz[1];
                    residue.atoms[i].xyz[2]=residue.atoms[j].xyz[2];
                    residue.atoms[j].name.clear();
                    break;
                }
            }
            for (i=0;i<nonstd_atom_vec.size();i++) residue.atoms.pop_back();
        }
    }
    
    map<string,vector<double> >::iterator it;
    bool found;
    for (it=ideal_rna[residue.resn].begin(); it!=ideal_rna[residue.resn].end(); it++)
    {
        found=false;
        for (a=0;a<residue.atoms.size();a++)
        {
            if (residue.atoms[a].name==it->first)
            {
                found=true;
                break;
            }
        }
        if (found==false) missing_atom_vec.push_back(it->first);
    }
    // terminal output that lists missing atoms
    if (missing_atom_vec.size())
    {
        cout<<residue.resn<<residue.resi<<residue.icode<<"\t"
            <<missing_atom_vec.size()<<" missing atoms\t";
        for (a=0;a<missing_atom_vec.size();a++)
            cout<<missing_atom_vec[a]<<' ';
        cout<<endl;
    }
    return true;
}

size_t countUniqAtom(vector<vector<double> >&xyz_list2)
{
    size_t atomNum=xyz_list2.size();
    double Extra=1.0e-6;
    size_t dup;
    size_t a,b;
    for (a=1;a<xyz_list2.size();a++)
    {
        dup=0;
        for (b=0;b<a;b++)
            if (fabs(xyz_list2[a][0]-xyz_list2[b][0])+
                fabs(xyz_list2[a][1]-xyz_list2[b][1])+
                fabs(xyz_list2[a][2]-xyz_list2[b][2])<Extra) dup++;
        if (dup) atomNum--;
    }
    return atomNum;
}

void fillMissingRNAatom(ResidueUnit &residue,
    map<string, map<string,vector<double> > >&ideal_rna,
    vector<string>missing_atom_vec)
{
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1(residue.atoms.size(),tmp);
    vector<vector<double> > xyz_list2(residue.atoms.size(),tmp);
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t

    size_t a,b;
    for (a=0;a<residue.atoms.size();a++)
    {
        xyz_list1[a][0]=ideal_rna[residue.resn][residue.atoms[a].name][0];
        xyz_list1[a][1]=ideal_rna[residue.resn][residue.atoms[a].name][1];
        xyz_list1[a][2]=ideal_rna[residue.resn][residue.atoms[a].name][2];

        xyz_list2[a][0]=residue.atoms[a].xyz[0];
        xyz_list2[a][1]=residue.atoms[a].xyz[1];
        xyz_list2[a][2]=residue.atoms[a].xyz[2];
    }

    size_t atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size())
        cerr<<"duplicated coordinate in "
            <<residue.resn<<residue.resi<<residue.icode<<endl;

    //if (residue.atoms.size()>=3)
    if (atomNum>=3)
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    else
    {
        cerr<<"too few ("<<atomNum<<") unique atoms to fill "
            <<residue.resn<<residue.resi<<residue.icode<<endl;
        RotMatix.push_back(tmp);
        RotMatix.push_back(tmp);
        RotMatix.push_back(tmp);
        RotMatix[0][0]=RotMatix[1][1]=RotMatix[2][2]=1;
        
        TranVect=tmp;
        for (a=0;a<xyz_list1.size();a++)
        {
            TranVect[0]+=xyz_list2[a][0]-xyz_list1[a][0];
            TranVect[1]+=xyz_list2[a][1]-xyz_list1[a][1];
            TranVect[2]+=xyz_list2[a][2]-xyz_list1[a][2];
        }
        TranVect[0]/=xyz_list1.size();
        TranVect[1]/=xyz_list1.size();
        TranVect[2]/=xyz_list1.size();
    }
    
    for (a=0;a<residue.atoms.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
                
    AtomUnit atom;
    atom.xyz.assign(3,0);
    atom.movable=1; // set movable to 1 if the atom is added

    for (a=0;a<missing_atom_vec.size();a++)
    {
        atom.name=missing_atom_vec[a];
        ChangeCoor(ideal_rna[residue.resn][atom.name],
            RotMatix, TranVect, atom.xyz);
        residue.atoms.push_back(atom);
    }

    atom.name.clear();
    atom.xyz.clear();
    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return;
}

void fillMissingRNAatom(size_t r, ChainUnit &chain, ChainUnit &fill_chain,
    map<string, map<string,vector<double> > >&ideal_rna,
    vector<string>missing_atom_vec)
{
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t

    size_t a,a0,a1,a2;
    string key,key0,key1,key2;
    size_t r0,r1,r2;
    size_t prevAtomNum=0;
    if (r==0) 
    {
        r0=0;
        r1=1;
        r2=2;
    }
    else if (r==chain.residues.size()-1)
    {
        r0=r-2;
        r1=r-1;
        r2=r;
        prevAtomNum=chain.residues[r0].atoms.size()+chain.residues[r1].atoms.size();
    }
    else
    {
        r0=r-1;
        r1=r;
        r2=r+1;
        prevAtomNum=chain.residues[r0].atoms.size();
    }
    key=chain.residues[r0].resn+chain.residues[r1].resn+chain.residues[r2].resn;
    key0=key+'0';
    key1=key+'1';
    key2=key+'2';
    if (r==0)                            key=key0;
    else if (r==chain.residues.size()-1) key=key2;
    else                                 key=key1;
    xyz_list1.assign(chain.residues[r0].atoms.size()+
                     chain.residues[r1].atoms.size()+
                     chain.residues[r2].atoms.size(),tmp);
    xyz_list2.assign(xyz_list1.size(),tmp);

    a=0;
    for (a0=0;a0<chain.residues[r0].atoms.size();a0++)
    {
        xyz_list1[a][0]=ideal_rna[key0][chain.residues[r0].atoms[a0].name][0];
        xyz_list1[a][1]=ideal_rna[key0][chain.residues[r0].atoms[a0].name][1];
        xyz_list1[a][2]=ideal_rna[key0][chain.residues[r0].atoms[a0].name][2];

        xyz_list2[a][0]=chain.residues[r0].atoms[a0].xyz[0];
        xyz_list2[a][1]=chain.residues[r0].atoms[a0].xyz[1];
        xyz_list2[a][2]=chain.residues[r0].atoms[a0].xyz[2];
        a++;
    }
    for (a1=0;a1<chain.residues[r1].atoms.size();a1++)
    {
        xyz_list1[a][0]=ideal_rna[key1][chain.residues[r1].atoms[a1].name][0];
        xyz_list1[a][1]=ideal_rna[key1][chain.residues[r1].atoms[a1].name][1];
        xyz_list1[a][2]=ideal_rna[key1][chain.residues[r1].atoms[a1].name][2];

        xyz_list2[a][0]=chain.residues[r1].atoms[a1].xyz[0];
        xyz_list2[a][1]=chain.residues[r1].atoms[a1].xyz[1];
        xyz_list2[a][2]=chain.residues[r1].atoms[a1].xyz[2];
        a++;
    }
    for (a2=0;a2<chain.residues[r2].atoms.size();a2++)
    {
        xyz_list1[a][0]=ideal_rna[key2][chain.residues[r2].atoms[a2].name][0];
        xyz_list1[a][1]=ideal_rna[key2][chain.residues[r2].atoms[a2].name][1];
        xyz_list1[a][2]=ideal_rna[key2][chain.residues[r2].atoms[a2].name][2];

        xyz_list2[a][0]=chain.residues[r2].atoms[a2].xyz[0];
        xyz_list2[a][1]=chain.residues[r2].atoms[a2].xyz[1];
        xyz_list2[a][2]=chain.residues[r2].atoms[a2].xyz[2];
        a++;
    }

    size_t atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size())
        cerr<<"duplicated coordinate when reconstructing "
            <<chain.residues[r].resn<<chain.residues[r].resi
            <<chain.residues[r].icode<<endl;

    if (atomNum>=3)
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    else
    {
        cerr<<"too few ("<<atomNum<<") unique atoms to reconstruct "
            <<chain.residues[r].resn<<chain.residues[r].resi
            <<chain.residues[r].icode<<endl;
        RotMatix.push_back(tmp);
        RotMatix.push_back(tmp);
        RotMatix.push_back(tmp);
        RotMatix[0][0]=RotMatix[1][1]=RotMatix[2][2]=1;
        
        TranVect=tmp;
        for (a=0;a<xyz_list1.size();a++)
        {
            TranVect[0]+=xyz_list2[a][0]-xyz_list1[a][0];
            TranVect[1]+=xyz_list2[a][1]-xyz_list1[a][1];
            TranVect[2]+=xyz_list2[a][2]-xyz_list1[a][2];
        }
        TranVect[0]/=xyz_list1.size();
        TranVect[1]/=xyz_list1.size();
        TranVect[2]/=xyz_list1.size();
    }
    
    /* adjust TranVect */
    for (a=0;a<chain.residues[r].atoms.size();a++)
    {
        ChangeCoor(ideal_rna[key][chain.residues[r].atoms[a].name],
            RotMatix,TranVect,xyz_list1[a+prevAtomNum]);
        tmp[0]+=(chain.residues[r].atoms[a].xyz[0]-xyz_list1[a+prevAtomNum][0]);
        tmp[1]+=(chain.residues[r].atoms[a].xyz[1]-xyz_list1[a+prevAtomNum][1]);
        tmp[2]+=(chain.residues[r].atoms[a].xyz[2]-xyz_list1[a+prevAtomNum][2]);
    }
    TranVect[0]+=tmp[0]/chain.residues[r].atoms.size();
    TranVect[1]+=tmp[1]/chain.residues[r].atoms.size();
    TranVect[2]+=tmp[2]/chain.residues[r].atoms.size();
                
    /* eventually assign the coordinate */
    for (a=0;a<missing_atom_vec.size();a++)
        ChangeCoor(ideal_rna[key][missing_atom_vec[a]],
            RotMatix, TranVect, fill_chain.residues[r].atoms[a].xyz);
    
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();

    key.clear();
    key0.clear();
    key1.clear();
    key2.clear();
    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return;
}

bool fillMissingRNAatomInPair(const size_t c1, const size_t r1, 
    const size_t c2, const size_t r2,
    const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    map<string, map<string,vector<double> > >&ideal_rna,
    vector<string>missing_atom_vec)
{
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    string key1=pdb_entry.chains[c1].residues[r1].resn;
    string key2=pdb_entry.chains[c2].residues[r2].resn;

    size_t a,a1,a2;
    size_t atomNum=0;
    for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
        atomNum+=pdb_entry.chains[c1].residues[r1].atoms[a1].movable==0;
    for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
        atomNum+=pdb_entry.chains[c2].residues[r2].atoms[a2].movable==0;
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);

    a=0;
    for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
    {
        if (pdb_entry.chains[c1].residues[r1].atoms[a1].movable) continue;
        xyz_list1[a][0]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a1].name][0];
        xyz_list1[a][1]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a1].name][1];
        xyz_list1[a][2]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a1].name][2];

        xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
        xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
        xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];
        a++;
    }
    for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
    {
        if (pdb_entry.chains[c2].residues[r2].atoms[a2].movable) continue;
        xyz_list1[a][0]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a2].name][0];
        xyz_list1[a][1]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a2].name][1];
        xyz_list1[a][2]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a2].name][2];

        xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[0];
        xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[1];
        xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[2];
        a++;
    }

    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size())
        cerr<<"duplicated coordinate when reconstructing "
            <<pdb_entry.chains[c1].residues[r1].resn
            <<pdb_entry.chains[c1].residues[r1].resi
            <<pdb_entry.chains[c1].residues[r1].icode<<" and "
            <<pdb_entry.chains[c2].residues[r2].resn
            <<pdb_entry.chains[c2].residues[r2].resi
            <<pdb_entry.chains[c2].residues[r2].icode<<endl;

    if (atomNum<3)
    {
        for (a=0;a<xyz_list1.size();a++)
        {
            xyz_list1[a].clear();
            xyz_list2[a].clear();
        }
        xyz_list1.clear();
        xyz_list2.clear();
        key1.clear();
        key2.clear();
        return false;
    }
    
    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    
    /* eventually assign the coordinate */
    for (a=0;a<missing_atom_vec.size();a++)
        ChangeCoor(ideal_rna[key1][missing_atom_vec[a]],
            RotMatix, TranVect, fill_chain.residues[r1].atoms[a].xyz);
    
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    key1.clear();
    key2.clear();
    return true;
}

bool fillMissingRNAatomInStack(const size_t c1, const size_t r1,
    const size_t c2, const size_t r2, const bool stack_prev,
    const bool stack_next, const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    map<string, map<string,vector<double> > >&ideal_rna,
    vector<string>missing_atom_vec)
{
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    string key, key1, key2;
    
    cout<<"reconstructing residue pair "
        <<pdb_entry.chains[c1].residues[r1].resn
        <<pdb_entry.chains[c1].residues[r1].resi
        <<pdb_entry.chains[c1].residues[r1].icode<<" and "
        <<pdb_entry.chains[c2].residues[r2].resn
        <<pdb_entry.chains[c2].residues[r2].resi
        <<pdb_entry.chains[c2].residues[r2].icode<<endl;

    size_t a,a1,a2;
    size_t atomNum=0;
    for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        atomNum+=pdb_entry.chains[c1].residues[r1].atoms[a].movable==0;
    for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        atomNum+=pdb_entry.chains[c2].residues[r2].atoms[a].movable==0;
    if (stack_prev)
    {
        for (a=0;a<pdb_entry.chains[c1].residues[r1-1].atoms.size();a++)
            atomNum+=pdb_entry.chains[c1].residues[r1-1].atoms[a].movable==0;
        for (a=0;a<pdb_entry.chains[c2].residues[r2+1].atoms.size();a++)
            atomNum+=pdb_entry.chains[c2].residues[r2+1].atoms[a].movable==0;
    }
    if (stack_next)
    {
        for (a=0;a<pdb_entry.chains[c1].residues[r1+1].atoms.size();a++)
            atomNum+=pdb_entry.chains[c1].residues[r1+1].atoms[a].movable==0;
        for (a=0;a<pdb_entry.chains[c2].residues[r2-1].atoms.size();a++)
            atomNum+=pdb_entry.chains[c2].residues[r2-1].atoms[a].movable==0;
    }
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);

    a=0;
    if (stack_prev && stack_next)
    {
        key=pdb_entry.chains[c1].residues[r1-1].resn+
            pdb_entry.chains[c1].residues[r1  ].resn+
            pdb_entry.chains[c1].residues[r1+1].resn;
        key1=key+"1";
        key2=key+"4";
    }
    else if (stack_prev)
    {
        key=pdb_entry.chains[c1].residues[r1-1].resn+
            pdb_entry.chains[c1].residues[r1  ].resn;
        key1=key+"1";
        key2=key+"2";
    }
    else if (stack_next)
    {
        key=pdb_entry.chains[c1].residues[r1  ].resn+
            pdb_entry.chains[c1].residues[r1+1].resn;
        key1=key+"0";
        key2=key+"3";
    }
    //cout<<"key1="<<key1<<". key2="<<key2<<"."<<endl;
    for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
    {
        if (pdb_entry.chains[c1].residues[r1].atoms[a1].movable) continue;
        xyz_list1[a][0]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a1].name][0];
        xyz_list1[a][1]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a1].name][1];
        xyz_list1[a][2]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a1].name][2];

        xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
        xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
        xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];
        a++;
    }
    for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
    {
        if (pdb_entry.chains[c2].residues[r2].atoms[a2].movable) continue;
        xyz_list1[a][0]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a2].name][0];
        xyz_list1[a][1]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a2].name][1];
        xyz_list1[a][2]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a2].name][2];

        xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[0];
        xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[1];
        xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[2];
        a++;
    }
    if (stack_prev)
    {
        key1=key+"0";
        key2=key+"3";
        if (stack_next) key2=key+"5";
        //cout<<"key1="<<key1<<". key2="<<key2<<"."<<endl;
        for (a1=0;a1<pdb_entry.chains[c1].residues[r1-1].atoms.size();a1++)
        {
            if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].movable) continue;
            xyz_list1[a][0]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1-1].atoms[a1].name][0];
            xyz_list1[a][1]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1-1].atoms[a1].name][1];
            xyz_list1[a][2]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1-1].atoms[a1].name][2];

            xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[2];
            a++;
        }
        for (a2=0;a2<pdb_entry.chains[c2].residues[r2+1].atoms.size();a2++)
        {
            if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].movable) continue;
            xyz_list1[a][0]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2+1].atoms[a2].name][0];
            xyz_list1[a][1]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2+1].atoms[a2].name][1];
            xyz_list1[a][2]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2+1].atoms[a2].name][2];

            xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz[2];
            a++;
        }
    }
    if (stack_next)
    {
        key1=key+"1";
        key2=key+"2";
        if (stack_prev)
        {
            key1=key+"2";
            key2=key+"3";
        }
        //cout<<"key1="<<key1<<". key2="<<key2<<"."<<endl;
        for (a1=0;a1<pdb_entry.chains[c1].residues[r1+1].atoms.size();a1++)
        {
            if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].movable) continue;
            xyz_list1[a][0]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1+1].atoms[a1].name][0];
            xyz_list1[a][1]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1+1].atoms[a1].name][1];
            xyz_list1[a][2]=ideal_rna[key1][pdb_entry.chains[c1].residues[r1+1].atoms[a1].name][2];

            xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[2];
            a++;
        }
        for (a2=0;a2<pdb_entry.chains[c2].residues[r2-1].atoms.size();a2++)
        {
            if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].movable) continue;
            xyz_list1[a][0]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2-1].atoms[a2].name][0];
            xyz_list1[a][1]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2-1].atoms[a2].name][1];
            xyz_list1[a][2]=ideal_rna[key2][pdb_entry.chains[c2].residues[r2-1].atoms[a2].name][2];

            xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz[2];
            a++;
        }
    }
    if (stack_prev && stack_next)
    {
        key1=key+"1";
        key2=key+"4";
    }
    else if (stack_prev)
    {
        key1=key+"1";
        key2=key+"2";
    }
    else if (stack_next)
    {
        key1=key+"0";
        key2=key+"3";
    }

    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size())
        cerr<<"duplicated coordinate when reconstructing "
            <<pdb_entry.chains[c1].residues[r1].resn
            <<pdb_entry.chains[c1].residues[r1].resi
            <<pdb_entry.chains[c1].residues[r1].icode<<" and "
            <<pdb_entry.chains[c2].residues[r2].resn
            <<pdb_entry.chains[c2].residues[r2].resi
            <<pdb_entry.chains[c2].residues[r2].icode<<endl;

    if (atomNum<3)
    {
        for (a=0;a<xyz_list1.size();a++)
        {
            xyz_list1[a].clear();
            xyz_list2[a].clear();
        }
        xyz_list1.clear();
        xyz_list2.clear();
        key.clear();
        key1.clear();
        key2.clear();
        return false;
    }
    
    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    
    /* adjust TranVect */
    atomNum=0;
    for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
    {
        if (pdb_entry.chains[c1].residues[r1].atoms[a].movable) continue;
        ChangeCoor(ideal_rna[key1][pdb_entry.chains[c1].residues[r1].atoms[a].name],
            RotMatix,TranVect,xyz_list1[atomNum]);
        tmp[0]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
        tmp[1]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
        tmp[2]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
        atomNum++;
    }
    for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
    {
        if (pdb_entry.chains[c2].residues[r2].atoms[a].movable) continue;
        ChangeCoor(ideal_rna[key2][pdb_entry.chains[c2].residues[r2].atoms[a].name],
            RotMatix,TranVect,xyz_list1[atomNum]);
        tmp[0]+=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
        tmp[1]+=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
        tmp[2]+=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
        atomNum++;
    }
    TranVect[0]+=tmp[0]/atomNum;
    TranVect[1]+=tmp[1]/atomNum;
    TranVect[2]+=tmp[2]/atomNum;
    
    /* eventually assign the coordinate */
    for (a=0;a<missing_atom_vec.size();a++)
        ChangeCoor(ideal_rna[key1][missing_atom_vec[a]],
            RotMatix, TranVect, fill_chain.residues[r1].atoms[a].xyz);
    
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    key.clear();
    key1.clear();
    key2.clear();
    return true;
}

bool check3atomInPair(ModelUnit &pdb_entry, const size_t c2,const size_t r2,
    const size_t c1, const size_t r1)
{
    size_t atomNum=0;
    size_t a=0;
    for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        atomNum+=(pdb_entry.chains[c1].residues[r1].atoms[a].movable==0);
    for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        atomNum+=(pdb_entry.chains[c2].residues[r2].atoms[a].movable==0);
    return atomNum>=3;
}

bool checkAtomInStack(const vector<pair<double, vector<size_t> > >&bp_vec,
    const size_t c1, const size_t r1, const size_t c2, const size_t r2,
    bool &stack_prev, bool &stack_next)
{
    size_t bp;
    stack_prev=false;
    stack_next=false;
    for (bp=0;bp<bp_vec.size();bp++)
    {
        if (bp_vec[bp].second[4]) continue; // Watson-Crick pair only
        if (c1!=bp_vec[bp].second[0]   || c2!=bp_vec[bp].second[2]) continue;
        if (r1==bp_vec[bp].second[1]-1 && r2==bp_vec[bp].second[3]+1)
            stack_next=true;
        if (r1==bp_vec[bp].second[1]+1 && r2==bp_vec[bp].second[3]-1)
            stack_prev=true;
    }
    return (stack_prev || stack_next);
}

bool MissingRNAatom(ModelUnit &pdb_entry, const size_t c,
    map<string, map<string,vector<double> > >&ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec,
    const int option=0)
{
    vector<string> missing_atom_vec;
    vector<size_t> nonstd_atom_vec;
    size_t r,i,j;
    vector<size_t> remove_residue_vec;
    map<size_t,vector<string> >fill_residue_map;
    for (r=0;r<pdb_entry.chains[c].residues.size();r++)
    {
        checkMissingRNAatom(pdb_entry.chains[c].residues[r], ideal_rna,
            missing_atom_vec, nonstd_atom_vec, option);
        if (option==2 && missing_atom_vec.size())
        {
            cout<<"remove residue "
                <<pdb_entry.chains[c].residues[r].resn
                <<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode
                <<" with "
                <<missing_atom_vec.size()
                <<" missing atoms"<<endl;
            pdb_entry.chains[c].residues[r].resn.clear();
            remove_residue_vec.push_back(r);
        }
        if (option>=3 && missing_atom_vec.size() && (
            pdb_entry.chains[c].residues[r].atoms.size()>=3
            || (option==5 && pdb_entry.chains[c].residues.size()<=2)))
        {
            cout<<"fill residue "
                <<pdb_entry.chains[c].residues[r].resn
                <<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode
                <<" with "
                <<missing_atom_vec.size()
                <<" missing atoms"<<endl;
            fillMissingRNAatom(pdb_entry.chains[c].residues[r], ideal_rna, missing_atom_vec);
            missing_atom_vec.clear();
        }
        else if (option==3 && pdb_entry.chains[c].residues[r].atoms.size()<3)
        {
            cout<<"failed to fill residue "
                <<pdb_entry.chains[c].residues[r].resn
                <<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode
                <<" with "
                <<pdb_entry.chains[c].residues[r].atoms.size()
                <<" atoms and "
                <<missing_atom_vec.size()
                <<" missing atoms"<<endl;
        }
        else if (option==4 && pdb_entry.chains[c].residues[r].atoms.size()<3)
        {
            cout<<"remove residue "
                <<pdb_entry.chains[c].residues[r].resn
                <<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode
                <<" with "
                <<pdb_entry.chains[c].residues[r].atoms.size()
                <<" atoms"<<endl;
            pdb_entry.chains[c].residues[r].resn.clear();
            remove_residue_vec.push_back(r);
        }
        else if (option==5)
        {
            if (pdb_entry.chains[c].residues[r].atoms.size()==0)
            {
                cout<<"remove residue "
                    <<pdb_entry.chains[c].residues[r].resn
                    <<pdb_entry.chains[c].residues[r].resi
                    <<pdb_entry.chains[c].residues[r].icode
                    <<" without standard atom"<<endl;
                pdb_entry.chains[c].residues[r].resn.clear();
                remove_residue_vec.push_back(r);
            }
            else if (missing_atom_vec.size())
                fill_residue_map[r]=missing_atom_vec;
        }
    }

    /* remove empty residues */
    if (remove_residue_vec.size())
    {
        for (i=0;i<pdb_entry.chains[c].residues.size();i++)
        {
            if (pdb_entry.chains[c].residues[i].resn.size()) continue;
            for (j=i+1;j<pdb_entry.chains[c].residues.size();j++)
            {
                if (pdb_entry.chains[c].residues[j].resn.size()==0) continue;
                pdb_entry.chains[c].residues[i].het  =pdb_entry.chains[c].residues[j].het;
                pdb_entry.chains[c].residues[i].resi =pdb_entry.chains[c].residues[j].resi;
                pdb_entry.chains[c].residues[i].icode=pdb_entry.chains[c].residues[j].icode;
                pdb_entry.chains[c].residues[i].resn =pdb_entry.chains[c].residues[j].resn;
                pdb_entry.chains[c].residues[i].atoms=pdb_entry.chains[c].residues[j].atoms;
                pdb_entry.chains[c].residues[j].resn.clear();
                break;
            }
        }
        for (i=0;i<remove_residue_vec.size();i++) pdb_entry.chains[c].residues.pop_back();
    }

    /* fill residues with 1 or 2 atoms for chains with >=3 residues */
    if (option==5 && fill_residue_map.size())
    {
        /* declare memory */
        ChainUnit fill_chain;
        ResidueUnit fill_residue;
        AtomUnit fill_atom;
        size_t a;
        
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            fill_chain.residues.push_back(fill_residue);
            if (fill_residue_map.count(r)==0) continue;
            for (a=0;a<fill_residue_map[r].size();a++)
            {
                fill_chain.residues[r].atoms.push_back(fill_atom);
                fill_chain.residues[r].atoms[a].name=fill_residue_map[r][a];
                fill_chain.residues[r].atoms[a].xyz.assign(3,0);
                fill_chain.residues[r].atoms[a].movable=1;
            }
        }

        /* reconstruct missing atoms in base pairs */
        size_t bp;
        size_t r1,r2,c1,c2;
        for (bp=0;bp<bp_vec.size();bp++)
        {
            if (bp_vec[bp].second[4]) continue; // Watson-Crick pair only
            c1= bp_vec[bp].second[0];
            r1= bp_vec[bp].second[1];
            c2= bp_vec[bp].second[2];
            r2= bp_vec[bp].second[3];
            if (c1==c && fill_residue_map.count(r1) && 
                check3atomInPair(pdb_entry,c1,r1,c2,r2) &&
                fillMissingRNAatomInPair(c1,r1,c2,r2, pdb_entry,
                    fill_chain, ideal_rna, fill_residue_map[r1]))
                    fill_residue_map[r1].clear();
            if (c2==c && fill_residue_map.count(r2) &&
                check3atomInPair(pdb_entry,c2,r2,c1,r1) &&
                fillMissingRNAatomInPair(c2,r2,c1,r1, pdb_entry,
                    fill_chain, ideal_rna, fill_residue_map[r2]))
                    fill_residue_map[r2].clear();
        }
        
        /* reconstruct missing atoms in base pair stack */
        bool stack_prev,stack_next;
        for (bp=0;bp<bp_vec.size();bp++)
        {
            if (bp_vec[bp].second[4]) continue; // Watson-Crick pair only
            c1= bp_vec[bp].second[0];
            r1= bp_vec[bp].second[1];
            c2= bp_vec[bp].second[2];
            r2= bp_vec[bp].second[3];
            if (c1==c && fill_residue_map.count(r1) && 
                checkAtomInStack(bp_vec,c1,r1,c2,r2, stack_prev, stack_next) &&
                fillMissingRNAatomInStack(c1,r1,c2,r2, stack_prev, stack_next,
                    pdb_entry, fill_chain, ideal_rna, fill_residue_map[r1]))
                    fill_residue_map[r1].clear();
            if (c2==c && fill_residue_map.count(r2) && 
                checkAtomInStack(bp_vec,c2,r2,c1,r1, stack_prev, stack_next) &&
                fillMissingRNAatomInStack(c2,r2,c1,r1, stack_prev, stack_next,
                    pdb_entry, fill_chain, ideal_rna, fill_residue_map[r2]))
                    fill_residue_map[r2].clear();
        }
        

        /* reconstruct missing atoms */
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            if (fill_residue_map.count(r)==0) continue;
            cout<<"reconstructing residue "
                <<pdb_entry.chains[c].residues[r].resn
                <<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode
                <<" from "<<pdb_entry.chains[c].residues[r].atoms.size()<<" atoms"<<endl;
            fillMissingRNAatom(
                r, pdb_entry.chains[c], fill_chain, ideal_rna, fill_residue_map[r]);
        }

        /* append missing atoms to chain*/
        for (r=0;r<pdb_entry.chains[c].residues.size();r++){
            for (a=0;a<fill_chain.residues[r].atoms.size();a++){
                pdb_entry.chains[c].residues[r].atoms.push_back(fill_chain.residues[r].atoms[a]);
        	}
        }

        /* clean up */
        fill_atom.name.clear();
        fill_atom.xyz.clear();
        fill_residue.resn.clear();
        vector<AtomUnit>   ().swap(fill_residue.atoms);
        vector<ResidueUnit>().swap(fill_chain.residues);
    }

    /* clean up */
    missing_atom_vec.clear();
    nonstd_atom_vec.clear();
    remove_residue_vec.clear();
    map<size_t,vector<string> >().swap(fill_residue_map);
    return true;
}

bool MissingRNAatom(ModelUnit &pdb_entry, 
    map<string, map<string,vector<double> > >&ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec,
    const int option=0)
{
    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
        MissingRNAatom(pdb_entry, c, ideal_rna, bp_vec, option);
    return true;
}

#endif
