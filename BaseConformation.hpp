#ifndef BaseConformation_HPP
#define BaseConformation_HPP 1

#include <iostream>
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "cssr.hpp"

int fixBaseConformation(ResidueUnit &residue,
    map<string, map<string,vector<double> > >&ideal_rna)
{
    vector<double> tmp(3,0);
    size_t atomNum=residue.atoms.size()-11;
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    size_t a,b;

    int immovableCount=0;
    for (a=11;a<residue.atoms.size();a++)
        immovableCount+=(residue.atoms[a].movable==0);

    if (immovableCount==atomNum)
    {
        tmp.clear();
        return 0;
    }

    if (immovableCount>=3) atomNum=immovableCount;
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);

    b=0;
    for (a=11;a<residue.atoms.size();a++)
    {
        if (immovableCount>=3 && residue.atoms[a].movable) continue;
        xyz_list1[b][0]=ideal_rna[residue.resn][residue.atoms[a].name][0];
        xyz_list1[b][1]=ideal_rna[residue.resn][residue.atoms[a].name][1];
        xyz_list1[b][2]=ideal_rna[residue.resn][residue.atoms[a].name][2];

        xyz_list2[b][0]=residue.atoms[a].xyz[0];
        xyz_list2[b][1]=residue.atoms[a].xyz[1];
        xyz_list2[b][2]=residue.atoms[a].xyz[2];
        b++;
    }

    atomNum=countUniqAtom(xyz_list2);
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

    double rmsd=0;
    for (a=0;a<xyz_list1.size();a++)
    {
        ChangeCoor(xyz_list1[a], RotMatix, TranVect, tmp);
        rmsd+=Points2Distance(tmp,xyz_list2[a]);
    }
    rmsd=sqrt(rmsd);
    
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    
    int moved=0;
    if (rmsd>=0.1)
    {
        //cout<<residue.resn<<residue.resi<<residue.icode<<" rmsd="<<rmsd<<endl;
        for (a=11;a<residue.atoms.size();a++)
        {
            if (!residue.atoms[a].movable) continue;
            moved++;
            ChangeCoor(
                ideal_rna[residue.resn][residue.atoms[a].name],
                RotMatix, TranVect, residue.atoms[a].xyz);
            if (immovableCount>=3) residue.atoms[a].movable=false;
        }
    }

    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return moved;
}

int fixBasePairConformation(ResidueUnit &residue1, ResidueUnit &residue2,
    map<string, map<string,vector<double> > >&ideal_rna)
{
    vector<double> tmp(3,0);
    size_t atomNum=residue1.atoms.size()-11;
    atomNum      +=residue2.atoms.size()-11;
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    size_t a,b;

    int immovableCount=0;
    for (a=11;a<residue1.atoms.size();a++)
        immovableCount+=(residue1.atoms[a].movable==0);
    for (a=11;a<residue2.atoms.size();a++)
        immovableCount+=(residue2.atoms[a].movable==0);

    if (immovableCount==atomNum)
    {
        tmp.clear();
        return 0;
    }

    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);

    b=0;
    for (a=11;a<residue1.atoms.size();a++)
    {
        xyz_list1[b][0]=ideal_rna[residue1.resn][residue1.atoms[a].name][0];
        xyz_list1[b][1]=ideal_rna[residue1.resn][residue1.atoms[a].name][1];
        xyz_list1[b][2]=ideal_rna[residue1.resn][residue1.atoms[a].name][2];

        xyz_list2[b][0]=residue1.atoms[a].xyz[0];
        xyz_list2[b][1]=residue1.atoms[a].xyz[1];
        xyz_list2[b][2]=residue1.atoms[a].xyz[2];
        b++;
    }
    for (a=11;a<residue2.atoms.size();a++)
    {
        xyz_list1[b][0]=ideal_rna[residue2.resn][residue2.atoms[a].name][0];
        xyz_list1[b][1]=ideal_rna[residue2.resn][residue2.atoms[a].name][1];
        xyz_list1[b][2]=ideal_rna[residue2.resn][residue2.atoms[a].name][2];

        xyz_list2[b][0]=residue2.atoms[a].xyz[0];
        xyz_list2[b][1]=residue2.atoms[a].xyz[1];
        xyz_list2[b][2]=residue2.atoms[a].xyz[2];
        b++;
    }

    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<3)
    {
        for (a=0;a<xyz_list1.size();a++)
        {
            xyz_list1[a].clear();
            xyz_list2[a].clear();
        }
        xyz_list1.clear();
        xyz_list2.clear();
        tmp.clear();
        return 0;
    }

    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    
    double rmsd=0;
    for (a=0;a<xyz_list1.size();a++)
    {
        ChangeCoor(xyz_list1[a], RotMatix, TranVect, tmp);
        rmsd+=Points2Distance(tmp,xyz_list2[a]);
    }
    rmsd=sqrt(rmsd);
    
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    
    int moved=0;
    if (rmsd>=0.28)
    {
        //cout<<residue1.resn<<residue1.resi<<residue1.icode<<" :: "
        //    <<residue2.resn<<residue2.resi<<residue2.icode<<" rmsd="<<rmsd<<endl;
        for (a=11;a<residue1.atoms.size();a++)
        {
            if (!residue1.atoms[a].movable) continue;
            moved++;
            ChangeCoor(
                ideal_rna[residue1.resn][residue1.atoms[a].name],
                RotMatix, TranVect, residue1.atoms[a].xyz);
        }
        for (a=11;a<residue2.atoms.size();a++)
        {
            if (!residue2.atoms[a].movable) continue;
            moved++;
            ChangeCoor(
                ideal_rna[residue2.resn][residue2.atoms[a].name],
                RotMatix, TranVect, residue2.atoms[a].xyz);
        }
    }

    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return moved;
}

int fixBaseConformation(ModelUnit &pdb_entry, 
    map<string, map<string,vector<double> > >&ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec)
{
    size_t c,r;
    int moved=0;
    for (c=0;c<pdb_entry.chains.size();c++)
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
            moved+=fixBaseConformation(pdb_entry.chains[c].residues[r], ideal_rna);
    
    size_t bp,c1,c2,r1,r2;
    string resn_pair;
    for (bp=0;bp<bp_vec.size();bp++)
    {
        c1= bp_vec[bp].second[0];
        r1= bp_vec[bp].second[1];
        c2= bp_vec[bp].second[2];
        r2= bp_vec[bp].second[3];
        resn_pair=pdb_entry.chains[c1].residues[r1].resn+
                  pdb_entry.chains[c2].residues[r2].resn;
        if (resn_pair!="  A  U" && resn_pair!="  U  A" && 
            resn_pair!="  C  G" && resn_pair!="  G  C") continue;
        moved+=fixBasePairConformation(pdb_entry.chains[c1].residues[r1],
            pdb_entry.chains[c2].residues[r2], ideal_rna);
    }
    
    resn_pair.clear();
    return moved;
}

#endif
