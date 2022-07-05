#!/usr/bin/python
''' create IdealRNA.hpp '''
from string import Template

atom_template=Template('        ${RESN}["$ATOM"]=tmp; ${RESN}["$ATOM"][0]=$X; ${RESN}["$ATOM"][1]=$Y; ${RESN}["$ATOM"][2]=$Z;\n')
atom2_template=Template('    ${RESN}["$ATOM"]=tmp; ${RESN}["$ATOM"][0]=$X; ${RESN}["$ATOM"][1]=$Y; ${RESN}["$ATOM"][2]=$Z;\n')

hpp_txt='''/* IdealRNA.hpp - Idealized RNA nucleotide conformations */
#include <string>
#include <map>
#include <vector>

using namespace std;

map<string, map<string,vector<float> > >parse_ideal_rna()
{
    vector<float> tmp(3,0);
    map<string, map<string,vector<float> > >ideal_rna;

'''

# generated using 
# fiber -rna -seq=ACGU ACGU.rna.pdb
fp=open("ACGU.rna.pdb",'r')
lines=fp.read().splitlines()
fp.close()


hpp_txt+="    map<string,vector<float> >P;\n"
for line in lines:
    if not line.startswith("ATOM  ") or line[21]!='A':
        continue
    atom=line[12:16]
    resi=line[22:26].strip()
    if (not (resi=="1" and atom in (" O3'")) and \
        not (resi=="2" and atom in (" P  "," OP1"," OP2"," O5'"))):
        continue
    x   =line[30:38]
    y   =line[38:46]
    z   =line[46:54]
    hpp_txt+=atom2_template.substitute(
        RESN="P",
        ATOM=atom,
        X=x,
        Y=y,
        Z=z,
    )
hpp_txt+='''    ideal_rna["P"]=P;
    map<string,vector<float> >().swap(P);

'''


chainID_dict={
    '  A':'A','  U':'B',
    '  C':'A','  G':'B',
}

for resn in ['  A','  C','  G','  U']:
    hpp_txt+="    map<string,vector<float> >%s;\n"%(resn.strip())
    for line in lines:
        if not line.startswith("ATOM  ") or line[77]=='H' or \
            line[17:20]!=resn or chainID_dict[resn]!=line[21]:
            continue
        atom=line[12:16]
        x   =line[30:38]
        y   =line[38:46]
        z   =line[46:54]
        hpp_txt+=atom2_template.substitute(
            RESN=resn.strip(),
            ATOM=atom,
            X=x,
            Y=y,
            Z=z,
        )
    hpp_txt+='''    ideal_rna["%s"]=%s;
    map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

# generated using
# fiber -rna -seq=AAACAGAUCACCCGCUGAGCGGGUUAUCUGUUAAAAACAAGAAUACAACCACGACUAGAAGCAGGAGUAUAAUCAUGAUUCAACACCAGCAUCCACCCCCGCCUCGACGCCGGCGUCUACUCCUGCUUGAAGACGAGGAUGCAGCCGCGGCUGGAGGCGGGGGUGUAGUCGUGGUUUAAUACUAGUAUUCAUCCUCGUCUUGAUGCUGGUGUUUAUUCUUGUUU 224.rna.pdb

fp=open("ACGU.rna.pdb",'r')
lines=fp.read().splitlines()
fp.close()

for c in range(2):
    chain="AB"[c]
    for r in range(3):
        key="%s%d"%(chain,r)
        hpp_txt+="    map<string,vector<float> >%s;\n"%key
        for line in lines:
            if not line.startswith("ATOM  ") or line[21]!=chain \
                or line[22:26]!="   %d"%(r+1+c*5):
                continue
            atom=line[12:16]
            if not atom in [" P  "," C4'"," C1'"]:
                continue
            x   =line[30:38]
            y   =line[38:46]
            z   =line[46:54]
            hpp_txt+=atom2_template.substitute(
                RESN=key,
                ATOM=atom,
                X=x,
                Y=y,
                Z=z,
            )
hpp_txt+='\n'
for resn1 in ['  A','  C','  G','  U']:
    for resn2 in ['  A','  C','  G','  U']:
        key=resn1+resn2
        hpp_txt+="    ideal_rna[\"%s0\"]=A0; ideal_rna[\"%s1\"]=A1; "%(
            key,key)
        hpp_txt+="ideal_rna[\"%s2\"]=B1; ideal_rna[\"%s3\"]=B2;\n"%(
            key,key)
for resn1 in ['  A','  C','  G','  U']:
    for resn2 in ['  A','  C','  G','  U']:
        for resn3 in ['  A','  C','  G','  U']:
            key=resn1+resn2+resn3
            hpp_txt+="    ideal_rna[\"%s0\"]=A0; ideal_rna[\"%s1\"]=A1; ideal_rna[\"%s2\"]=A2; "%(
                key,key,key)
            hpp_txt+=" ideal_rna[\"%s3\"]=B0; ideal_rna[\"%s4\"]=B1; ideal_rna[\"%s5\"]=B2;\n"%(
                key,key,key)
hpp_txt+='\n'
for r in range(3):
    hpp_txt+="    map<string,vector<float> >().swap(A%d);\n"%(r)
    hpp_txt+="    map<string,vector<float> >().swap(B%d);\n"%(r)

hpp_txt+='''
    vector<vector<float> > xyz_list1(3,tmp);
    vector<vector<float> > xyz_list2(3,tmp);
    vector<vector<float> > RotMatix;  // U
    vector<float> TranVect;  // t
    for (map<string, map<string,vector<float> > >::iterator iter = ideal_rna.begin();
        iter != ideal_rna.end(); iter++)
    {
        string key =  iter->first;
        if (key.size()<=3) continue;
        int idx=(char)(key[key.size()-1])-'0';
        bool reverse=false;
        if (key.size()==7 && idx>=2)
        {
            reverse=true;
            idx=(3-idx);
        }
        else if (key.size()==10 && idx>=3)
        {
            reverse=true;
            idx=(5-idx);
        }
        string nt=key.substr(idx*3,3);
        if (reverse)
        {
            if      (nt=="  A") nt="  U";
            else if (nt=="  C") nt="  G";
            else if (nt=="  G") nt="  C";
            else if (nt=="  U") nt="  A";
        }
    
        xyz_list1[0][0]=ideal_rna[nt][" P  "][0];
        xyz_list1[0][1]=ideal_rna[nt][" P  "][1];
        xyz_list1[0][2]=ideal_rna[nt][" P  "][2];
        xyz_list1[1][0]=ideal_rna[nt][" C4'"][0];
        xyz_list1[1][1]=ideal_rna[nt][" C4'"][1];
        xyz_list1[1][2]=ideal_rna[nt][" C4'"][2];
        xyz_list1[2][0]=ideal_rna[nt][" C1'"][0];
        xyz_list1[2][1]=ideal_rna[nt][" C1'"][1];
        xyz_list1[2][2]=ideal_rna[nt][" C1'"][2];
        
        xyz_list2[0][0]=ideal_rna[key][" P  "][0];
        xyz_list2[0][1]=ideal_rna[key][" P  "][1];
        xyz_list2[0][2]=ideal_rna[key][" P  "][2];
        xyz_list2[1][0]=ideal_rna[key][" C4'"][0];
        xyz_list2[1][1]=ideal_rna[key][" C4'"][1];
        xyz_list2[1][2]=ideal_rna[key][" C4'"][2];
        xyz_list2[2][0]=ideal_rna[key][" C1'"][0];
        xyz_list2[2][1]=ideal_rna[key][" C1'"][1];
        xyz_list2[2][2]=ideal_rna[key][" C1'"][2];
        
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
        
        for (map<string, vector<float> >::iterator i = ideal_rna[nt].begin();
            i != ideal_rna[nt].end(); i++)
        {
            string k=i-> first;
            if (ideal_rna[key].count(k)) continue;
            ideal_rna[key][k]=tmp;
            ChangeCoor(ideal_rna[nt][k], RotMatix, TranVect, ideal_rna[key][k]);
        }
    }

    vector<vector<float> >().swap(xyz_list1);
    vector<vector<float> >().swap(xyz_list2);
    vector<vector<float> >().swap(RotMatix);  // U
    vector<float> ().swap(TranVect);  // t
    vector<float> ().swap(tmp);
    return ideal_rna;
}
'''

fp=open("IdealRNA.hpp",'w')
fp.write(hpp_txt)
fp.close()
