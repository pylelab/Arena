/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef cssr_HPP
#define cssr_HPP 1
#include <cstring>
#include <set>
#include <map>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;

const double Pm_mu    = 53.755; const double Pm_sd    =42.773; // P[i-1]-P[i]-P[j]-P[j-1]
const double Pp_mu    = 80.233; const double Pp_sd    =20.466; // P[i+1]-P[i]-P[j]-P[j+1]
const double O5m_mu   = 60.047; const double O5m_sd   =38.731; // O5'[i-1]-O5'[i]-O5'[j]-O5'[j-1]
const double O5p_mu   = 96.645; const double O5p_sd   =20.185; // O5'[i+1]-O5'[i]-O5'[j]-O5'[j+1]
const double C5m_mu   = 58.526; const double C5m_sd   =36.222; // C5'[i-1]-C5'[i]-C5'[j]-C5'[j-1]
const double C5p_mu   =100.129; const double C5p_sd   =23.505; // C5'[i+1]-C5'[i]-C5'[j]-C5'[j+1]
const double C4m_mu   = 64.830; const double C4m_sd   =31.640; // C4'[i-1]-C4'[i]-C4'[j]-C4'[j-1]
const double C4p_mu   =115.546; const double C4p_sd   =23.989; // C4'[i+1]-C4'[i]-C4'[j]-C4'[j+1]
const double C3m_mu   = 80.193; const double C3m_sd   =26.156; // C3'[i-1]-C3'[i]-C3'[j]-C3'[j-1]
const double C3p_mu   =133.291; const double C3p_sd   =26.435; // C3'[i+1]-C3'[i]-C3'[j]-C3'[j+1]
const double C2m_mu   = 87.048; const double C2m_sd   =31.557; // C2'[i-1]-C2'[i]-C2'[j]-C2'[j-1]
const double C2p_mu   =147.550; const double C2p_sd   =29.374; // C2'[i+1]-C2'[i]-C2'[j]-C2'[j+1]
const double C1m_mu   = 74.283; const double C1m_sd   =38.789; // C1'[i-1]-C1'[i]-C1'[j]-C1'[j-1]
const double C1p_mu   =134.728; const double C1p_sd   =27.557; // C1'[i+1]-C1'[i]-C1'[j]-C1'[j+1]
const double O4m_mu   = 62.955; const double O4m_sd   =38.254; // O4'[i-1]-O4'[i]-O4'[j]-O4'[j-1]
const double O4p_mu   =113.863; const double O4p_sd   =24.303; // O4'[i+1]-O4'[i]-O4'[j]-O4'[j+1]
const double O3m_mu   = 80.432; const double O3m_sd   =22.356; // O3'[i-1]-O3'[i]-O3'[j]-O3'[j-1]
const double O3p_mu   =137.697; const double O3p_sd   =28.655; // O3'[i+1]-O3'[i]-O3'[j]-O3'[j+1]
const double P44P_mu  = -0.166; const double P44P_sd  =28.901; // P[i]-C4'[i]-C4'[j]-P[j]
const double C4114C_mu=172.313; const double C4114C_sd=43.508; // C4'[i]-C1'[i]-C1'[j]-C4'[j]
const double PP_mu    = 18.564; const double PP_sd    = 0.896; // P[i]-P[j]
const double O5O5_mu  = 16.808; const double O5O5_sd  = 0.774; // O5'[i]-O5'[j]
const double C5C5_mu  = 17.160; const double C5C5_sd  = 0.513; // C5'[i]-C5'[j]
const double C4C4_mu  = 15.018; const double C4C4_sd  = 0.419; // C4'[i]-C4'[j]
const double C3C3_mu  = 13.676; const double C3C3_sd  = 0.498; // C3'[i]-C3'[j]
const double C2C2_mu  = 11.046; const double C2C2_sd  = 0.515; // C2'[i]-C2'[j]
const double C1C1_mu  = 10.659; const double C1C1_sd  = 0.336; // C1'[i]-C1'[j]
const double O4O4_mu  = 13.232; const double O4O4_sd  = 0.399; // O4'[i]-O4'[j]
const double O3O3_mu  = 15.266; const double O3O3_sd  = 0.641; // O3'[i]-O3'[j]
const double NN_mu    =  8.929; const double NN_sd    = 0.250; // N[i]-N[j]
const double aPm_mu   =109.192; const double aPm_sd   =21.448; // <P[i-1]P[i],P[j+1]P[j]>
const double aPp_mu   =109.326; const double aPp_sd   =20.496; // <P[i+1]P[i],P[j-1]P[j]>
const double aO5m_mu  = 97.545; const double aO5m_sd  =18.629; // <O5'[i-1]O5'[i],O5'[j+1]O5'[j]>
const double aO5p_mu  = 97.816; const double aO5p_sd  =18.227; // <O5'[i+1]O5'[i],O5'[j-1]O5'[j]>
const double aC5m_mu  = 94.309; const double aC5m_sd  =19.724; // <C5'[i-1]C5'[i],C5'[j+1]C5'[j]>
const double aC5p_mu  = 94.640; const double aC5p_sd  =19.843; // <C5'[i+1]C5'[i],C5'[j-1]C5'[j]>
const double aC4m_mu  = 81.590; const double aC4m_sd  =17.368; // <C4'[i-1]C4'[i],C4'[j+1]C4'[j]>
const double aC4p_mu  = 81.352; const double aC4p_sd  =18.082; // <C4'[i+1]C4'[i],C4'[j-1]C4'[j]>
const double aC3m_mu  = 71.338; const double aC3m_sd  =17.754; // <C3'[i-1]C3'[i],C3'[j+1]C3'[j]>
const double aC3p_mu  = 70.590; const double aC3p_sd  =18.100; // <C3'[i+1]C3'[i],C3'[j-1]C3'[j]>
const double aC2m_mu  = 59.411; const double aC2m_sd  =20.249; // <C2'[i-1]C2'[i],C2'[j+1]C2'[j]>
const double aC2p_mu  = 58.352; const double aC2p_sd  =20.141; // <C2'[i+1]C2'[i],C2'[j-1]C2'[j]>
const double aC1m_mu  = 64.542; const double aC1m_sd  =19.817; // <C1'[i-1]C1'[i],C1'[j+1]C1'[j]>
const double aC1p_mu  = 63.616; const double aC1p_sd  =20.319; // <C1'[i+1]C1'[i],C1'[j-1]C1'[j]>
const double aO4m_mu  = 78.198; const double aO4m_sd  =17.830; // <O4'[i-1]O4'[i],O4'[j+1]O4'[j]>
const double aO4p_mu  = 77.672; const double aO4p_sd  =18.859; // <O4'[i+1]O4'[i],O4'[j-1]O4'[j]>
const double aO3m_mu  = 70.536; const double aO3m_sd  =17.393; // <O3'[i-1]O3'[i],O3'[j+1]O3'[j]>
const double aO3p_mu  = 69.919; const double aO3p_sd  =17.786; // <O3'[i+1]O3'[i],O3'[j-1]O3'[j]>
const double aPC_mu   = 59.764; const double aPC_sd   =18.397; // <P[i]C4'[i],P[j]C4'[j]>
const double aCC_mu   =165.213; const double aCC_sd   =13.607; // <C4'[i]C1'[i],C4'[j]C1'[j]>

/* if a residue is paired with multiple other residues, just keep the
 * highest scoring residue pair */
void filter_bp(vector<pair<double,vector<size_t> > >&bp_vec)
{
    vector<pair<size_t,size_t> > paired_res_vec;
    vector<pair<double,size_t> > bp_idx_vec;
    vector<size_t> filtered_bp_vec;
    size_t bp;
    for (bp=0;bp<bp_vec.size();bp++)
        bp_idx_vec.push_back(pair<double,size_t>(-bp_vec[bp].first,bp));
    sort(bp_idx_vec.begin(), bp_idx_vec.end());
    
    for (bp=0;bp<bp_idx_vec.size();bp++)
    {
        if (-bp_idx_vec[bp].first<0.5) break;
        if (find(paired_res_vec.begin(), paired_res_vec.end(), 
            pair<size_t,size_t>(bp_vec[bp_idx_vec[bp].second].second[0],
            bp_vec[bp_idx_vec[bp].second].second[1]))!=paired_res_vec.end())
            continue;
        if (find(paired_res_vec.begin(), paired_res_vec.end(),
            pair<size_t,size_t>(bp_vec[bp_idx_vec[bp].second].second[2],
            bp_vec[bp_idx_vec[bp].second].second[3]))!=paired_res_vec.end())
            continue;
        filtered_bp_vec.push_back(bp_idx_vec[bp].second);
        paired_res_vec.push_back(pair<size_t,size_t>(
                    bp_vec[bp_idx_vec[bp].second].second[0],
                    bp_vec[bp_idx_vec[bp].second].second[1]));
        paired_res_vec.push_back(pair<size_t,size_t>(
                    bp_vec[bp_idx_vec[bp].second].second[2],
                    bp_vec[bp_idx_vec[bp].second].second[3]));
    }
    vector<pair<size_t,size_t> >().swap(paired_res_vec);
    vector<pair<double,size_t> >().swap(bp_idx_vec);

    vector<pair<double,vector<size_t> > >bp_vec_tmp;
    for (bp=0;bp<bp_vec.size();bp++)
    {
        if (find(filtered_bp_vec.begin(), filtered_bp_vec.end(),bp
            )!=filtered_bp_vec.end()) bp_vec_tmp.push_back(bp_vec[bp]);
        bp_vec[bp].second.clear();
    }
    bp_vec.clear();
    for (bp=0;bp<bp_vec_tmp.size();bp++)
        bp_vec.push_back(bp_vec_tmp[bp]);

    vector<pair<double,vector<size_t> > >().swap(bp_vec_tmp);
    vector<size_t>().swap(filtered_bp_vec);
    return;
}

inline bool bp_tor_score(
    const vector<double> &c1, const vector<double> &c2,
    const vector<double> &c3, const vector<double> &c4,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    double tor=rad2deg(Points2Dihedral(c1,c2,c3,c4));
    if(tor<-180) return false;
    denominator+=weight;
    double diff=fabs(tor-mu); 
    if (diff>180) diff-=180;
    nominator+=weight*(1-diff/(tol* sd));
    return true;
}

inline bool bp_ang_score(
    const vector<double> &c1, const vector<double> &c2,
    const vector<double> &c3, const vector<double> &c4,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    double ang=rad2deg(Points4Angle(c1,c2,c3,c4));
    if(ang<-180) return false;
    denominator+=weight;
    nominator+=weight*(1-fabs(ang-mu)/(tol* sd));
    return true;
}

inline bool bp_len_score(const vector<double>&c1, const vector<double>&c2,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    denominator+=weight;
    nominator+=weight*(1-fabs(Points2Distance(c1,c2)-mu)/(tol*sd));
    return true;
}

inline bool bp_nn_score(const bool previnextj, const bool nextiprevj,
    const bool has_prev_ci,       const bool has_next_ci,
    const bool has_prev_cj,       const bool has_next_cj,
    const vector<double> &prev_ci, const vector<double> &next_ci,
    const vector<double> &prev_cj, const vector<double> &next_cj,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    if (!(has_prev_ci && has_next_cj) && !(has_next_ci && has_prev_cj))
        return false;
    denominator+=weight;
    double nominator_previnextj=0;
    double nominator_nextiprevj=0;
    if (has_prev_ci  && has_next_cj )
    {
        nominator_previnextj=1-fabs(Points2Distance(prev_ci,next_cj)-mu)/(tol*sd);
        if (!previnextj && nominator_previnextj>0) nominator_previnextj=0;
    }
    if (has_next_ci  && has_prev_cj )
    {
        nominator_nextiprevj=1-fabs(Points2Distance(next_ci,prev_cj)-mu)/(tol*sd);
        if (!nextiprevj && nominator_nextiprevj>0) nominator_nextiprevj=0;
    }

    if (nominator_previnextj>nominator_nextiprevj)
         nominator+=weight*nominator_previnextj;
    else nominator+=weight*nominator_nextiprevj;
    return true;
}

inline double mean_vec(const vector<double>&vec)
{
    double sum=0;
    int i;
    for (i=0;i<vec.size();i++) sum+=vec[i];
    return sum/vec.size();
}

void cssr(const ModelUnit &pdb_entry, vector<vector<size_t> >&res_str_vec,
    vector<pair<double,vector<size_t> > >&bp_vec, 
    const bool interchain=false)
{
    /* pre-trained parameters */
    double weight_len=1;
    double weight_nn =1;
    double weight_tor=1;
    double weight_ang=1;
    double tol=3; // tolerance, in the unit of standard deviation
    double adjust1=0.1;  // adjust for varying number of tests
    double adjust2=0;    // adjust for varying std
    double adjust3=-0.07;// baseline value
    double totaltest=(10-1)*(weight_len>0)+
                    (10-1)*(weight_nn >0)+
                    ( 9-1)*(weight_tor>0)+
                    (10-1)*(weight_ang>0);
    
    /* other variables */
    vector<size_t> tmp_bp(5,0);
    vector<size_t> tmp_nt(2,0);
    double ang;
    size_t c1, c2, r1, r2, a1, a2;
    char base1,base2;
    char base1prev,base1next;
    char base2prev,base2next;
    char icode,chainID;
    stringstream ss;
    double nominator, denominator;
    double nominator_previnextj;
    double nominator_nextiprevj;
    int   successtest;
    bool previnextj,nextiprevj;

    /* coordinates of previous residue */
    vector<double>prev_Pi(3,0);  bool has_prev_Pi;  vector<double>prev_Pj(3,0);  bool has_prev_Pj; 
    vector<double>prev_O5i(3,0); bool has_prev_O5i; vector<double>prev_O5j(3,0); bool has_prev_O5j;
    vector<double>prev_C5i(3,0); bool has_prev_C5i; vector<double>prev_C5j(3,0); bool has_prev_C5j;
    vector<double>prev_C4i(3,0); bool has_prev_C4i; vector<double>prev_C4j(3,0); bool has_prev_C4j;
    vector<double>prev_C3i(3,0); bool has_prev_C3i; vector<double>prev_C3j(3,0); bool has_prev_C3j;
    vector<double>prev_C2i(3,0); bool has_prev_C2i; vector<double>prev_C2j(3,0); bool has_prev_C2j;
    vector<double>prev_C1i(3,0); bool has_prev_C1i; vector<double>prev_C1j(3,0); bool has_prev_C1j;
    vector<double>prev_O4i(3,0); bool has_prev_O4i; vector<double>prev_O4j(3,0); bool has_prev_O4j;
    vector<double>prev_O3i(3,0); bool has_prev_O3i; vector<double>prev_O3j(3,0); bool has_prev_O3j;
    vector<double>prev_Nxi(3,0); bool has_prev_Nxi; vector<double>prev_Nxj(3,0); bool has_prev_Nxj;
    /* coordinates of current residue */
    vector<double>Pi(3,0);  bool has_Pi;  vector<double>Pj(3,0);  bool has_Pj; 
    vector<double>O5i(3,0); bool has_O5i; vector<double>O5j(3,0); bool has_O5j;
    vector<double>C5i(3,0); bool has_C5i; vector<double>C5j(3,0); bool has_C5j;
    vector<double>C4i(3,0); bool has_C4i; vector<double>C4j(3,0); bool has_C4j;
    vector<double>C3i(3,0); bool has_C3i; vector<double>C3j(3,0); bool has_C3j;
    vector<double>C2i(3,0); bool has_C2i; vector<double>C2j(3,0); bool has_C2j;
    vector<double>C1i(3,0); bool has_C1i; vector<double>C1j(3,0); bool has_C1j;
    vector<double>O4i(3,0); bool has_O4i; vector<double>O4j(3,0); bool has_O4j;
    vector<double>O3i(3,0); bool has_O3i; vector<double>O3j(3,0); bool has_O3j;
    vector<double>Nxi(3,0); bool has_Nxi; vector<double>Nxj(3,0); bool has_Nxj;
    /* coordinates of next residue */
    vector<double>next_Pi(3,0);  bool has_next_Pi;  vector<double>next_Pj(3,0);  bool has_next_Pj; 
    vector<double>next_O5i(3,0); bool has_next_O5i; vector<double>next_O5j(3,0); bool has_next_O5j;
    vector<double>next_C5i(3,0); bool has_next_C5i; vector<double>next_C5j(3,0); bool has_next_C5j;
    vector<double>next_C4i(3,0); bool has_next_C4i; vector<double>next_C4j(3,0); bool has_next_C4j;
    vector<double>next_C3i(3,0); bool has_next_C3i; vector<double>next_C3j(3,0); bool has_next_C3j;
    vector<double>next_C2i(3,0); bool has_next_C2i; vector<double>next_C2j(3,0); bool has_next_C2j;
    vector<double>next_C1i(3,0); bool has_next_C1i; vector<double>next_C1j(3,0); bool has_next_C1j;
    vector<double>next_O4i(3,0); bool has_next_O4i; vector<double>next_O4j(3,0); bool has_next_O4j;
    vector<double>next_O3i(3,0); bool has_next_O3i; vector<double>next_O3j(3,0); bool has_next_O3j;
    vector<double>next_Nxi(3,0); bool has_next_Nxi; vector<double>next_Nxj(3,0); bool has_next_Nxj;

    /* loop over base pairs */
    for (c1=0; c1<pdb_entry.chains.size(); c1++)
    {
        tmp_nt[0]=tmp_bp[0]=c1;
        for (r1=0; r1<pdb_entry.chains[c1].residues.size(); r1++)
        {
            base1=pdb_entry.chains[c1].residues[r1].resn[2];
            base1prev=base1next=0;
            if (!pdb_entry.chains[c1].residues[r1].atoms.size()) continue;
            //if ((chainID=pdb_entry.chains[c1].chainID)!=' ') 
                //ss<<chainID<<".";
            //ss<<base1<<pdb_entry.chains[c1].residues[r1].resi;
            //if ((icode=pdb_entry.chains[c1].residues[r1].icode)!=' ')
                //ss<<'^'<<icode;
            //tmp_bp[0]=ss.str();
            tmp_nt[1]=tmp_bp[1]=r1;
            res_str_vec.push_back(tmp_nt);
            ss.str(string());

            has_prev_Pi  = false;
            has_prev_O5i = false;
            has_prev_C5i = false;
            has_prev_C4i = false;
            has_prev_C3i = false;
            has_prev_C2i = false;
            has_prev_C1i = false;
            has_prev_O4i = false;
            has_prev_O3i = false;
            has_prev_Nxi = false;
            has_Pi       = false;
            has_O5i      = false;
            has_C5i      = false;
            has_C4i      = false;
            has_C3i      = false;
            has_C2i      = false;
            has_C1i      = false;
            has_O4i      = false;
            has_O3i      = false;
            has_Nxi      = false;
            has_next_Pi  = false;
            has_next_O5i = false;
            has_next_C5i = false;
            has_next_C4i = false;
            has_next_C3i = false;
            has_next_C2i = false;
            has_next_C1i = false;
            has_next_O4i = false;
            has_next_O3i = false;
            has_next_Nxi = false;

            for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
            {
                if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" P  ")
                {
                    has_Pi=true;
                    Pi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O5'")
                {
                    has_O5i=true;
                    O5i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C5'")
                {
                    has_C5i=true;
                    C5i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C4'")
                {
                    has_C4i=true;
                    C4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C3'")
                {
                    has_C3i=true;
                    C3i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C2'")
                {
                    has_C2i=true;
                    C2i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C1'")
                {
                    has_C1i=true;
                    C1i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O4'")
                {
                    has_O4i=true;
                    O4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O3'")
                {
                    has_O3i=true;
                    O3i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if ((pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N9 "
                        && (base1=='A' || base1=='G')) ||
                        (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N1 "
                        && (base1=='C' || base1=='T' || base1=='U')))
                {
                    has_Nxi=true;
                    Nxi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
            }

            if (r1>0 && (weight_tor>0 || weight_ang>0 || weight_nn>0))
            {
                base1prev=pdb_entry.chains[c1].residues[r1-1].resn[2];
                for (a1=0;a1<pdb_entry.chains[c1].residues[r1-1].atoms.size();a1++)
                {
                    if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" P  ")
                    {
                        has_prev_Pi=true;
                        prev_Pi=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O5'")
                    {
                        has_prev_O5i=true;
                        prev_O5i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C5'")
                    {
                        has_prev_C5i=true;
                        prev_C5i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C4'")
                    {
                        has_prev_C4i=true;
                        prev_C4i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C3'")
                    {
                        has_prev_C3i=true;
                        prev_C3i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C2'")
                    {
                        has_prev_C2i=true;
                        prev_C2i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C1'")
                    {
                        has_prev_C1i=true;
                        prev_C1i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O4'")
                    {
                        has_prev_O4i=true;
                        prev_O4i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O3'")
                    {
                        has_prev_O3i=true;
                        prev_O3i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if ((pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" N9 "
                            && (base1prev=='A' || base1prev=='G')) ||
                             (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" N1 "
                            && (base1prev=='C' || base1prev=='T' || base1prev=='U')))
                    {
                        has_prev_Nxi=true;
                        prev_Nxi=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                }
            }

            if (r1<pdb_entry.chains[c1].residues.size()-1 && (weight_tor>0 || weight_ang>0 || weight_nn>0))
            {
                base1next=pdb_entry.chains[c1].residues[r1+1].resn[2];
                for (a1=0;a1<pdb_entry.chains[c1].residues[r1+1].atoms.size();a1++)
                {
                    if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" P  ")
                    {
                        has_next_Pi=true;
                        next_Pi=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O5'")
                    {
                        has_next_O5i=true;
                        next_O5i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C5'")
                    {
                        has_next_C5i=true;
                        next_C5i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C4'")
                    {
                        has_next_C4i=true;
                        next_C4i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C3'")
                    {
                        has_next_C3i=true;
                        next_C3i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C2'")
                    {
                        has_next_C2i=true;
                        next_C2i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C1'")
                    {
                        has_next_C1i=true;
                        next_C1i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O4'")
                    {
                        has_next_O4i=true;
                        next_O4i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O3'")
                    {
                        has_next_O3i=true;
                        next_O3i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if ((pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" N9 "
                            && (base1next=='A' || base1next=='G')) ||
                             (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" N1 "
                            && (base1next=='C' || base1next=='T' || base1next=='U')))
                    {
                        has_next_Nxi=true;
                        next_Nxi=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                }
            }


            for (c2=c1+interchain; c2<pdb_entry.chains.size(); c2++)
            {
                tmp_bp[2]=c2;
                for (r2=(c1==c2)*(r1+1); r2<pdb_entry.chains[c2].residues.size(); r2++)
                {
                    if (!pdb_entry.chains[c2].residues[r2].atoms.size() ||
                        Points2Distance2(
                        pdb_entry.chains[c1].residues[r1].atoms[0].xyz,
                        pdb_entry.chains[c2].residues[r2].atoms[0].xyz)>530)
                        continue;
                    base2=pdb_entry.chains[c2].residues[r2].resn[2];
                    base2prev=base2next=0;
                    if ((base1=='A' &&(base2=='U' || base2=='T')) ||
                        (base1=='C' && base2=='G')||
                       ((base1=='U' || base1=='T')&& base2=='A')  ||
                        (base1=='G' && base2=='C'))
                        tmp_bp[4]=0;
                        //tmp_bp[2]=base1+string("-")+base2+" WC";
                    else if ((base1=='G' && base2=='U') ||
                             (base1=='U' && base2=='G')) 
                        tmp_bp[4]=1;
                        //tmp_bp[2]=base1+string("-")+base2+" Wobble";
                    else continue;

                    //if ((chainID=pdb_entry.chains[c2].chainID)!=' ') 
                        //ss<<chainID<<".";
                    //ss<<base2<<pdb_entry.chains[c2].residues[r2].resi;
                    //if ((icode=pdb_entry.chains[c2].residues[r2].icode)!=' ')
                        //ss<<'^'<<icode;
                    //tmp_bp[1]=ss.str();
                    tmp_bp[3]=r2;
                    //ss.str(string());

                    has_prev_Pj  = false;
                    has_prev_O5j = false;
                    has_prev_C5j = false;
                    has_prev_C4j = false;
                    has_prev_C3j = false;
                    has_prev_C2j = false;
                    has_prev_C1j = false;
                    has_prev_O4j = false;
                    has_prev_O3j = false;
                    has_Pj       = false;
                    has_O5j      = false;
                    has_C5j      = false;
                    has_C4j      = false;
                    has_C3j      = false;
                    has_C2j      = false;
                    has_C1j      = false;
                    has_O4j      = false;
                    has_O3j      = false;
                    has_Nxj      = false;
                    has_next_Pj  = false;
                    has_next_O5j = false;
                    has_next_C5j = false;
                    has_next_C4j = false;
                    has_next_C3j = false;
                    has_next_C2j = false;
                    has_next_C1j = false;
                    has_next_O4j = false;
                    has_next_O3j = false;

                    for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
                    {
                        if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" P  ")
                        {
                            has_Pj=true;
                            Pj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O5'")
                        {
                            has_O5j=true;
                            O5j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C5'")
                        {
                            has_C5j=true;
                            C5j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C4'")
                        {
                            has_C4j=true;
                            C4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C3'")
                        {
                            has_C3j=true;
                            C3j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C2'")
                        {
                            has_C2j=true;
                            C2j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C1'")
                        {
                            has_C1j=true;
                            C1j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O4'")
                        {
                            has_O4j=true;
                            O4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O3'")
                        {
                            has_O3j=true;
                            O3j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if ((pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N9 "
                                && (base2=='A' || base2=='G')) ||
                                (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N1 "
                                && (base2=='C' || base2=='T' || base2=='U')))
                        {
                            has_Nxj=true;
                            Nxj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                    }

                    if (r2>0 && (weight_tor>0 || weight_ang>0 || weight_nn>0) && c2-c1+r2-r1>1)
                    {
                        base2prev=pdb_entry.chains[c2].residues[r2-1].resn[2];
                        for (a2=0;a2<pdb_entry.chains[c2].residues[r2-1].atoms.size();a2++)
                        {
                            if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" P  ")
                            {
                                has_prev_Pj=true;
                                prev_Pj=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O5'")
                            {
                                has_prev_O5j=true;
                                prev_O5j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C5'")
                            {
                                has_prev_C5j=true;
                                prev_C5j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C4'")
                            {
                                has_prev_C4j=true;
                                prev_C4j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C3'")
                            {
                                has_prev_C3j=true;
                                prev_C3j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C2'")
                            {
                                has_prev_C2j=true;
                                prev_C2j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C1'")
                            {
                                has_prev_C1j=true;
                                prev_C1j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O4'")
                            {
                                has_prev_O4j=true;
                                prev_O4j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O3'")
                            {
                                has_prev_O3j=true;
                                prev_O3j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if ((pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" N9 "
                                  && (base2prev=='A' || base2prev=='G')) ||
                                     (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" N1 "
                                  && (base2prev=='C' || base2prev=='T' || base2prev=='U')))
                            {
                                has_prev_Nxj=true;
                                prev_Nxj=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                        }
                    }

                    if (r2<pdb_entry.chains[c2].residues.size()-1 && 
                       (weight_tor>0 || weight_ang>0 || weight_nn>0) && c2-c1+r2-r1>1)
                    {
                        base2next=pdb_entry.chains[c2].residues[r2+1].resn[2];
                        for (a2=0;a2<pdb_entry.chains[c2].residues[r2+1].atoms.size();a2++)
                        {
                            if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" P  ")
                            {
                                has_next_Pj=true;
                                next_Pj=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O5'")
                            {
                                has_next_O5j=true;
                                next_O5j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C5'")
                            {
                                has_next_C5j=true;
                                next_C5j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C4'")
                            {
                                has_next_C4j=true;
                                next_C4j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C3'")
                            {
                                has_next_C3j=true;
                                next_C3j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C2'")
                            {
                                has_next_C2j=true;
                                next_C2j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C1'")
                            {
                                has_next_C1j=true;
                                next_C1j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O4'")
                            {
                                has_next_O4j=true;
                                next_O4j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O3'")
                            {
                                has_next_O3j=true;
                                next_O3j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if ((pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" N9 "
                                  && (base2next=='A' || base2next=='G')) ||
                                     (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" N1 "
                                  && (base2next=='C' || base2next=='T' || base2next=='U')))
                            {
                                has_next_Nxj=true;
                                next_Nxj=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                        }
                    }

                    nominator=denominator=successtest=0;
                    if (weight_len>0)
                    {
                        successtest--;
                        if (has_Pi  && has_Pj ) successtest+=bp_len_score( Pi, Pj,  PP_mu,  PP_sd,tol,weight_len,nominator,denominator); //   PP: P[i]-P[j]
                        if (has_O5i && has_O5j) successtest+=bp_len_score(O5i,O5j,O5O5_mu,O5O5_sd,tol,weight_len,nominator,denominator); // O5O5: O5'[i]-O5'[j]
                        if (has_C5i && has_C5j) successtest+=bp_len_score(C5i,C5j,C5C5_mu,C5C5_sd,tol,weight_len,nominator,denominator); // C5C5: C5'[i]-C5'[j]
                        if (has_C4i && has_C4j) successtest+=bp_len_score(C4i,C4j,C4C4_mu,C4C4_sd,tol,weight_len,nominator,denominator); // C4C4: C4'[i]-C4'[j]
                        if (has_C3i && has_C3j) successtest+=bp_len_score(C3i,C3j,C3C3_mu,C3C3_sd,tol,weight_len,nominator,denominator); // C3C3: C3'[i]-C3'[j]
                        if (has_C2i && has_C2j) successtest+=bp_len_score(C2i,C2j,C2C2_mu,C2C2_sd,tol,weight_len,nominator,denominator); // C2C2: C2'[i]-C2'[j]
                        if (has_C1i && has_C1j) successtest+=bp_len_score(C1i,C1j,C1C1_mu,C1C1_sd,tol,weight_len,nominator,denominator); // C1C1: C1'[i]-C1'[j]
                        if (has_O4i && has_O4j) successtest+=bp_len_score(O4i,O4j,O4O4_mu,O4O4_sd,tol,weight_len,nominator,denominator); // O4O4: O4'[i]-O4'[j]
                        if (has_O3i && has_O3j) successtest+=bp_len_score(O3i,O3j,O3O3_mu,O3O3_sd,tol,weight_len,nominator,denominator); // O3O3: O3'[i]-O3'[j]
                        if (has_Nxi && has_Nxj) successtest+=bp_len_score(Nxi,Nxj,  NN_mu,  NN_sd,tol,weight_len,nominator,denominator); //   NN: N[i]-N[j]
                    }
                    if (weight_nn>0)
                    {
                        successtest--;
                        previnextj=((base1prev=='A' &&(base2next=='U' || base2next=='T')) ||
                                    (base1prev=='C' && base2next=='G')||
                                   ((base1prev=='U' || base1prev=='T')&& base2next=='A')  ||
                                    (base1prev=='G' && base2next=='C'));
                        nextiprevj=((base1next=='A' &&(base2prev=='U' || base2prev=='T')) ||
                                    (base1next=='C' && base2prev=='G')||
                                   ((base1next=='U' || base1next=='T')&& base2prev=='A')  ||
                                    (base1next=='G' && base2prev=='C'));
                        if (has_Pi  && has_Pj) //   PP: P[i]-P[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_Pi, has_next_Pi, has_prev_Pj, has_next_Pj,
                                prev_Pi, next_Pi, prev_Pj, next_Pj,   PP_mu,  PP_sd,tol,weight_nn,nominator,denominator);
                        if (has_O5i && has_O5j) // O5O5: O5'[i]-O5'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O5i,has_next_O5i,has_prev_O5j,has_next_O5j,
                                prev_O5i,next_O5i,prev_O5j,next_O5j,O5O5_mu,O5O5_sd,tol,weight_nn,nominator,denominator);
                        if (has_C5i && has_C5j) // C5C5: C5'[i]-C5'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C5i,has_next_C5i,has_prev_C5j,has_next_C5j,
                                prev_C5i,next_C5i,prev_C5j,next_C5j,C5C5_mu,C5C5_sd,tol,weight_nn,nominator,denominator);
                        if (has_C4i && has_C4j) // C4C4: C4'[i]-C4'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C4i,has_next_C4i,has_prev_C4j,has_next_C4j,
                                prev_C4i,next_C4i,prev_C4j,next_C4j,C4C4_mu,C4C4_sd,tol,weight_nn,nominator,denominator);
                        if (has_C3i && has_C3j) // C3C3: C3'[i]-C3'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C3i,has_next_C3i,has_prev_C3j,has_next_C3j,
                                prev_C3i,next_C3i,prev_C3j,next_C3j,C3C3_mu,C3C3_sd,tol,weight_nn,nominator,denominator);
                        if (has_C2i && has_C2j) // C2C2: C2'[i]-C2'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C2i,has_next_C2i,has_prev_C2j,has_next_C2j,
                                prev_C2i,next_C2i,prev_C2j,next_C2j,C2C2_mu,C2C2_sd,tol,weight_nn,nominator,denominator);
                        if (has_C1i && has_C1j) // C1C1: C1'[i]-C1'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C1i,has_next_C1i,has_prev_C1j,has_next_C1j,
                                prev_C1i,next_C1i,prev_C1j,next_C1j,C1C1_mu,C1C1_sd,tol,weight_nn,nominator,denominator);
                        if (has_O4i && has_O4j) // O4O4: O4'[i]-O4'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O4i,has_next_O4i,has_prev_O4j,has_next_O4j,
                                prev_O4i,next_O4i,prev_O4j,next_O4j,O4O4_mu,O4O4_sd,tol,weight_nn,nominator,denominator);
                        if (has_O3i && has_O3j) // O3O3: O3'[i]-O3'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O3i,has_next_O3i,has_prev_O3j,has_next_O3j,
                                prev_O3i,next_O3i,prev_O3j,next_O3j,O3O3_mu,O3O3_sd,tol,weight_nn,nominator,denominator);
                        if (has_Nxi && has_Nxj) //   NN: N[i]-N[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_Nxi,has_next_Nxi,has_prev_Nxj,has_next_Nxj,
                                prev_Nxi,next_Nxi,prev_Nxj,next_Nxj,  NN_mu,  NN_sd,tol,weight_nn,nominator,denominator);
                    }
                    if (weight_tor>0)
                    {
                        successtest--;
                        if (has_next_Pi  && has_Pi  && has_Pj  && has_next_Pj ) successtest+=bp_tor_score( next_Pi, Pi, Pj, next_Pj, Pp_mu, Pp_sd, tol,weight_tor,nominator,denominator); //  Pp:  P[i+1]-P[i]-P[j]-P[j-1]
                        if (has_next_O5i && has_O5i && has_O5j && has_next_O5j) successtest+=bp_tor_score(next_O5i,O5i,O5j,next_O5j,O5p_mu,O5p_sd, tol,weight_tor,nominator,denominator); // O5p: O5'[i+1]-O5'[i]-O5'[j]-O5'[j-1]
                        if (has_next_C5i && has_C5i && has_C5j && has_next_C5j) successtest+=bp_tor_score(next_C5i,C5i,C5j,next_C5j,C5p_mu,C5p_sd, tol,weight_tor,nominator,denominator); // C5p: C5'[i+1]-C5'[i]-C5'[j]-C5'[j-1]
                        if (has_next_C4i && has_C4i && has_C4j && has_next_C4j) successtest+=bp_tor_score(next_C4i,C4i,C4j,next_C4j,C4p_mu,C4p_sd, tol,weight_tor,nominator,denominator); // C4p: C4'[i+1]-C4'[i]-C4'[j]-C4'[j-1]
                        if (has_next_C3i && has_C3i && has_C3j && has_next_C3j) successtest+=bp_tor_score(next_C3i,C3i,C3j,next_C3j,C3p_mu,C3p_sd, tol,weight_tor,nominator,denominator); // C3p: C3'[i+1]-C3'[i]-C3'[j]-C3'[j-1]
                        if (has_next_C2i && has_C2i && has_C2j && has_next_C2j) successtest+=bp_tor_score(next_C2i,C2i,C2j,next_C2j,C2p_mu,C2p_sd, tol,weight_tor,nominator,denominator); // C2p: C2'[i+1]-C2'[i]-C2'[j]-C2'[j-1]
                        if (has_next_C1i && has_C1i && has_C1j && has_next_C1j) successtest+=bp_tor_score(next_C1i,C1i,C1j,next_C1j,C1p_mu,C1p_sd, tol,weight_tor,nominator,denominator); // C1p: C1'[i+1]-C1'[i]-C1'[j]-C1'[j-1]
                        if (has_next_O4i && has_O4i && has_O4j && has_next_O4j) successtest+=bp_tor_score(next_O4i,O4i,O4j,next_O4j,O4p_mu,O4p_sd, tol,weight_tor,nominator,denominator); // O4p: O4'[i+1]-O4'[i]-O4'[j]-O4'[j-1]
                        if (has_prev_O3i && has_O3i && has_O3j && has_prev_O3j) successtest+=bp_tor_score(prev_O3i,O3i,O3j,prev_O3j,O3m_mu,O3m_sd, tol,weight_tor,nominator,denominator); // O3m: O3'[i-1]-O3'[i]-O3'[j]-O3'[j+1]
                    }
                    if (weight_ang>0)
                    {
                        successtest--;
                        if (has_next_Pi  && has_Pi  && has_prev_Pj  && has_Pj ) successtest+=bp_ang_score( next_Pi, Pi, prev_Pj, Pj, aPp_mu, aPp_sd, tol,weight_ang,nominator,denominator); //  aPp:      <P[i+1]P[i],P[j-1]P[j]>
                        if (has_next_O5i && has_O5i && has_prev_O5j && has_O5j) successtest+=bp_ang_score(next_O5i,O5i,prev_O5j,O5j,aO5p_mu,aO5p_sd, tol,weight_ang,nominator,denominator); // aO5p: <O5'[i+1]O5'[i],O5'[j-1]O5'[j]> 
                        if (has_next_C5i && has_C5i && has_prev_C5j && has_C5j) successtest+=bp_ang_score(next_C5i,C5i,prev_C5j,C5j,aC5p_mu,aC5p_sd, tol,weight_ang,nominator,denominator); // aC5p: <C5'[i+1]C5'[i],C5'[j-1]C5'[j]> 
                        if (has_next_C4i && has_C4i && has_prev_C4j && has_C4j) successtest+=bp_ang_score(next_C4i,C4i,prev_C4j,C4j,aC4p_mu,aC4p_sd, tol,weight_ang,nominator,denominator); // aC4p: <C4'[i+1]C4'[i],C4'[j-1]C4'[j]> 
                        if (has_next_C3i && has_C3i && has_prev_C3j && has_C3j) successtest+=bp_ang_score(next_C3i,C3i,prev_C3j,C3j,aC3p_mu,aC3p_sd, tol,weight_ang,nominator,denominator); // aC3p: <C3'[i+1]C3'[i],C3'[j-1]C3'[j]> 
                        if (has_next_C2i && has_C2i && has_prev_C2j && has_C2j) successtest+=bp_ang_score(next_C2i,C2i,prev_C2j,C2j,aC2p_mu,aC2p_sd, tol,weight_ang,nominator,denominator); // aC2p: <C2'[i+1]C2'[i],C2'[j-1]C2'[j]> 
                        if (has_next_C1i && has_C1i && has_prev_C1j && has_C1j) successtest+=bp_ang_score(next_C1i,C1i,prev_C1j,C1j,aC1p_mu,aC1p_sd, tol,weight_ang,nominator,denominator); // aC1p: <C1'[i+1]C1'[i],C1'[j-1]C1'[j]> 
                        if (has_next_O4i && has_O4i && has_prev_O4j && has_O4j) successtest+=bp_ang_score(next_O4i,O4i,prev_O4j,O4j,aO4p_mu,aO4p_sd, tol,weight_ang,nominator,denominator); // aO4p: <O4'[i+1]O4'[i],O4'[j-1]O4'[j]> 
                        if (has_next_O3i && has_O3i && has_prev_O3j && has_O3j) successtest+=bp_ang_score(next_O3i,O3i,prev_O3j,O3j,aO3p_mu,aO3p_sd, tol,weight_ang,nominator,denominator); // aO3p: <O3'[i+1]O3'[i],O3'[j-1]O3'[j]> 
                        if (has_C4i      && has_C1i && has_C4j      && has_C1j)//successtest+=bp_ang_score(    C4i,C1i,     C4j,C1j, aCC_mu, aCC_sd, tol,weight_ang,nominator,denominator); //  aCC:    <C4'[i]C1'[i],C4'[j]C1'[j]>
                        {
                            ang=rad2deg(Points4Angle(C4i,C1i,C4j,C1j));
                            if(ang>=-180)
                            {
                                successtest++;
                                denominator+=weight_ang; 
                                if (ang>aCC_mu) nominator+=weight_ang;
                                else nominator+=weight_ang*(1-fabs(ang-aCC_mu)/(tol*aCC_sd));
                            }
                        }
                    }
                    
                    if (denominator>0) bp_vec.push_back(
                        pair<double,vector<size_t> >(nominator/denominator+
                        adjust1*successtest/totaltest+adjust3, tmp_bp));
                }
            }
        }
        //cerr<<endl;
    }

    /* clean up */
    vector<size_t>().swap(tmp_bp);

    prev_Pi.clear();  prev_Pj.clear();
    prev_O5i.clear(); prev_O5j.clear();
    prev_C5i.clear(); prev_C5j.clear();
    prev_C4i.clear(); prev_C4j.clear();
    prev_C3i.clear(); prev_C3j.clear();
    prev_C2i.clear(); prev_C2j.clear();
    prev_C1i.clear(); prev_C1j.clear();
    prev_O4i.clear(); prev_O4j.clear();
    prev_O3i.clear(); prev_O3j.clear();
    prev_Nxi.clear(); prev_Nxj.clear();

    Pi.clear();  Pj.clear();
    O5i.clear(); O5j.clear();
    C5i.clear(); C5j.clear();
    C4i.clear(); C4j.clear();
    C3i.clear(); C3j.clear();
    C2i.clear(); C2j.clear();
    C1i.clear(); C1j.clear();
    O4i.clear(); O4j.clear();
    O3i.clear(); O3j.clear();
    Nxi.clear(); Nxj.clear();

    next_Pi.clear();  next_Pj.clear();  
    next_O5i.clear(); next_O5j.clear();
    next_C5i.clear(); next_C5j.clear();
    next_C4i.clear(); next_C4j.clear();
    next_C3i.clear(); next_C3j.clear();
    next_C2i.clear(); next_C2j.clear();
    next_C1i.clear(); next_C1j.clear();
    next_O4i.clear(); next_O4j.clear();
    next_O3i.clear(); next_O3j.clear();
    next_Nxi.clear(); next_Nxj.clear();
    return;
}

#endif
