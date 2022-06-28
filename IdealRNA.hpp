/* IdealRNA.hpp - Idealized RNA nucleotide conformations */
#include <string>
#include <map>
#include <vector>

#include "Superpose.hpp"
using namespace std;

/* option: 0 - ideal pdb; 1 - model pdb; 2 - 3dna fiber */
map<string, map<string,vector<float> > >parse_ideal_rna(int option=2)
{
    vector<float> tmp(3,0);
    map<string, map<string,vector<float> > >ideal_rna;
    
    map<string,vector<float> >A;
    A[" P  "]=tmp; A[" P  "][0]=   3.063; A[" P  "][1]=   8.025; A[" P  "][2]=  -4.135;
    A[" OP1"]=tmp; A[" OP1"][0]=   3.223; A[" OP1"][1]=   8.856; A[" OP1"][2]=  -5.350;
    A[" OP2"]=tmp; A[" OP2"][0]=   1.891; A[" OP2"][1]=   7.121; A[" OP2"][2]=  -4.118;
    A[" O5'"]=tmp; A[" O5'"][0]=   4.396; A[" O5'"][1]=   7.181; A[" O5'"][2]=  -3.881;
    A[" C5'"]=tmp; A[" C5'"][0]=   5.621; A[" C5'"][1]=   7.881; A[" C5'"][2]=  -3.587;
    A[" C4'"]=tmp; A[" C4'"][0]=   6.719; A[" C4'"][1]=   6.889; A[" C4'"][2]=  -3.258;
    A[" O4'"]=tmp; A[" O4'"][0]=   6.486; A[" O4'"][1]=   6.364; A[" O4'"][2]=  -1.919;
    A[" C3'"]=tmp; A[" C3'"][0]=   6.776; A[" C3'"][1]=   5.642; A[" C3'"][2]=  -4.140;
    A[" O3'"]=tmp; A[" O3'"][0]=   7.397; A[" O3'"][1]=   5.908; A[" O3'"][2]=  -5.390;
    A[" C2'"]=tmp; A[" C2'"][0]=   7.567; A[" C2'"][1]=   4.725; A[" C2'"][2]=  -3.208;
    A[" O2'"]=tmp; A[" O2'"][0]=   8.962; A[" O2'"][1]=   4.911; A[" O2'"][2]=  -3.068;
    A[" C1'"]=tmp; A[" C1'"][0]=   6.874; A[" C1'"][1]=   4.999; A[" C1'"][2]=  -1.877;
    A[" N9 "]=tmp; A[" N9 "][0]=   5.666; A[" N9 "][1]=   4.158; A[" N9 "][2]=  -1.647;
    A[" C8 "]=tmp; A[" C8 "][0]=   4.345; A[" C8 "][1]=   4.449; A[" C8 "][2]=  -1.900;
    A[" N7 "]=tmp; A[" N7 "][0]=   3.524; A[" N7 "][1]=   3.496; A[" N7 "][2]=  -1.584;
    A[" C5 "]=tmp; A[" C5 "][0]=   4.348; A[" C5 "][1]=   2.496; A[" C5 "][2]=  -1.086;
    A[" C6 "]=tmp; A[" C6 "][0]=   4.082; A[" C6 "][1]=   1.215; A[" C6 "][2]=  -0.578;
    A[" N6 "]=tmp; A[" N6 "][0]=   2.849; A[" N6 "][1]=   0.697; A[" N6 "][2]=  -0.485;
    A[" N1 "]=tmp; A[" N1 "][0]=   5.134; A[" N1 "][1]=   0.481; A[" N1 "][2]=  -0.167;
    A[" C2 "]=tmp; A[" C2 "][0]=   6.356; A[" C2 "][1]=   1.001; A[" C2 "][2]=  -0.262;
    A[" N3 "]=tmp; A[" N3 "][0]=   6.725; A[" N3 "][1]=   2.179; A[" N3 "][2]=  -0.716;
    A[" C4 "]=tmp; A[" C4 "][0]=   5.655; A[" C4 "][1]=   2.893; A[" C4 "][2]=  -1.121;
    ideal_rna["  A"]=A;
    map<string,vector<float> >().swap(A);

    map<string,vector<float> >C;
    C[" P  "]=tmp; C[" P  "][0]=   6.913; C[" P  "][1]=   5.099; C[" P  "][2]=  -6.683;
    C[" OP1"]=tmp; C[" OP1"][0]=   7.496; C[" OP1"][1]=   5.711; C[" OP1"][2]=  -7.898;
    C[" OP2"]=tmp; C[" OP2"][0]=   5.439; C[" OP2"][1]=   4.971; C[" OP2"][2]=  -6.666;
    C[" O5'"]=tmp; C[" O5'"][0]=   7.579; C[" O5'"][1]=   3.667; C[" O5'"][2]=  -6.429;
    C[" C5'"]=tmp; C[" C5'"][0]=   8.987; C[" C5'"][1]=   3.596; C[" C5'"][2]=  -6.135;
    C[" C4'"]=tmp; C[" C4'"][0]=   9.376; C[" C4'"][1]=   2.167; C[" C4'"][2]=  -5.806;
    C[" O4'"]=tmp; C[" O4'"][0]=   8.896; C[" O4'"][1]=   1.851; C[" O4'"][2]=  -4.467;
    C[" C3'"]=tmp; C[" C3'"][0]=   8.750; C[" C3'"][1]=   1.087; C[" C3'"][2]=  -6.688;
    C[" O3'"]=tmp; C[" O3'"][0]=   9.418; C[" O3'"][1]=   0.975; C[" O3'"][2]=  -7.938;
    C[" C2'"]=tmp; C[" C2'"][0]=   8.920; C[" C2'"][1]=  -0.112; C[" C2'"][2]=  -5.756;
    C[" O2'"]=tmp; C[" O2'"][0]=  10.194; C[" O2'"][1]=  -0.710; C[" O2'"][2]=  -5.617;
    C[" C1'"]=tmp; C[" C1'"][0]=   8.486; C[" C1'"][1]=   0.494; C[" C1'"][2]=  -4.425;
    C[" N1 "]=tmp; C[" N1 "][0]=   7.014; C[" N1 "][1]=   0.438; C[" N1 "][2]=  -4.195;
    C[" C2 "]=tmp; C[" C2 "][0]=   6.487; C[" C2 "][1]=  -0.729; C[" C2 "][2]=  -3.648;
    C[" O2 "]=tmp; C[" O2 "][0]=   7.253; C[" O2 "][1]=  -1.661; C[" O2 "][2]=  -3.379;
    C[" N3 "]=tmp; C[" N3 "][0]=   5.148; C[" N3 "][1]=  -0.800; C[" N3 "][2]=  -3.431;
    C[" C4 "]=tmp; C[" C4 "][0]=   4.352; C[" C4 "][1]=   0.233; C[" C4 "][2]=  -3.736;
    C[" N4 "]=tmp; C[" N4 "][0]=   3.054; C[" N4 "][1]=   0.114; C[" N4 "][2]=  -3.504;
    C[" C5 "]=tmp; C[" C5 "][0]=   4.874; C[" C5 "][1]=   1.442; C[" C5 "][2]=  -4.299;
    C[" C6 "]=tmp; C[" C6 "][0]=   6.214; C[" C6 "][1]=   1.492; C[" C6 "][2]=  -4.509;
    ideal_rna["  C"]=C;
    map<string,vector<float> >().swap(C);

    map<string,vector<float> >G;
    G[" P  "]=tmp; G[" P  "][0]=   8.572; G[" P  "][1]=   0.556; G[" P  "][2]=  -9.231;
    G[" OP1"]=tmp; G[" OP1"][0]=   9.394; G[" OP1"][1]=   0.756; G[" OP1"][2]= -10.446;
    G[" OP2"]=tmp; G[" OP2"][0]=   7.262; G[" OP2"][1]=   1.245; G[" OP2"][2]=  -9.214;
    G[" O5'"]=tmp; G[" O5'"][0]=   8.359; G[" O5'"][1]=  -1.009; G[" O5'"][2]=  -8.977;
    G[" C5'"]=tmp; G[" C5'"][0]=   9.505; G[" C5'"][1]=  -1.830; G[" C5'"][2]=  -8.683;
    G[" C4'"]=tmp; G[" C4'"][0]=   9.061; G[" C4'"][1]=  -3.241; G[" C4'"][2]=  -8.354;
    G[" O4'"]=tmp; G[" O4'"][0]=   8.485; G[" O4'"][1]=  -3.248; G[" O4'"][2]=  -7.015;
    G[" C3'"]=tmp; G[" C3'"][0]=   7.951; G[" C3'"][1]=  -3.812; G[" C3'"][2]=  -9.236;
    G[" O3'"]=tmp; G[" O3'"][0]=   8.452; G[" O3'"][1]=  -4.267; G[" O3'"][2]= -10.486;
    G[" C2'"]=tmp; G[" C2'"][0]=   7.447; G[" C2'"][1]=  -4.913; G[" C2'"][2]=  -8.304;
    G[" O2'"]=tmp; G[" O2'"][0]=   8.195; G[" O2'"][1]=  -6.103; G[" O2'"][2]=  -8.164;
    G[" C1'"]=tmp; G[" C1'"][0]=   7.407; G[" C1'"][1]=  -4.169; G[" C1'"][2]=  -6.973;
    G[" N9 "]=tmp; G[" N9 "][0]=   6.139; G[" N9 "][1]=  -3.421; G[" N9 "][2]=  -6.743;
    G[" C8 "]=tmp; G[" C8 "][0]=   5.851; G[" C8 "][1]=  -2.097; G[" C8 "][2]=  -6.995;
    G[" N7 "]=tmp; G[" N7 "][0]=   4.628; G[" N7 "][1]=  -1.748; G[" N7 "][2]=  -6.675;
    G[" C5 "]=tmp; G[" C5 "][0]=   4.069; G[" C5 "][1]=  -2.923; G[" C5 "][2]=  -6.175;
    G[" C6 "]=tmp; G[" C6 "][0]=   2.765; G[" C6 "][1]=  -3.169; G[" C6 "][2]=  -5.670;
    G[" O6 "]=tmp; G[" O6 "][0]=   1.822; G[" O6 "][1]=  -2.391; G[" O6 "][2]=  -5.558;
    G[" N1 "]=tmp; G[" N1 "][0]=   2.618; G[" N1 "][1]=  -4.505; G[" N1 "][2]=  -5.268;
    G[" C2 "]=tmp; G[" C2 "][0]=   3.600; G[" C2 "][1]=  -5.472; G[" C2 "][2]=  -5.344;
    G[" N2 "]=tmp; G[" N2 "][0]=   3.260; G[" N2 "][1]=  -6.688; G[" N2 "][2]=  -4.908;
    G[" N3 "]=tmp; G[" N3 "][0]=   4.821; G[" N3 "][1]=  -5.239; G[" N3 "][2]=  -5.818;
    G[" C4 "]=tmp; G[" C4 "][0]=   4.982; G[" C4 "][1]=  -3.949; G[" C4 "][2]=  -6.213;
    ideal_rna["  G"]=G;
    map<string,vector<float> >().swap(G);

    map<string,vector<float> >U;
    U[" P  "]=tmp; U[" P  "][0]=   7.514; U[" P  "][1]=  -4.163; U[" P  "][2]= -11.779;
    U[" OP1"]=tmp; U[" OP1"][0]=   8.313; U[" OP1"][1]=  -4.439; U[" OP1"][2]= -12.994;
    U[" OP2"]=tmp; U[" OP2"][0]=   6.784; U[" OP2"][1]=  -2.876; U[" OP2"][2]= -11.762;
    U[" O5'"]=tmp; U[" O5'"][0]=   6.490; U[" O5'"][1]=  -5.365; U[" O5'"][2]= -11.525;
    U[" C5'"]=tmp; U[" C5'"][0]=   7.010; U[" C5'"][1]=  -6.675; U[" C5'"][2]= -11.231;
    U[" C4'"]=tmp; U[" C4'"][0]=   5.874; U[" C4'"][1]=  -7.623; U[" C4'"][2]= -10.902;
    U[" O4'"]=tmp; U[" O4'"][0]=   5.387; U[" O4'"][1]=  -7.318; U[" O4'"][2]=  -9.563;
    U[" C3'"]=tmp; U[" C3'"][0]=   4.631; U[" C3'"][1]=  -7.503; U[" C3'"][2]= -11.784;
    U[" O3'"]=tmp; U[" O3'"][0]=   4.807; U[" O3'"][1]=  -8.156; U[" O3'"][2]= -13.034;
    U[" C2'"]=tmp; U[" C2'"][0]=   3.612; U[" C2'"][1]=  -8.157; U[" C2'"][2]= -10.852;
    U[" O2'"]=tmp; U[" O2'"][0]=   3.598; U[" O2'"][1]=  -9.565; U[" O2'"][2]= -10.712;
    U[" C1'"]=tmp; U[" C1'"][0]=   3.981; U[" C1'"][1]=  -7.510; U[" C1'"][2]=  -9.521;
    U[" N1 "]=tmp; U[" N1 "][0]=   3.318; U[" N1 "][1]=  -6.195; U[" N1 "][2]=  -9.291;
    U[" C2 "]=tmp; U[" C2 "][0]=   2.055; U[" C2 "][1]=  -6.218; U[" C2 "][2]=  -8.751;
    U[" O2 "]=tmp; U[" O2 "][0]=   1.475; U[" O2 "][1]=  -7.251; U[" O2 "][2]=  -8.463;
    U[" N3 "]=tmp; U[" N3 "][0]=   1.473; U[" N3 "][1]=  -4.982; U[" N3 "][2]=  -8.552;
    U[" C4 "]=tmp; U[" C4 "][0]=   2.034; U[" C4 "][1]=  -3.755; U[" C4 "][2]=  -8.841;
    U[" O4 "]=tmp; U[" O4 "][0]=   1.415; U[" O4 "][1]=  -2.713; U[" O4 "][2]=  -8.618;
    U[" C5 "]=tmp; U[" C5 "][0]=   3.362; U[" C5 "][1]=  -3.834; U[" C5 "][2]=  -9.404;
    U[" C6 "]=tmp; U[" C6 "][0]=   3.953; U[" C6 "][1]=  -5.023; U[" C6 "][2]=  -9.608;
    ideal_rna["  U"]=U;
    map<string,vector<float> >().swap(U);

    map<string,vector<float> >AA0;
    AA0[" P  "]=tmp; AA0[" P  "][0]=   3.063; AA0[" P  "][1]=   8.025; AA0[" P  "][2]=  -4.135;
    AA0[" C4'"]=tmp; AA0[" C4'"][0]=   6.719; AA0[" C4'"][1]=   6.889; AA0[" C4'"][2]=  -3.258;
    AA0[" C1'"]=tmp; AA0[" C1'"][0]=   6.874; AA0[" C1'"][1]=   4.999; AA0[" C1'"][2]=  -1.877;
    ideal_rna["  A  A0"]=AA0;
    map<string,vector<float> >().swap(AA0);

    map<string,vector<float> >AA1;
    AA1[" P  "]=tmp; AA1[" P  "][0]=   6.913; AA1[" P  "][1]=   5.099; AA1[" P  "][2]=  -6.683;
    AA1[" C4'"]=tmp; AA1[" C4'"][0]=   9.376; AA1[" C4'"][1]=   2.167; AA1[" C4'"][2]=  -5.806;
    AA1[" C1'"]=tmp; AA1[" C1'"][0]=   8.486; AA1[" C1'"][1]=   0.493; AA1[" C1'"][2]=  -4.425;
    ideal_rna["  A  A1"]=AA1;
    map<string,vector<float> >().swap(AA1);

    map<string,vector<float> >AC0;
    AC0[" P  "]=tmp; AC0[" P  "][0]=   8.572; AC0[" P  "][1]=   0.556; AC0[" P  "][2]=  -9.231;
    AC0[" C4'"]=tmp; AC0[" C4'"][0]=   9.061; AC0[" C4'"][1]=  -3.241; AC0[" C4'"][2]=  -8.354;
    AC0[" C1'"]=tmp; AC0[" C1'"][0]=   7.407; AC0[" C1'"][1]=  -4.169; AC0[" C1'"][2]=  -6.973;
    ideal_rna["  A  C0"]=AC0;
    map<string,vector<float> >().swap(AC0);

    map<string,vector<float> >AC1;
    AC1[" P  "]=tmp; AC1[" P  "][0]=   7.514; AC1[" P  "][1]=  -4.163; AC1[" P  "][2]= -11.779;
    AC1[" C4'"]=tmp; AC1[" C4'"][0]=   5.874; AC1[" C4'"][1]=  -7.623; AC1[" C4'"][2]= -10.902;
    AC1[" C1'"]=tmp; AC1[" C1'"][0]=   3.981; AC1[" C1'"][1]=  -7.510; AC1[" C1'"][2]=  -9.521;
    ideal_rna["  A  C1"]=AC1;
    map<string,vector<float> >().swap(AC1);

    map<string,vector<float> >AG0;
    AG0[" P  "]=tmp; AG0[" P  "][0]=   4.074; AG0[" P  "][1]=  -7.563; AG0[" P  "][2]= -14.327;
    AG0[" C4'"]=tmp; AG0[" C4'"][0]=   0.825; AG0[" C4'"][1]=  -9.588; AG0[" C4'"][2]= -13.450;
    AG0[" C1'"]=tmp; AG0[" C1'"][0]=  -0.707; AG0[" C1'"][1]=  -8.471; AG0[" C1'"][2]= -12.069;
    ideal_rna["  A  G0"]=AG0;
    map<string,vector<float> >().swap(AG0);

    map<string,vector<float> >AG1;
    AG1[" P  "]=tmp; AG1[" P  "][0]=  -0.658; AG1[" P  "][1]=  -8.565; AG1[" P  "][2]= -16.875;
    AG1[" C4'"]=tmp; AG1[" C4'"][0]=  -4.486; AG1[" C4'"][1]=  -8.514; AG1[" C4'"][2]= -15.998;
    AG1[" C1'"]=tmp; AG1[" C1'"][0]=  -5.171; AG1[" C1'"][1]=  -6.746; AG1[" C1'"][2]= -14.617;
    ideal_rna["  A  G1"]=AG1;
    map<string,vector<float> >().swap(AG1);

    map<string,vector<float> >AU0;
    AU0[" P  "]=tmp; AU0[" P  "][0]=  -5.180; AU0[" P  "][1]=  -6.852; AU0[" P  "][2]= -19.423;
    AU0[" C4'"]=tmp; AU0[" C4'"][0]=  -8.374; AU0[" C4'"][1]=  -4.741; AU0[" C4'"][2]= -18.546;
    AU0[" C1'"]=tmp; AU0[" C1'"][0]=  -7.996; AU0[" C1'"][1]=  -2.883; AU0[" C1'"][2]= -17.165;
    ideal_rna["  A  U0"]=AU0;
    map<string,vector<float> >().swap(AU0);

    map<string,vector<float> >AU1;
    AU1[" P  "]=tmp; AU1[" P  "][0]=  -8.061; AU1[" P  "][1]=  -2.967; AU1[" P  "][2]= -21.971;
    AU1[" C4'"]=tmp; AU1[" C4'"][0]=  -9.608; AU1[" C4'"][1]=   0.534; AU1[" C4'"][2]= -21.094;
    AU1[" C1'"]=tmp; AU1[" C1'"][0]=  -8.286; AU1[" C1'"][1]=   1.894; AU1[" C1'"][2]= -19.713;
    ideal_rna["  A  U1"]=AU1;
    map<string,vector<float> >().swap(AU1);

    map<string,vector<float> >CA0;
    CA0[" P  "]=tmp; CA0[" P  "][0]=  -8.387; CA0[" P  "][1]=   1.858; CA0[" P  "][2]= -24.519;
    CA0[" C4'"]=tmp; CA0[" C4'"][0]=  -7.797; CA0[" C4'"][1]=   5.640; CA0[" C4'"][2]= -23.642;
    CA0[" C1'"]=tmp; CA0[" C1'"][0]=  -5.950; CA0[" C1'"][1]=   6.070; CA0[" C1'"][2]= -22.261;
    ideal_rna["  C  A0"]=CA0;
    map<string,vector<float> >().swap(CA0);

    map<string,vector<float> >CA1;
    CA1[" P  "]=tmp; CA1[" P  "][0]=  -6.054; CA1[" P  "][1]=   6.094; CA1[" P  "][2]= -27.067;
    CA1[" C4'"]=tmp; CA1[" C4'"][0]=  -3.514; CA1[" C4'"][1]=   8.959; CA1[" C4'"][2]= -26.190;
    CA1[" C1'"]=tmp; CA1[" C1'"][0]=  -1.727; CA1[" C1'"][1]=   8.323; CA1[" C1'"][2]= -24.809;
    ideal_rna["  C  A1"]=CA1;
    map<string,vector<float> >().swap(CA1);

    map<string,vector<float> >CC0;
    CC0[" P  "]=tmp; CC0[" P  "][0]=  -1.802; CC0[" P  "][1]=   8.399; CC0[" P  "][2]= -29.615;
    CC0[" C4'"]=tmp; CC0[" C4'"][0]=   1.883; CC0[" C4'"][1]=   9.437; CC0[" C4'"][2]= -28.738;
    CC0[" C1'"]=tmp; CC0[" C1'"][0]=   3.042; CC0[" C1'"][1]=   7.937; CC0[" C1'"][2]= -27.357;
    ideal_rna["  C  C0"]=CC0;
    map<string,vector<float> >().swap(CC0);

    map<string,vector<float> >CC1;
    CC1[" P  "]=tmp; CC1[" P  "][0]=   3.021; CC1[" P  "][1]=   8.041; CC1[" P  "][2]= -32.163;
    CC1[" C4'"]=tmp; CC1[" C4'"][0]=   6.683; CC1[" C4'"][1]=   6.924; CC1[" C4'"][2]= -31.286;
    CC1[" C1'"]=tmp; CC1[" C1'"][0]=   6.848; CC1[" C1'"][1]=   5.036; CC1[" C1'"][2]= -29.905;
    ideal_rna["  C  C1"]=CC1;
    map<string,vector<float> >().swap(CC1);

    map<string,vector<float> >CG0;
    CG0[" P  "]=tmp; CG0[" P  "][0]=   6.886; CG0[" P  "][1]=   5.135; CG0[" P  "][2]= -34.711;
    CG0[" C4'"]=tmp; CG0[" C4'"][0]=   9.364; CG0[" C4'"][1]=   2.216; CG0[" C4'"][2]= -33.834;
    CG0[" C1'"]=tmp; CG0[" C1'"][0]=   8.483; CG0[" C1'"][1]=   0.538; CG0[" C1'"][2]= -32.453;
    ideal_rna["  C  G0"]=CG0;
    map<string,vector<float> >().swap(CG0);

    map<string,vector<float> >CG1;
    CG1[" P  "]=tmp; CG1[" P  "][0]=   8.569; CG1[" P  "][1]=   0.601; CG1[" P  "][2]= -37.259;
    CG1[" C4'"]=tmp; CG1[" C4'"][0]=   9.078; CG1[" C4'"][1]=  -3.194; CG1[" C4'"][2]= -36.382;
    CG1[" C1'"]=tmp; CG1[" C1'"][0]=   7.429; CG1[" C1'"][1]=  -4.130; CG1[" C1'"][2]= -35.001;
    ideal_rna["  C  G1"]=CG1;
    map<string,vector<float> >().swap(CG1);

    map<string,vector<float> >CU0;
    CU0[" P  "]=tmp; CU0[" P  "][0]=   7.535; CU0[" P  "][1]=  -4.124; CU0[" P  "][2]= -39.807;
    CU0[" C4'"]=tmp; CU0[" C4'"][0]=   5.913; CU0[" C4'"][1]=  -7.592; CU0[" C4'"][2]= -38.930;
    CU0[" C1'"]=tmp; CU0[" C1'"][0]=   4.021; CU0[" C1'"][1]=  -7.489; CU0[" C1'"][2]= -37.549;
    ideal_rna["  C  U0"]=CU0;
    map<string,vector<float> >().swap(CU0);

    map<string,vector<float> >CU1;
    CU1[" P  "]=tmp; CU1[" P  "][0]=   4.113; CU1[" P  "][1]=  -7.541; CU1[" P  "][2]= -42.355;
    CU1[" C4'"]=tmp; CU1[" C4'"][0]=   0.875; CU1[" C4'"][1]=  -9.583; CU1[" C4'"][2]= -41.478;
    CU1[" C1'"]=tmp; CU1[" C1'"][0]=  -0.663; CU1[" C1'"][1]=  -8.474; CU1[" C1'"][2]= -40.097;
    ideal_rna["  C  U1"]=CU1;
    map<string,vector<float> >().swap(CU1);

    map<string,vector<float> >GA0;
    GA0[" P  "]=tmp; GA0[" P  "][0]=  -0.613; GA0[" P  "][1]=  -8.568; GA0[" P  "][2]= -44.903;
    GA0[" C4'"]=tmp; GA0[" C4'"][0]=  -4.441; GA0[" C4'"][1]=  -8.537; GA0[" C4'"][2]= -44.026;
    GA0[" C1'"]=tmp; GA0[" C1'"][0]=  -5.136; GA0[" C1'"][1]=  -6.773; GA0[" C1'"][2]= -42.645;
    ideal_rna["  G  A0"]=GA0;
    map<string,vector<float> >().swap(GA0);

    map<string,vector<float> >GA1;
    GA1[" P  "]=tmp; GA1[" P  "][0]=  -5.145; GA1[" P  "][1]=  -6.879; GA1[" P  "][2]= -47.451;
    GA1[" C4'"]=tmp; GA1[" C4'"][0]=  -8.349; GA1[" C4'"][1]=  -4.785; GA1[" C4'"][2]= -46.574;
    GA1[" C1'"]=tmp; GA1[" C1'"][0]=  -7.981; GA1[" C1'"][1]=  -2.925; GA1[" C1'"][2]= -45.193;
    ideal_rna["  G  A1"]=GA1;
    map<string,vector<float> >().swap(GA1);

    map<string,vector<float> >GC0;
    GC0[" P  "]=tmp; GC0[" P  "][0]=  -8.046; GC0[" P  "][1]=  -3.010; GC0[" P  "][2]= -49.999;
    GC0[" C4'"]=tmp; GC0[" C4'"][0]=  -9.611; GC0[" C4'"][1]=   0.484; GC0[" C4'"][2]= -49.122;
    GC0[" C1'"]=tmp; GC0[" C1'"][0]=  -8.296; GC0[" C1'"][1]=   1.850; GC0[" C1'"][2]= -47.741;
    ideal_rna["  G  C0"]=GC0;
    map<string,vector<float> >().swap(GC0);

    map<string,vector<float> >GC1;
    GC1[" P  "]=tmp; GC1[" P  "][0]=  -8.396; GC1[" P  "][1]=   1.814; GC1[" P  "][2]= -52.547;
    GC1[" C4'"]=tmp; GC1[" C4'"][0]=  -7.826; GC1[" C4'"][1]=   5.600; GC1[" C4'"][2]= -51.670;
    GC1[" C1'"]=tmp; GC1[" C1'"][0]=  -5.982; GC1[" C1'"][1]=   6.039; GC1[" C1'"][2]= -50.289;
    ideal_rna["  G  C1"]=GC1;
    map<string,vector<float> >().swap(GC1);

    map<string,vector<float> >GG0;
    GG0[" P  "]=tmp; GG0[" P  "][0]=  -6.086; GG0[" P  "][1]=   6.062; GG0[" P  "][2]= -55.095;
    GG0[" C4'"]=tmp; GG0[" C4'"][0]=  -3.561; GG0[" C4'"][1]=   8.940; GG0[" C4'"][2]= -54.218;
    GG0[" C1'"]=tmp; GG0[" C1'"][0]=  -1.771; GG0[" C1'"][1]=   8.313; GG0[" C1'"][2]= -52.837;
    ideal_rna["  G  G0"]=GG0;
    map<string,vector<float> >().swap(GG0);

    map<string,vector<float> >GG1;
    GG1[" P  "]=tmp; GG1[" P  "][0]=  -1.846; GG1[" P  "][1]=   8.389; GG1[" P  "][2]= -57.643;
    GG1[" C4'"]=tmp; GG1[" C4'"][0]=   1.834; GG1[" C4'"][1]=   9.447; GG1[" C4'"][2]= -56.766;
    GG1[" C1'"]=tmp; GG1[" C1'"][0]=   3.001; GG1[" C1'"][1]=   7.953; GG1[" C1'"][2]= -55.385;
    ideal_rna["  G  G1"]=GG1;
    map<string,vector<float> >().swap(GG1);

    map<string,vector<float> >GU0;
    GU0[" P  "]=tmp; GU0[" P  "][0]=   2.979; GU0[" P  "][1]=   8.057; GU0[" P  "][2]= -60.191;
    GU0[" C4'"]=tmp; GU0[" C4'"][0]=   6.646; GU0[" C4'"][1]=   6.959; GU0[" C4'"][2]= -59.314;
    GU0[" C1'"]=tmp; GU0[" C1'"][0]=   6.822; GU0[" C1'"][1]=   5.071; GU0[" C1'"][2]= -57.933;
    ideal_rna["  G  U0"]=GU0;
    map<string,vector<float> >().swap(GU0);

    map<string,vector<float> >GU1;
    GU1[" P  "]=tmp; GU1[" P  "][0]=   6.859; GU1[" P  "][1]=   5.171; GU1[" P  "][2]= -62.739;
    GU1[" C4'"]=tmp; GU1[" C4'"][0]=   9.353; GU1[" C4'"][1]=   2.265; GU1[" C4'"][2]= -61.862;
    GU1[" C1'"]=tmp; GU1[" C1'"][0]=   8.480; GU1[" C1'"][1]=   0.582; GU1[" C1'"][2]= -60.481;
    ideal_rna["  G  U1"]=GU1;
    map<string,vector<float> >().swap(GU1);

    map<string,vector<float> >UA0;
    UA0[" P  "]=tmp; UA0[" P  "][0]=   8.566; UA0[" P  "][1]=   0.645; UA0[" P  "][2]= -65.287;
    UA0[" C4'"]=tmp; UA0[" C4'"][0]=   9.094; UA0[" C4'"][1]=  -3.146; UA0[" C4'"][2]= -64.410;
    UA0[" C1'"]=tmp; UA0[" C1'"][0]=   7.450; UA0[" C1'"][1]=  -4.092; UA0[" C1'"][2]= -63.029;
    ideal_rna["  U  A0"]=UA0;
    map<string,vector<float> >().swap(UA0);

    map<string,vector<float> >UA1;
    UA1[" P  "]=tmp; UA1[" P  "][0]=   7.557; UA1[" P  "][1]=  -4.084; UA1[" P  "][2]= -67.835;
    UA1[" C4'"]=tmp; UA1[" C4'"][0]=   5.953; UA1[" C4'"][1]=  -7.561; UA1[" C4'"][2]= -66.958;
    UA1[" C1'"]=tmp; UA1[" C1'"][0]=   4.059; UA1[" C1'"][1]=  -7.468; UA1[" C1'"][2]= -65.577;
    ideal_rna["  U  A1"]=UA1;
    map<string,vector<float> >().swap(UA1);

    map<string,vector<float> >UC0;
    UC0[" P  "]=tmp; UC0[" P  "][0]=   4.153; UC0[" P  "][1]=  -7.520; UC0[" P  "][2]= -70.383;
    UC0[" C4'"]=tmp; UC0[" C4'"][0]=   0.925; UC0[" C4'"][1]=  -9.578; UC0[" C4'"][2]= -69.506;
    UC0[" C1'"]=tmp; UC0[" C1'"][0]=  -0.619; UC0[" C1'"][1]=  -8.477; UC0[" C1'"][2]= -68.125;
    ideal_rna["  U  C0"]=UC0;
    map<string,vector<float> >().swap(UC0);

    map<string,vector<float> >UC1;
    UC1[" P  "]=tmp; UC1[" P  "][0]=  -0.568; UC1[" P  "][1]=  -8.571; UC1[" P  "][2]= -72.931;
    UC1[" C4'"]=tmp; UC1[" C4'"][0]=  -4.396; UC1[" C4'"][1]=  -8.560; UC1[" C4'"][2]= -72.054;
    UC1[" C1'"]=tmp; UC1[" C1'"][0]=  -5.100; UC1[" C1'"][1]=  -6.800; UC1[" C1'"][2]= -70.673;
    ideal_rna["  U  C1"]=UC1;
    map<string,vector<float> >().swap(UC1);

    map<string,vector<float> >UG0;
    UG0[" P  "]=tmp; UG0[" P  "][0]=  -5.108; UG0[" P  "][1]=  -6.906; UG0[" P  "][2]= -75.479;
    UG0[" C4'"]=tmp; UG0[" C4'"][0]=  -8.324; UG0[" C4'"][1]=  -4.828; UG0[" C4'"][2]= -74.602;
    UG0[" C1'"]=tmp; UG0[" C1'"][0]=  -7.966; UG0[" C1'"][1]=  -2.966; UG0[" C1'"][2]= -73.221;
    ideal_rna["  U  G0"]=UG0;
    map<string,vector<float> >().swap(UG0);

    map<string,vector<float> >UG1;
    UG1[" P  "]=tmp; UG1[" P  "][0]=  -8.030; UG1[" P  "][1]=  -3.052; UG1[" P  "][2]= -78.027;
    UG1[" C4'"]=tmp; UG1[" C4'"][0]=  -9.613; UG1[" C4'"][1]=   0.434; UG1[" C4'"][2]= -77.150;
    UG1[" C1'"]=tmp; UG1[" C1'"][0]=  -8.306; UG1[" C1'"][1]=   1.807; UG1[" C1'"][2]= -75.769;
    ideal_rna["  U  G1"]=UG1;
    map<string,vector<float> >().swap(UG1);

    map<string,vector<float> >UU0;
    UU0[" P  "]=tmp; UU0[" P  "][0]=  -8.406; UU0[" P  "][1]=   1.770; UU0[" P  "][2]= -80.575;
    UU0[" C4'"]=tmp; UU0[" C4'"][0]=  -7.855; UU0[" C4'"][1]=   5.559; UU0[" C4'"][2]= -79.698;
    UU0[" C1'"]=tmp; UU0[" C1'"][0]=  -6.013; UU0[" C1'"][1]=   6.008; UU0[" C1'"][2]= -78.317;
    ideal_rna["  U  U0"]=UU0;
    map<string,vector<float> >().swap(UU0);

    map<string,vector<float> >UU1;
    UU1[" P  "]=tmp; UU1[" P  "][0]=  -6.117; UU1[" P  "][1]=   6.031; UU1[" P  "][2]= -83.123;
    UU1[" C4'"]=tmp; UU1[" C4'"][0]=  -3.607; UU1[" C4'"][1]=   8.921; UU1[" C4'"][2]= -82.246;
    UU1[" C1'"]=tmp; UU1[" C1'"][0]=  -1.815; UU1[" C1'"][1]=   8.304; UU1[" C1'"][2]= -80.865;
    ideal_rna["  U  U1"]=UU1;
    map<string,vector<float> >().swap(UU1);

    map<string,vector<float> >AAA0;
    AAA0[" P  "]=tmp; AAA0[" P  "][0]=  -1.890; AAA0[" P  "][1]=   8.380; AAA0[" P  "][2]= -85.671;
    AAA0[" C4'"]=tmp; AAA0[" C4'"][0]=   1.784; AAA0[" C4'"][1]=   9.456; AAA0[" C4'"][2]= -84.794;
    AAA0[" C1'"]=tmp; AAA0[" C1'"][0]=   2.959; AAA0[" C1'"][1]=   7.968; AAA0[" C1'"][2]= -83.413;
    ideal_rna["  A  A  A0"]=AAA0;
    map<string,vector<float> >().swap(AAA0);

    map<string,vector<float> >AAA1;
    AAA1[" P  "]=tmp; AAA1[" P  "][0]=   2.937; AAA1[" P  "][1]=   8.072; AAA1[" P  "][2]= -88.219;
    AAA1[" C4'"]=tmp; AAA1[" C4'"][0]=   6.610; AAA1[" C4'"][1]=   6.994; AAA1[" C4'"][2]= -87.342;
    AAA1[" C1'"]=tmp; AAA1[" C1'"][0]=   6.795; AAA1[" C1'"][1]=   5.107; AAA1[" C1'"][2]= -85.961;
    ideal_rna["  A  A  A1"]=AAA1;
    map<string,vector<float> >().swap(AAA1);

    map<string,vector<float> >AAA2;
    AAA2[" P  "]=tmp; AAA2[" P  "][0]=   6.832; AAA2[" P  "][1]=   5.207; AAA2[" P  "][2]= -90.767;
    AAA2[" C4'"]=tmp; AAA2[" C4'"][0]=   9.341; AAA2[" C4'"][1]=   2.314; AAA2[" C4'"][2]= -89.890;
    AAA2[" C1'"]=tmp; AAA2[" C1'"][0]=   8.477; AAA2[" C1'"][1]=   0.626; AAA2[" C1'"][2]= -88.509;
    ideal_rna["  A  A  A2"]=AAA2;
    map<string,vector<float> >().swap(AAA2);

    map<string,vector<float> >AAC0;
    AAC0[" P  "]=tmp; AAC0[" P  "][0]=   8.562; AAC0[" P  "][1]=   0.690; AAC0[" P  "][2]= -93.315;
    AAC0[" C4'"]=tmp; AAC0[" C4'"][0]=   9.110; AAC0[" C4'"][1]=  -3.099; AAC0[" C4'"][2]= -92.438;
    AAC0[" C1'"]=tmp; AAC0[" C1'"][0]=   7.472; AAC0[" C1'"][1]=  -4.052; AAC0[" C1'"][2]= -91.057;
    ideal_rna["  A  A  C0"]=AAC0;
    map<string,vector<float> >().swap(AAC0);

    map<string,vector<float> >AAC1;
    AAC1[" P  "]=tmp; AAC1[" P  "][0]=   7.578; AAC1[" P  "][1]=  -4.045; AAC1[" P  "][2]= -95.863;
    AAC1[" C4'"]=tmp; AAC1[" C4'"][0]=   5.993; AAC1[" C4'"][1]=  -7.529; AAC1[" C4'"][2]= -94.986;
    AAC1[" C1'"]=tmp; AAC1[" C1'"][0]=   4.098; AAC1[" C1'"][1]=  -7.447; AAC1[" C1'"][2]= -93.605;
    ideal_rna["  A  A  C1"]=AAC1;
    map<string,vector<float> >().swap(AAC1);

    map<string,vector<float> >AAC2;
    AAC2[" P  "]=tmp; AAC2[" P  "][0]=   4.192; AAC2[" P  "][1]=  -7.498; AAC2[" P  "][2]= -98.411;
    AAC2[" C4'"]=tmp; AAC2[" C4'"][0]=   0.975; AAC2[" C4'"][1]=  -9.573; AAC2[" C4'"][2]= -97.534;
    AAC2[" C1'"]=tmp; AAC2[" C1'"][0]=  -0.574; AAC2[" C1'"][1]=  -8.481; AAC2[" C1'"][2]= -96.153;
    ideal_rna["  A  A  C2"]=AAC2;
    map<string,vector<float> >().swap(AAC2);

    map<string,vector<float> >AAG0;
    AAG0[" P  "]=tmp; AAG0[" P  "][0]=  -0.523; AAG0[" P  "][1]=  -8.574; AAG0[" P  "][2]=-100.959;
    AAG0[" C4'"]=tmp; AAG0[" C4'"][0]=  -4.351; AAG0[" C4'"][1]=  -8.583; AAG0[" C4'"][2]=-100.082;
    AAG0[" C1'"]=tmp; AAG0[" C1'"][0]=  -5.065; AAG0[" C1'"][1]=  -6.826; AAG0[" C1'"][2]= -98.701;
    ideal_rna["  A  A  G0"]=AAG0;
    map<string,vector<float> >().swap(AAG0);

    map<string,vector<float> >AAG1;
    AAG1[" P  "]=tmp; AAG1[" P  "][0]=  -5.072; AAG1[" P  "][1]=  -6.933; AAG1[" P  "][2]=-103.507;
    AAG1[" C4'"]=tmp; AAG1[" C4'"][0]=  -8.299; AAG1[" C4'"][1]=  -4.872; AAG1[" C4'"][2]=-102.630;
    AAG1[" C1'"]=tmp; AAG1[" C1'"][0]=  -7.950; AAG1[" C1'"][1]=  -3.008; AAG1[" C1'"][2]=-101.249;
    ideal_rna["  A  A  G1"]=AAG1;
    map<string,vector<float> >().swap(AAG1);

    map<string,vector<float> >AAG2;
    AAG2[" P  "]=tmp; AAG2[" P  "][0]=  -8.014; AAG2[" P  "][1]=  -3.094; AAG2[" P  "][2]=-106.055;
    AAG2[" C4'"]=tmp; AAG2[" C4'"][0]=  -9.615; AAG2[" C4'"][1]=   0.384; AAG2[" C4'"][2]=-105.178;
    AAG2[" C1'"]=tmp; AAG2[" C1'"][0]=  -8.315; AAG2[" C1'"][1]=   1.763; AAG2[" C1'"][2]=-103.797;
    ideal_rna["  A  A  G2"]=AAG2;
    map<string,vector<float> >().swap(AAG2);

    map<string,vector<float> >AAU0;
    AAU0[" P  "]=tmp; AAU0[" P  "][0]=  -8.415; AAU0[" P  "][1]=   1.726; AAU0[" P  "][2]=-108.603;
    AAU0[" C4'"]=tmp; AAU0[" C4'"][0]=  -7.884; AAU0[" C4'"][1]=   5.517; AAU0[" C4'"][2]=-107.726;
    AAU0[" C1'"]=tmp; AAU0[" C1'"][0]=  -6.045; AAU0[" C1'"][1]=   5.976; AAU0[" C1'"][2]=-106.345;
    ideal_rna["  A  A  U0"]=AAU0;
    map<string,vector<float> >().swap(AAU0);

    map<string,vector<float> >AAU1;
    AAU1[" P  "]=tmp; AAU1[" P  "][0]=  -6.149; AAU1[" P  "][1]=   5.998; AAU1[" P  "][2]=-111.151;
    AAU1[" C4'"]=tmp; AAU1[" C4'"][0]=  -3.654; AAU1[" C4'"][1]=   8.902; AAU1[" C4'"][2]=-110.274;
    AAU1[" C1'"]=tmp; AAU1[" C1'"][0]=  -1.858; AAU1[" C1'"][1]=   8.294; AAU1[" C1'"][2]=-108.893;
    ideal_rna["  A  A  U1"]=AAU1;
    map<string,vector<float> >().swap(AAU1);

    map<string,vector<float> >AAU2;
    AAU2[" P  "]=tmp; AAU2[" P  "][0]=  -1.934; AAU2[" P  "][1]=   8.370; AAU2[" P  "][2]=-113.699;
    AAU2[" C4'"]=tmp; AAU2[" C4'"][0]=   1.734; AAU2[" C4'"][1]=   9.465; AAU2[" C4'"][2]=-112.822;
    AAU2[" C1'"]=tmp; AAU2[" C1'"][0]=   2.917; AAU2[" C1'"][1]=   7.984; AAU2[" C1'"][2]=-111.441;
    ideal_rna["  A  A  U2"]=AAU2;
    map<string,vector<float> >().swap(AAU2);

    map<string,vector<float> >ACA0;
    ACA0[" P  "]=tmp; ACA0[" P  "][0]=   2.894; ACA0[" P  "][1]=   8.088; ACA0[" P  "][2]=-116.247;
    ACA0[" C4'"]=tmp; ACA0[" C4'"][0]=   6.573; ACA0[" C4'"][1]=   7.028; ACA0[" C4'"][2]=-115.370;
    ACA0[" C1'"]=tmp; ACA0[" C1'"][0]=   6.768; ACA0[" C1'"][1]=   5.142; ACA0[" C1'"][2]=-113.989;
    ideal_rna["  A  C  A0"]=ACA0;
    map<string,vector<float> >().swap(ACA0);

    map<string,vector<float> >ACA1;
    ACA1[" P  "]=tmp; ACA1[" P  "][0]=   6.805; ACA1[" P  "][1]=   5.242; ACA1[" P  "][2]=-118.795;
    ACA1[" C4'"]=tmp; ACA1[" C4'"][0]=   9.328; ACA1[" C4'"][1]=   2.363; ACA1[" C4'"][2]=-117.918;
    ACA1[" C1'"]=tmp; ACA1[" C1'"][0]=   8.473; ACA1[" C1'"][1]=   0.671; ACA1[" C1'"][2]=-116.537;
    ideal_rna["  A  C  A1"]=ACA1;
    map<string,vector<float> >().swap(ACA1);

    map<string,vector<float> >ACA2;
    ACA2[" P  "]=tmp; ACA2[" P  "][0]=   8.558; ACA2[" P  "][1]=   0.735; ACA2[" P  "][2]=-121.343;
    ACA2[" C4'"]=tmp; ACA2[" C4'"][0]=   9.127; ACA2[" C4'"][1]=  -3.051; ACA2[" C4'"][2]=-120.466;
    ACA2[" C1'"]=tmp; ACA2[" C1'"][0]=   7.493; ACA2[" C1'"][1]=  -4.013; ACA2[" C1'"][2]=-119.085;
    ideal_rna["  A  C  A2"]=ACA2;
    map<string,vector<float> >().swap(ACA2);

    map<string,vector<float> >ACC0;
    ACC0[" P  "]=tmp; ACC0[" P  "][0]=   7.599; ACC0[" P  "][1]=  -4.005; ACC0[" P  "][2]=-123.891;
    ACC0[" C4'"]=tmp; ACC0[" C4'"][0]=   6.032; ACC0[" C4'"][1]=  -7.498; ACC0[" C4'"][2]=-123.014;
    ACC0[" C1'"]=tmp; ACC0[" C1'"][0]=   4.137; ACC0[" C1'"][1]=  -7.425; ACC0[" C1'"][2]=-121.633;
    ideal_rna["  A  C  C0"]=ACC0;
    map<string,vector<float> >().swap(ACC0);

    map<string,vector<float> >ACC1;
    ACC1[" P  "]=tmp; ACC1[" P  "][0]=   4.231; ACC1[" P  "][1]=  -7.476; ACC1[" P  "][2]=-126.439;
    ACC1[" C4'"]=tmp; ACC1[" C4'"][0]=   1.025; ACC1[" C4'"][1]=  -9.568; ACC1[" C4'"][2]=-125.562;
    ACC1[" C1'"]=tmp; ACC1[" C1'"][0]=  -0.529; ACC1[" C1'"][1]=  -8.484; ACC1[" C1'"][2]=-124.181;
    ideal_rna["  A  C  C1"]=ACC1;
    map<string,vector<float> >().swap(ACC1);

    map<string,vector<float> >ACC2;
    ACC2[" P  "]=tmp; ACC2[" P  "][0]=  -0.478; ACC2[" P  "][1]=  -8.577; ACC2[" P  "][2]=-128.987;
    ACC2[" C4'"]=tmp; ACC2[" C4'"][0]=  -4.306; ACC2[" C4'"][1]=  -8.606; ACC2[" C4'"][2]=-128.110;
    ACC2[" C1'"]=tmp; ACC2[" C1'"][0]=  -5.028; ACC2[" C1'"][1]=  -6.853; ACC2[" C1'"][2]=-126.729;
    ideal_rna["  A  C  C2"]=ACC2;
    map<string,vector<float> >().swap(ACC2);

    map<string,vector<float> >ACG0;
    ACG0[" P  "]=tmp; ACG0[" P  "][0]=  -5.036; ACG0[" P  "][1]=  -6.959; ACG0[" P  "][2]=-131.535;
    ACG0[" C4'"]=tmp; ACG0[" C4'"][0]=  -8.273; ACG0[" C4'"][1]=  -4.915; ACG0[" C4'"][2]=-130.658;
    ACG0[" C1'"]=tmp; ACG0[" C1'"][0]=  -7.934; ACG0[" C1'"][1]=  -3.050; ACG0[" C1'"][2]=-129.277;
    ideal_rna["  A  C  G0"]=ACG0;
    map<string,vector<float> >().swap(ACG0);

    map<string,vector<float> >ACG1;
    ACG1[" P  "]=tmp; ACG1[" P  "][0]=  -7.997; ACG1[" P  "][1]=  -3.136; ACG1[" P  "][2]=-134.083;
    ACG1[" C4'"]=tmp; ACG1[" C4'"][0]=  -9.617; ACG1[" C4'"][1]=   0.333; ACG1[" C4'"][2]=-133.206;
    ACG1[" C1'"]=tmp; ACG1[" C1'"][0]=  -8.324; ACG1[" C1'"][1]=   1.719; ACG1[" C1'"][2]=-131.825;
    ideal_rna["  A  C  G1"]=ACG1;
    map<string,vector<float> >().swap(ACG1);

    map<string,vector<float> >ACG2;
    ACG2[" P  "]=tmp; ACG2[" P  "][0]=  -8.424; ACG2[" P  "][1]=   1.682; ACG2[" P  "][2]=-136.631;
    ACG2[" C4'"]=tmp; ACG2[" C4'"][0]=  -7.913; ACG2[" C4'"][1]=   5.476; ACG2[" C4'"][2]=-135.754;
    ACG2[" C1'"]=tmp; ACG2[" C1'"][0]=  -6.076; ACG2[" C1'"][1]=   5.944; ACG2[" C1'"][2]=-134.373;
    ideal_rna["  A  C  G2"]=ACG2;
    map<string,vector<float> >().swap(ACG2);

    map<string,vector<float> >ACU0;
    ACU0[" P  "]=tmp; ACU0[" P  "][0]=  -6.180; ACU0[" P  "][1]=   5.966; ACU0[" P  "][2]=-139.179;
    ACU0[" C4'"]=tmp; ACU0[" C4'"][0]=  -3.701; ACU0[" C4'"][1]=   8.883; ACU0[" C4'"][2]=-138.302;
    ACU0[" C1'"]=tmp; ACU0[" C1'"][0]=  -1.901; ACU0[" C1'"][1]=   8.285; ACU0[" C1'"][2]=-136.921;
    ideal_rna["  A  C  U0"]=ACU0;
    map<string,vector<float> >().swap(ACU0);

    map<string,vector<float> >ACU1;
    ACU1[" P  "]=tmp; ACU1[" P  "][0]=  -1.977; ACU1[" P  "][1]=   8.359; ACU1[" P  "][2]=-141.727;
    ACU1[" C4'"]=tmp; ACU1[" C4'"][0]=   1.685; ACU1[" C4'"][1]=   9.474; ACU1[" C4'"][2]=-140.850;
    ACU1[" C1'"]=tmp; ACU1[" C1'"][0]=   2.875; ACU1[" C1'"][1]=   7.999; ACU1[" C1'"][2]=-139.469;
    ideal_rna["  A  C  U1"]=ACU1;
    map<string,vector<float> >().swap(ACU1);

    map<string,vector<float> >ACU2;
    ACU2[" P  "]=tmp; ACU2[" P  "][0]=   2.852; ACU2[" P  "][1]=   8.103; ACU2[" P  "][2]=-144.275;
    ACU2[" C4'"]=tmp; ACU2[" C4'"][0]=   6.536; ACU2[" C4'"][1]=   7.063; ACU2[" C4'"][2]=-143.398;
    ACU2[" C1'"]=tmp; ACU2[" C1'"][0]=   6.741; ACU2[" C1'"][1]=   5.178; ACU2[" C1'"][2]=-142.017;
    ideal_rna["  A  C  U2"]=ACU2;
    map<string,vector<float> >().swap(ACU2);

    map<string,vector<float> >AGA0;
    AGA0[" P  "]=tmp; AGA0[" P  "][0]=   6.777; AGA0[" P  "][1]=   5.278; AGA0[" P  "][2]=-146.823;
    AGA0[" C4'"]=tmp; AGA0[" C4'"][0]=   9.316; AGA0[" C4'"][1]=   2.412; AGA0[" C4'"][2]=-145.946;
    AGA0[" C1'"]=tmp; AGA0[" C1'"][0]=   8.470; AGA0[" C1'"][1]=   0.715; AGA0[" C1'"][2]=-144.565;
    ideal_rna["  A  G  A0"]=AGA0;
    map<string,vector<float> >().swap(AGA0);

    map<string,vector<float> >AGA1;
    AGA1[" P  "]=tmp; AGA1[" P  "][0]=   8.555; AGA1[" P  "][1]=   0.780; AGA1[" P  "][2]=-149.371;
    AGA1[" C4'"]=tmp; AGA1[" C4'"][0]=   9.142; AGA1[" C4'"][1]=  -3.003; AGA1[" C4'"][2]=-148.494;
    AGA1[" C1'"]=tmp; AGA1[" C1'"][0]=   7.514; AGA1[" C1'"][1]=  -3.974; AGA1[" C1'"][2]=-147.113;
    ideal_rna["  A  G  A1"]=AGA1;
    map<string,vector<float> >().swap(AGA1);

    map<string,vector<float> >AGA2;
    AGA2[" P  "]=tmp; AGA2[" P  "][0]=   7.620; AGA2[" P  "][1]=  -3.965; AGA2[" P  "][2]=-151.919;
    AGA2[" C4'"]=tmp; AGA2[" C4'"][0]=   6.071; AGA2[" C4'"][1]=  -7.466; AGA2[" C4'"][2]=-151.042;
    AGA2[" C1'"]=tmp; AGA2[" C1'"][0]=   4.176; AGA2[" C1'"][1]=  -7.403; AGA2[" C1'"][2]=-149.661;
    ideal_rna["  A  G  A2"]=AGA2;
    map<string,vector<float> >().swap(AGA2);

    map<string,vector<float> >AGC0;
    AGC0[" P  "]=tmp; AGC0[" P  "][0]=   4.270; AGC0[" P  "][1]=  -7.453; AGC0[" P  "][2]=-154.467;
    AGC0[" C4'"]=tmp; AGC0[" C4'"][0]=   1.075; AGC0[" C4'"][1]=  -9.563; AGC0[" C4'"][2]=-153.590;
    AGC0[" C1'"]=tmp; AGC0[" C1'"][0]=  -0.485; AGC0[" C1'"][1]=  -8.486; AGC0[" C1'"][2]=-152.209;
    ideal_rna["  A  G  C0"]=AGC0;
    map<string,vector<float> >().swap(AGC0);

    map<string,vector<float> >AGC1;
    AGC1[" P  "]=tmp; AGC1[" P  "][0]=  -0.433; AGC1[" P  "][1]=  -8.579; AGC1[" P  "][2]=-157.015;
    AGC1[" C4'"]=tmp; AGC1[" C4'"][0]=  -4.261; AGC1[" C4'"][1]=  -8.628; AGC1[" C4'"][2]=-156.138;
    AGC1[" C1'"]=tmp; AGC1[" C1'"][0]=  -4.993; AGC1[" C1'"][1]=  -6.879; AGC1[" C1'"][2]=-154.757;
    ideal_rna["  A  G  C1"]=AGC1;
    map<string,vector<float> >().swap(AGC1);

    map<string,vector<float> >AGC2;
    AGC2[" P  "]=tmp; AGC2[" P  "][0]=  -4.999; AGC2[" P  "][1]=  -6.985; AGC2[" P  "][2]=-159.563;
    AGC2[" C4'"]=tmp; AGC2[" C4'"][0]=  -8.247; AGC2[" C4'"][1]=  -4.959; AGC2[" C4'"][2]=-158.686;
    AGC2[" C1'"]=tmp; AGC2[" C1'"][0]=  -7.918; AGC2[" C1'"][1]=  -3.092; AGC2[" C1'"][2]=-157.305;
    ideal_rna["  A  G  C2"]=AGC2;
    map<string,vector<float> >().swap(AGC2);

    map<string,vector<float> >AGG0;
    AGG0[" P  "]=tmp; AGG0[" P  "][0]=  -7.981; AGG0[" P  "][1]=  -3.177; AGG0[" P  "][2]=-162.111;
    AGG0[" C4'"]=tmp; AGG0[" C4'"][0]=  -9.619; AGG0[" C4'"][1]=   0.283; AGG0[" C4'"][2]=-161.234;
    AGG0[" C1'"]=tmp; AGG0[" C1'"][0]=  -8.333; AGG0[" C1'"][1]=   1.676; AGG0[" C1'"][2]=-159.853;
    ideal_rna["  A  G  G0"]=AGG0;
    map<string,vector<float> >().swap(AGG0);

    map<string,vector<float> >AGG1;
    AGG1[" P  "]=tmp; AGG1[" P  "][0]=  -8.432; AGG1[" P  "][1]=   1.638; AGG1[" P  "][2]=-164.659;
    AGG1[" C4'"]=tmp; AGG1[" C4'"][0]=  -7.942; AGG1[" C4'"][1]=   5.434; AGG1[" C4'"][2]=-163.782;
    AGG1[" C1'"]=tmp; AGG1[" C1'"][0]=  -6.107; AGG1[" C1'"][1]=   5.912; AGG1[" C1'"][2]=-162.401;
    ideal_rna["  A  G  G1"]=AGG1;
    map<string,vector<float> >().swap(AGG1);

    map<string,vector<float> >AGG2;
    AGG2[" P  "]=tmp; AGG2[" P  "][0]=  -6.211; AGG2[" P  "][1]=   5.934; AGG2[" P  "][2]=-167.207;
    AGG2[" C4'"]=tmp; AGG2[" C4'"][0]=  -3.747; AGG2[" C4'"][1]=   8.864; AGG2[" C4'"][2]=-166.330;
    AGG2[" C1'"]=tmp; AGG2[" C1'"][0]=  -1.945; AGG2[" C1'"][1]=   8.275; AGG2[" C1'"][2]=-164.949;
    ideal_rna["  A  G  G2"]=AGG2;
    map<string,vector<float> >().swap(AGG2);

    map<string,vector<float> >AGU0;
    AGU0[" P  "]=tmp; AGU0[" P  "][0]=  -2.021; AGU0[" P  "][1]=   8.349; AGU0[" P  "][2]=-169.755;
    AGU0[" C4'"]=tmp; AGU0[" C4'"][0]=   1.635; AGU0[" C4'"][1]=   9.483; AGU0[" C4'"][2]=-168.878;
    AGU0[" C1'"]=tmp; AGU0[" C1'"][0]=   2.834; AGU0[" C1'"][1]=   8.014; AGU0[" C1'"][2]=-167.497;
    ideal_rna["  A  G  U0"]=AGU0;
    map<string,vector<float> >().swap(AGU0);

    map<string,vector<float> >AGU1;
    AGU1[" P  "]=tmp; AGU1[" P  "][0]=   2.810; AGU1[" P  "][1]=   8.118; AGU1[" P  "][2]=-172.303;
    AGU1[" C4'"]=tmp; AGU1[" C4'"][0]=   6.499; AGU1[" C4'"][1]=   7.097; AGU1[" C4'"][2]=-171.426;
    AGU1[" C1'"]=tmp; AGU1[" C1'"][0]=   6.714; AGU1[" C1'"][1]=   5.213; AGU1[" C1'"][2]=-170.045;
    ideal_rna["  A  G  U1"]=AGU1;
    map<string,vector<float> >().swap(AGU1);

    map<string,vector<float> >AGU2;
    AGU2[" P  "]=tmp; AGU2[" P  "][0]=   6.750; AGU2[" P  "][1]=   5.313; AGU2[" P  "][2]=-174.851;
    AGU2[" C4'"]=tmp; AGU2[" C4'"][0]=   9.303; AGU2[" C4'"][1]=   2.461; AGU2[" C4'"][2]=-173.974;
    AGU2[" C1'"]=tmp; AGU2[" C1'"][0]=   8.466; AGU2[" C1'"][1]=   0.759; AGU2[" C1'"][2]=-172.593;
    ideal_rna["  A  G  U2"]=AGU2;
    map<string,vector<float> >().swap(AGU2);

    map<string,vector<float> >AUA0;
    AUA0[" P  "]=tmp; AUA0[" P  "][0]=   8.550; AUA0[" P  "][1]=   0.825; AUA0[" P  "][2]=-177.399;
    AUA0[" C4'"]=tmp; AUA0[" C4'"][0]=   9.158; AUA0[" C4'"][1]=  -2.955; AUA0[" C4'"][2]=-176.522;
    AUA0[" C1'"]=tmp; AUA0[" C1'"][0]=   7.535; AUA0[" C1'"][1]=  -3.935; AUA0[" C1'"][2]=-175.141;
    ideal_rna["  A  U  A0"]=AUA0;
    map<string,vector<float> >().swap(AUA0);

    map<string,vector<float> >AUA1;
    AUA1[" P  "]=tmp; AUA1[" P  "][0]=   7.641; AUA1[" P  "][1]=  -3.925; AUA1[" P  "][2]=-179.947;
    AUA1[" C4'"]=tmp; AUA1[" C4'"][0]=   6.110; AUA1[" C4'"][1]=  -7.434; AUA1[" C4'"][2]=-179.070;
    AUA1[" C1'"]=tmp; AUA1[" C1'"][0]=   4.215; AUA1[" C1'"][1]=  -7.381; AUA1[" C1'"][2]=-177.689;
    ideal_rna["  A  U  A1"]=AUA1;
    map<string,vector<float> >().swap(AUA1);

    map<string,vector<float> >AUA2;
    AUA2[" P  "]=tmp; AUA2[" P  "][0]=   4.309; AUA2[" P  "][1]=  -7.431; AUA2[" P  "][2]=-182.495;
    AUA2[" C4'"]=tmp; AUA2[" C4'"][0]=   1.125; AUA2[" C4'"][1]=  -9.557; AUA2[" C4'"][2]=-181.618;
    AUA2[" C1'"]=tmp; AUA2[" C1'"][0]=  -0.441; AUA2[" C1'"][1]=  -8.489; AUA2[" C1'"][2]=-180.237;
    ideal_rna["  A  U  A2"]=AUA2;
    map<string,vector<float> >().swap(AUA2);

    map<string,vector<float> >AUC0;
    AUC0[" P  "]=tmp; AUC0[" P  "][0]=  -0.388; AUC0[" P  "][1]=  -8.581; AUC0[" P  "][2]=-185.043;
    AUC0[" C4'"]=tmp; AUC0[" C4'"][0]=  -4.216; AUC0[" C4'"][1]=  -8.650; AUC0[" C4'"][2]=-184.166;
    AUC0[" C1'"]=tmp; AUC0[" C1'"][0]=  -4.957; AUC0[" C1'"][1]=  -6.905; AUC0[" C1'"][2]=-182.785;
    ideal_rna["  A  U  C0"]=AUC0;
    map<string,vector<float> >().swap(AUC0);

    map<string,vector<float> >AUC1;
    AUC1[" P  "]=tmp; AUC1[" P  "][0]=  -4.963; AUC1[" P  "][1]=  -7.011; AUC1[" P  "][2]=-187.591;
    AUC1[" C4'"]=tmp; AUC1[" C4'"][0]=  -8.221; AUC1[" C4'"][1]=  -5.002; AUC1[" C4'"][2]=-186.714;
    AUC1[" C1'"]=tmp; AUC1[" C1'"][0]=  -7.902; AUC1[" C1'"][1]=  -3.133; AUC1[" C1'"][2]=-185.333;
    ideal_rna["  A  U  C1"]=AUC1;
    map<string,vector<float> >().swap(AUC1);

    map<string,vector<float> >AUC2;
    AUC2[" P  "]=tmp; AUC2[" P  "][0]=  -7.964; AUC2[" P  "][1]=  -3.219; AUC2[" P  "][2]=-190.139;
    AUC2[" C4'"]=tmp; AUC2[" C4'"][0]=  -9.620; AUC2[" C4'"][1]=   0.232; AUC2[" C4'"][2]=-189.262;
    AUC2[" C1'"]=tmp; AUC2[" C1'"][0]=  -8.342; AUC2[" C1'"][1]=   1.632; AUC2[" C1'"][2]=-187.881;
    ideal_rna["  A  U  C2"]=AUC2;
    map<string,vector<float> >().swap(AUC2);

    map<string,vector<float> >AUG0;
    AUG0[" P  "]=tmp; AUG0[" P  "][0]=  -8.441; AUG0[" P  "][1]=   1.594; AUG0[" P  "][2]=-192.687;
    AUG0[" C4'"]=tmp; AUG0[" C4'"][0]=  -7.970; AUG0[" C4'"][1]=   5.393; AUG0[" C4'"][2]=-191.810;
    AUG0[" C1'"]=tmp; AUG0[" C1'"][0]=  -6.138; AUG0[" C1'"][1]=   5.880; AUG0[" C1'"][2]=-190.429;
    ideal_rna["  A  U  G0"]=AUG0;
    map<string,vector<float> >().swap(AUG0);

    map<string,vector<float> >AUG1;
    AUG1[" P  "]=tmp; AUG1[" P  "][0]=  -6.242; AUG1[" P  "][1]=   5.901; AUG1[" P  "][2]=-195.235;
    AUG1[" C4'"]=tmp; AUG1[" C4'"][0]=  -3.793; AUG1[" C4'"][1]=   8.844; AUG1[" C4'"][2]=-194.358;
    AUG1[" C1'"]=tmp; AUG1[" C1'"][0]=  -1.988; AUG1[" C1'"][1]=   8.264; AUG1[" C1'"][2]=-192.977;
    ideal_rna["  A  U  G1"]=AUG1;
    map<string,vector<float> >().swap(AUG1);

    map<string,vector<float> >AUG2;
    AUG2[" P  "]=tmp; AUG2[" P  "][0]=  -2.065; AUG2[" P  "][1]=   8.338; AUG2[" P  "][2]=-197.783;
    AUG2[" C4'"]=tmp; AUG2[" C4'"][0]=   1.586; AUG2[" C4'"][1]=   9.491; AUG2[" C4'"][2]=-196.906;
    AUG2[" C1'"]=tmp; AUG2[" C1'"][0]=   2.792; AUG2[" C1'"][1]=   8.028; AUG2[" C1'"][2]=-195.525;
    ideal_rna["  A  U  G2"]=AUG2;
    map<string,vector<float> >().swap(AUG2);

    map<string,vector<float> >AUU0;
    AUU0[" P  "]=tmp; AUU0[" P  "][0]=   2.767; AUU0[" P  "][1]=   8.132; AUU0[" P  "][2]=-200.331;
    AUU0[" C4'"]=tmp; AUU0[" C4'"][0]=   6.462; AUU0[" C4'"][1]=   7.131; AUU0[" C4'"][2]=-199.454;
    AUU0[" C1'"]=tmp; AUU0[" C1'"][0]=   6.687; AUU0[" C1'"][1]=   5.248; AUU0[" C1'"][2]=-198.073;
    ideal_rna["  A  U  U0"]=AUU0;
    map<string,vector<float> >().swap(AUU0);

    map<string,vector<float> >AUU1;
    AUU1[" P  "]=tmp; AUU1[" P  "][0]=   6.722; AUU1[" P  "][1]=   5.348; AUU1[" P  "][2]=-202.879;
    AUU1[" C4'"]=tmp; AUU1[" C4'"][0]=   9.290; AUU1[" C4'"][1]=   2.509; AUU1[" C4'"][2]=-202.002;
    AUU1[" C1'"]=tmp; AUU1[" C1'"][0]=   8.462; AUU1[" C1'"][1]=   0.804; AUU1[" C1'"][2]=-200.621;
    ideal_rna["  A  U  U1"]=AUU1;
    map<string,vector<float> >().swap(AUU1);

    map<string,vector<float> >AUU2;
    AUU2[" P  "]=tmp; AUU2[" P  "][0]=   8.546; AUU2[" P  "][1]=   0.869; AUU2[" P  "][2]=-205.427;
    AUU2[" C4'"]=tmp; AUU2[" C4'"][0]=   9.173; AUU2[" C4'"][1]=  -2.907; AUU2[" C4'"][2]=-204.550;
    AUU2[" C1'"]=tmp; AUU2[" C1'"][0]=   7.555; AUU2[" C1'"][1]=  -3.895; AUU2[" C1'"][2]=-203.169;
    ideal_rna["  A  U  U2"]=AUU2;
    map<string,vector<float> >().swap(AUU2);

    map<string,vector<float> >CAA0;
    CAA0[" P  "]=tmp; CAA0[" P  "][0]=   7.661; CAA0[" P  "][1]=  -3.885; CAA0[" P  "][2]=-207.975;
    CAA0[" C4'"]=tmp; CAA0[" C4'"][0]=   6.149; CAA0[" C4'"][1]=  -7.402; CAA0[" C4'"][2]=-207.098;
    CAA0[" C1'"]=tmp; CAA0[" C1'"][0]=   4.254; CAA0[" C1'"][1]=  -7.359; CAA0[" C1'"][2]=-205.717;
    ideal_rna["  C  A  A0"]=CAA0;
    map<string,vector<float> >().swap(CAA0);

    map<string,vector<float> >CAA1;
    CAA1[" P  "]=tmp; CAA1[" P  "][0]=   4.348; CAA1[" P  "][1]=  -7.408; CAA1[" P  "][2]=-210.523;
    CAA1[" C4'"]=tmp; CAA1[" C4'"][0]=   1.175; CAA1[" C4'"][1]=  -9.551; CAA1[" C4'"][2]=-209.646;
    CAA1[" C1'"]=tmp; CAA1[" C1'"][0]=  -0.397; CAA1[" C1'"][1]=  -8.491; CAA1[" C1'"][2]=-208.265;
    ideal_rna["  C  A  A1"]=CAA1;
    map<string,vector<float> >().swap(CAA1);

    map<string,vector<float> >CAA2;
    CAA2[" P  "]=tmp; CAA2[" P  "][0]=  -0.343; CAA2[" P  "][1]=  -8.583; CAA2[" P  "][2]=-213.071;
    CAA2[" C4'"]=tmp; CAA2[" C4'"][0]=  -4.171; CAA2[" C4'"][1]=  -8.672; CAA2[" C4'"][2]=-212.194;
    CAA2[" C1'"]=tmp; CAA2[" C1'"][0]=  -4.921; CAA2[" C1'"][1]=  -6.931; CAA2[" C1'"][2]=-210.813;
    ideal_rna["  C  A  A2"]=CAA2;
    map<string,vector<float> >().swap(CAA2);

    map<string,vector<float> >CAC0;
    CAC0[" P  "]=tmp; CAC0[" P  "][0]=  -4.926; CAC0[" P  "][1]=  -7.037; CAC0[" P  "][2]=-215.619;
    CAC0[" C4'"]=tmp; CAC0[" C4'"][0]=  -8.195; CAC0[" C4'"][1]=  -5.045; CAC0[" C4'"][2]=-214.742;
    CAC0[" C1'"]=tmp; CAC0[" C1'"][0]=  -7.885; CAC0[" C1'"][1]=  -3.175; CAC0[" C1'"][2]=-213.361;
    ideal_rna["  C  A  C0"]=CAC0;
    map<string,vector<float> >().swap(CAC0);

    map<string,vector<float> >CAC1;
    CAC1[" P  "]=tmp; CAC1[" P  "][0]=  -7.947; CAC1[" P  "][1]=  -3.261; CAC1[" P  "][2]=-218.167;
    CAC1[" C4'"]=tmp; CAC1[" C4'"][0]=  -9.621; CAC1[" C4'"][1]=   0.182; CAC1[" C4'"][2]=-217.290;
    CAC1[" C1'"]=tmp; CAC1[" C1'"][0]=  -8.350; CAC1[" C1'"][1]=   1.589; CAC1[" C1'"][2]=-215.909;
    ideal_rna["  C  A  C1"]=CAC1;
    map<string,vector<float> >().swap(CAC1);

    map<string,vector<float> >CAC2;
    CAC2[" P  "]=tmp; CAC2[" P  "][0]=  -8.449; CAC2[" P  "][1]=   1.549; CAC2[" P  "][2]=-220.715;
    CAC2[" C4'"]=tmp; CAC2[" C4'"][0]=  -7.998; CAC2[" C4'"][1]=   5.351; CAC2[" C4'"][2]=-219.838;
    CAC2[" C1'"]=tmp; CAC2[" C1'"][0]=  -6.169; CAC2[" C1'"][1]=   5.848; CAC2[" C1'"][2]=-218.457;
    ideal_rna["  C  A  C2"]=CAC2;
    map<string,vector<float> >().swap(CAC2);

    map<string,vector<float> >CAG0;
    CAG0[" P  "]=tmp; CAG0[" P  "][0]=  -6.273; CAG0[" P  "][1]=   5.868; CAG0[" P  "][2]=-223.263;
    CAG0[" C4'"]=tmp; CAG0[" C4'"][0]=  -3.840; CAG0[" C4'"][1]=   8.824; CAG0[" C4'"][2]=-222.386;
    CAG0[" C1'"]=tmp; CAG0[" C1'"][0]=  -2.032; CAG0[" C1'"][1]=   8.254; CAG0[" C1'"][2]=-221.005;
    ideal_rna["  C  A  G0"]=CAG0;
    map<string,vector<float> >().swap(CAG0);

    map<string,vector<float> >CAG1;
    CAG1[" P  "]=tmp; CAG1[" P  "][0]=  -2.108; CAG1[" P  "][1]=   8.327; CAG1[" P  "][2]=-225.811;
    CAG1[" C4'"]=tmp; CAG1[" C4'"][0]=   1.536; CAG1[" C4'"][1]=   9.500; CAG1[" C4'"][2]=-224.934;
    CAG1[" C1'"]=tmp; CAG1[" C1'"][0]=   2.750; CAG1[" C1'"][1]=   8.043; CAG1[" C1'"][2]=-223.553;
    ideal_rna["  C  A  G1"]=CAG1;
    map<string,vector<float> >().swap(CAG1);

    map<string,vector<float> >CAG2;
    CAG2[" P  "]=tmp; CAG2[" P  "][0]=   2.724; CAG2[" P  "][1]=   8.147; CAG2[" P  "][2]=-228.359;
    CAG2[" C4'"]=tmp; CAG2[" C4'"][0]=   6.425; CAG2[" C4'"][1]=   7.164; CAG2[" C4'"][2]=-227.482;
    CAG2[" C1'"]=tmp; CAG2[" C1'"][0]=   6.659; CAG2[" C1'"][1]=   5.283; CAG2[" C1'"][2]=-226.101;
    ideal_rna["  C  A  G2"]=CAG2;
    map<string,vector<float> >().swap(CAG2);

    map<string,vector<float> >CAU0;
    CAU0[" P  "]=tmp; CAU0[" P  "][0]=   6.694; CAU0[" P  "][1]=   5.384; CAU0[" P  "][2]=-230.907;
    CAU0[" C4'"]=tmp; CAU0[" C4'"][0]=   9.277; CAU0[" C4'"][1]=   2.558; CAU0[" C4'"][2]=-230.030;
    CAU0[" C1'"]=tmp; CAU0[" C1'"][0]=   8.458; CAU0[" C1'"][1]=   0.849; CAU0[" C1'"][2]=-228.649;
    ideal_rna["  C  A  U0"]=CAU0;
    map<string,vector<float> >().swap(CAU0);

    map<string,vector<float> >CAU1;
    CAU1[" P  "]=tmp; CAU1[" P  "][0]=   8.541; CAU1[" P  "][1]=   0.914; CAU1[" P  "][2]=-233.455;
    CAU1[" C4'"]=tmp; CAU1[" C4'"][0]=   9.188; CAU1[" C4'"][1]=  -2.859; CAU1[" C4'"][2]=-232.578;
    CAU1[" C1'"]=tmp; CAU1[" C1'"][0]=   7.575; CAU1[" C1'"][1]=  -3.855; CAU1[" C1'"][2]=-231.197;
    ideal_rna["  C  A  U1"]=CAU1;
    map<string,vector<float> >().swap(CAU1);

    map<string,vector<float> >CAU2;
    CAU2[" P  "]=tmp; CAU2[" P  "][0]=   7.681; CAU2[" P  "][1]=  -3.845; CAU2[" P  "][2]=-236.003;
    CAU2[" C4'"]=tmp; CAU2[" C4'"][0]=   6.188; CAU2[" C4'"][1]=  -7.370; CAU2[" C4'"][2]=-235.126;
    CAU2[" C1'"]=tmp; CAU2[" C1'"][0]=   4.292; CAU2[" C1'"][1]=  -7.337; CAU2[" C1'"][2]=-233.745;
    ideal_rna["  C  A  U2"]=CAU2;
    map<string,vector<float> >().swap(CAU2);

    map<string,vector<float> >CCA0;
    CCA0[" P  "]=tmp; CCA0[" P  "][0]=   4.387; CCA0[" P  "][1]=  -7.385; CCA0[" P  "][2]=-238.551;
    CCA0[" C4'"]=tmp; CCA0[" C4'"][0]=   1.225; CCA0[" C4'"][1]=  -9.545; CCA0[" C4'"][2]=-237.674;
    CCA0[" C1'"]=tmp; CCA0[" C1'"][0]=  -0.351; CCA0[" C1'"][1]=  -8.493; CCA0[" C1'"][2]=-236.293;
    ideal_rna["  C  C  A0"]=CCA0;
    map<string,vector<float> >().swap(CCA0);

    map<string,vector<float> >CCA1;
    CCA1[" P  "]=tmp; CCA1[" P  "][0]=  -0.298; CCA1[" P  "][1]=  -8.585; CCA1[" P  "][2]=-241.099;
    CCA1[" C4'"]=tmp; CCA1[" C4'"][0]=  -4.125; CCA1[" C4'"][1]=  -8.694; CCA1[" C4'"][2]=-240.222;
    CCA1[" C1'"]=tmp; CCA1[" C1'"][0]=  -4.884; CCA1[" C1'"][1]=  -6.957; CCA1[" C1'"][2]=-238.841;
    ideal_rna["  C  C  A1"]=CCA1;
    map<string,vector<float> >().swap(CCA1);

    map<string,vector<float> >CCA2;
    CCA2[" P  "]=tmp; CCA2[" P  "][0]=  -4.889; CCA2[" P  "][1]=  -7.063; CCA2[" P  "][2]=-243.647;
    CCA2[" C4'"]=tmp; CCA2[" C4'"][0]=  -8.168; CCA2[" C4'"][1]=  -5.087; CCA2[" C4'"][2]=-242.770;
    CCA2[" C1'"]=tmp; CCA2[" C1'"][0]=  -7.868; CCA2[" C1'"][1]=  -3.215; CCA2[" C1'"][2]=-241.389;
    ideal_rna["  C  C  A2"]=CCA2;
    map<string,vector<float> >().swap(CCA2);

    map<string,vector<float> >CCC0;
    CCC0[" P  "]=tmp; CCC0[" P  "][0]=  -7.930; CCC0[" P  "][1]=  -3.302; CCC0[" P  "][2]=-246.195;
    CCC0[" C4'"]=tmp; CCC0[" C4'"][0]=  -9.622; CCC0[" C4'"][1]=   0.132; CCC0[" C4'"][2]=-245.318;
    CCC0[" C1'"]=tmp; CCC0[" C1'"][0]=  -8.359; CCC0[" C1'"][1]=   1.544; CCC0[" C1'"][2]=-243.937;
    ideal_rna["  C  C  C0"]=CCC0;
    map<string,vector<float> >().swap(CCC0);

    map<string,vector<float> >CCC1;
    CCC1[" P  "]=tmp; CCC1[" P  "][0]=  -8.457; CCC1[" P  "][1]=   1.505; CCC1[" P  "][2]=-248.743;
    CCC1[" C4'"]=tmp; CCC1[" C4'"][0]=  -8.026; CCC1[" C4'"][1]=   5.309; CCC1[" C4'"][2]=-247.866;
    CCC1[" C1'"]=tmp; CCC1[" C1'"][0]=  -6.199; CCC1[" C1'"][1]=   5.815; CCC1[" C1'"][2]=-246.485;
    ideal_rna["  C  C  C1"]=CCC1;
    map<string,vector<float> >().swap(CCC1);

    map<string,vector<float> >CCC2;
    CCC2[" P  "]=tmp; CCC2[" P  "][0]=  -6.304; CCC2[" P  "][1]=   5.835; CCC2[" P  "][2]=-251.291;
    CCC2[" C4'"]=tmp; CCC2[" C4'"][0]=  -3.886; CCC2[" C4'"][1]=   8.804; CCC2[" C4'"][2]=-250.414;
    CCC2[" C1'"]=tmp; CCC2[" C1'"][0]=  -2.075; CCC2[" C1'"][1]=   8.243; CCC2[" C1'"][2]=-249.033;
    ideal_rna["  C  C  C2"]=CCC2;
    map<string,vector<float> >().swap(CCC2);

    map<string,vector<float> >CCG0;
    CCG0[" P  "]=tmp; CCG0[" P  "][0]=  -2.152; CCG0[" P  "][1]=   8.316; CCG0[" P  "][2]=-253.839;
    CCG0[" C4'"]=tmp; CCG0[" C4'"][0]=   1.486; CCG0[" C4'"][1]=   9.508; CCG0[" C4'"][2]=-252.962;
    CCG0[" C1'"]=tmp; CCG0[" C1'"][0]=   2.707; CCG0[" C1'"][1]=   8.057; CCG0[" C1'"][2]=-251.581;
    ideal_rna["  C  C  G0"]=CCG0;
    map<string,vector<float> >().swap(CCG0);

    map<string,vector<float> >CCG1;
    CCG1[" P  "]=tmp; CCG1[" P  "][0]=   2.682; CCG1[" P  "][1]=   8.161; CCG1[" P  "][2]=-256.387;
    CCG1[" C4'"]=tmp; CCG1[" C4'"][0]=   6.387; CCG1[" C4'"][1]=   7.198; CCG1[" C4'"][2]=-255.510;
    CCG1[" C1'"]=tmp; CCG1[" C1'"][0]=   6.631; CCG1[" C1'"][1]=   5.318; CCG1[" C1'"][2]=-254.129;
    ideal_rna["  C  C  G1"]=CCG1;
    map<string,vector<float> >().swap(CCG1);

    map<string,vector<float> >CCG2;
    CCG2[" P  "]=tmp; CCG2[" P  "][0]=   6.665; CCG2[" P  "][1]=   5.419; CCG2[" P  "][2]=-258.935;
    CCG2[" C4'"]=tmp; CCG2[" C4'"][0]=   9.263; CCG2[" C4'"][1]=   2.607; CCG2[" C4'"][2]=-258.058;
    CCG2[" C1'"]=tmp; CCG2[" C1'"][0]=   8.453; CCG2[" C1'"][1]=   0.892; CCG2[" C1'"][2]=-256.677;
    ideal_rna["  C  C  G2"]=CCG2;
    map<string,vector<float> >().swap(CCG2);

    map<string,vector<float> >CCU0;
    CCU0[" P  "]=tmp; CCU0[" P  "][0]=   8.536; CCU0[" P  "][1]=   0.959; CCU0[" P  "][2]=-261.483;
    CCU0[" C4'"]=tmp; CCU0[" C4'"][0]=   9.203; CCU0[" C4'"][1]=  -2.811; CCU0[" C4'"][2]=-260.606;
    CCU0[" C1'"]=tmp; CCU0[" C1'"][0]=   7.596; CCU0[" C1'"][1]=  -3.815; CCU0[" C1'"][2]=-259.225;
    ideal_rna["  C  C  U0"]=CCU0;
    map<string,vector<float> >().swap(CCU0);

    map<string,vector<float> >CCU1;
    CCU1[" P  "]=tmp; CCU1[" P  "][0]=   7.701; CCU1[" P  "][1]=  -3.805; CCU1[" P  "][2]=-264.031;
    CCU1[" C4'"]=tmp; CCU1[" C4'"][0]=   6.226; CCU1[" C4'"][1]=  -7.337; CCU1[" C4'"][2]=-263.154;
    CCU1[" C1'"]=tmp; CCU1[" C1'"][0]=   4.331; CCU1[" C1'"][1]=  -7.314; CCU1[" C1'"][2]=-261.773;
    ideal_rna["  C  C  U1"]=CCU1;
    map<string,vector<float> >().swap(CCU1);

    map<string,vector<float> >CCU2;
    CCU2[" P  "]=tmp; CCU2[" P  "][0]=   4.425; CCU2[" P  "][1]=  -7.362; CCU2[" P  "][2]=-266.579;
    CCU2[" C4'"]=tmp; CCU2[" C4'"][0]=   1.275; CCU2[" C4'"][1]=  -9.538; CCU2[" C4'"][2]=-265.702;
    CCU2[" C1'"]=tmp; CCU2[" C1'"][0]=  -0.308; CCU2[" C1'"][1]=  -8.494; CCU2[" C1'"][2]=-264.321;
    ideal_rna["  C  C  U2"]=CCU2;
    map<string,vector<float> >().swap(CCU2);

    map<string,vector<float> >CGA0;
    CGA0[" P  "]=tmp; CGA0[" P  "][0]=  -0.253; CGA0[" P  "][1]=  -8.586; CGA0[" P  "][2]=-269.127;
    CGA0[" C4'"]=tmp; CGA0[" C4'"][0]=  -4.080; CGA0[" C4'"][1]=  -8.715; CGA0[" C4'"][2]=-268.250;
    CGA0[" C1'"]=tmp; CGA0[" C1'"][0]=  -4.847; CGA0[" C1'"][1]=  -6.982; CGA0[" C1'"][2]=-266.869;
    ideal_rna["  C  G  A0"]=CGA0;
    map<string,vector<float> >().swap(CGA0);

    map<string,vector<float> >CGA1;
    CGA1[" P  "]=tmp; CGA1[" P  "][0]=  -4.852; CGA1[" P  "][1]=  -7.088; CGA1[" P  "][2]=-271.675;
    CGA1[" C4'"]=tmp; CGA1[" C4'"][0]=  -8.141; CGA1[" C4'"][1]=  -5.130; CGA1[" C4'"][2]=-270.798;
    CGA1[" C1'"]=tmp; CGA1[" C1'"][0]=  -7.851; CGA1[" C1'"][1]=  -3.256; CGA1[" C1'"][2]=-269.417;
    ideal_rna["  C  G  A1"]=CGA1;
    map<string,vector<float> >().swap(CGA1);

    map<string,vector<float> >CGA2;
    CGA2[" P  "]=tmp; CGA2[" P  "][0]=  -7.912; CGA2[" P  "][1]=  -3.344; CGA2[" P  "][2]=-274.223;
    CGA2[" C4'"]=tmp; CGA2[" C4'"][0]=  -9.623; CGA2[" C4'"][1]=   0.081; CGA2[" C4'"][2]=-273.346;
    CGA2[" C1'"]=tmp; CGA2[" C1'"][0]=  -8.366; CGA2[" C1'"][1]=   1.501; CGA2[" C1'"][2]=-271.965;
    ideal_rna["  C  G  A2"]=CGA2;
    map<string,vector<float> >().swap(CGA2);

    map<string,vector<float> >CGC0;
    CGC0[" P  "]=tmp; CGC0[" P  "][0]=  -8.465; CGC0[" P  "][1]=   1.461; CGC0[" P  "][2]=-276.771;
    CGC0[" C4'"]=tmp; CGC0[" C4'"][0]=  -8.054; CGC0[" C4'"][1]=   5.267; CGC0[" C4'"][2]=-275.894;
    CGC0[" C1'"]=tmp; CGC0[" C1'"][0]=  -6.230; CGC0[" C1'"][1]=   5.783; CGC0[" C1'"][2]=-274.513;
    ideal_rna["  C  G  C0"]=CGC0;
    map<string,vector<float> >().swap(CGC0);

    map<string,vector<float> >CGC1;
    CGC1[" P  "]=tmp; CGC1[" P  "][0]=  -6.334; CGC1[" P  "][1]=   5.802; CGC1[" P  "][2]=-279.319;
    CGC1[" C4'"]=tmp; CGC1[" C4'"][0]=  -3.932; CGC1[" C4'"][1]=   8.783; CGC1[" C4'"][2]=-278.442;
    CGC1[" C1'"]=tmp; CGC1[" C1'"][0]=  -2.118; CGC1[" C1'"][1]=   8.232; CGC1[" C1'"][2]=-277.061;
    ideal_rna["  C  G  C1"]=CGC1;
    map<string,vector<float> >().swap(CGC1);

    map<string,vector<float> >CGC2;
    CGC2[" P  "]=tmp; CGC2[" P  "][0]=  -2.196; CGC2[" P  "][1]=   8.305; CGC2[" P  "][2]=-281.867;
    CGC2[" C4'"]=tmp; CGC2[" C4'"][0]=   1.436; CGC2[" C4'"][1]=   9.515; CGC2[" C4'"][2]=-280.990;
    CGC2[" C1'"]=tmp; CGC2[" C1'"][0]=   2.665; CGC2[" C1'"][1]=   8.072; CGC2[" C1'"][2]=-279.609;
    ideal_rna["  C  G  C2"]=CGC2;
    map<string,vector<float> >().swap(CGC2);

    map<string,vector<float> >CGG0;
    CGG0[" P  "]=tmp; CGG0[" P  "][0]=   2.639; CGG0[" P  "][1]=   8.175; CGG0[" P  "][2]=-284.415;
    CGG0[" C4'"]=tmp; CGG0[" C4'"][0]=   6.349; CGG0[" C4'"][1]=   7.231; CGG0[" C4'"][2]=-283.538;
    CGG0[" C1'"]=tmp; CGG0[" C1'"][0]=   6.603; CGG0[" C1'"][1]=   5.353; CGG0[" C1'"][2]=-282.157;
    ideal_rna["  C  G  G0"]=CGG0;
    map<string,vector<float> >().swap(CGG0);

    map<string,vector<float> >CGG1;
    CGG1[" P  "]=tmp; CGG1[" P  "][0]=   6.637; CGG1[" P  "][1]=   5.453; CGG1[" P  "][2]=-286.963;
    CGG1[" C4'"]=tmp; CGG1[" C4'"][0]=   9.249; CGG1[" C4'"][1]=   2.655; CGG1[" C4'"][2]=-286.086;
    CGG1[" C1'"]=tmp; CGG1[" C1'"][0]=   8.448; CGG1[" C1'"][1]=   0.937; CGG1[" C1'"][2]=-284.705;
    ideal_rna["  C  G  G1"]=CGG1;
    map<string,vector<float> >().swap(CGG1);

    map<string,vector<float> >CGG2;
    CGG2[" P  "]=tmp; CGG2[" P  "][0]=   8.531; CGG2[" P  "][1]=   1.004; CGG2[" P  "][2]=-289.511;
    CGG2[" C4'"]=tmp; CGG2[" C4'"][0]=   9.218; CGG2[" C4'"][1]=  -2.763; CGG2[" C4'"][2]=-288.634;
    CGG2[" C1'"]=tmp; CGG2[" C1'"][0]=   7.615; CGG2[" C1'"][1]=  -3.776; CGG2[" C1'"][2]=-287.253;
    ideal_rna["  C  G  G2"]=CGG2;
    map<string,vector<float> >().swap(CGG2);

    map<string,vector<float> >CGU0;
    CGU0[" P  "]=tmp; CGU0[" P  "][0]=   7.721; CGU0[" P  "][1]=  -3.764; CGU0[" P  "][2]=-292.059;
    CGU0[" C4'"]=tmp; CGU0[" C4'"][0]=   6.264; CGU0[" C4'"][1]=  -7.305; CGU0[" C4'"][2]=-291.182;
    CGU0[" C1'"]=tmp; CGU0[" C1'"][0]=   4.369; CGU0[" C1'"][1]=  -7.291; CGU0[" C1'"][2]=-289.801;
    ideal_rna["  C  G  U0"]=CGU0;
    map<string,vector<float> >().swap(CGU0);

    map<string,vector<float> >CGU1;
    CGU1[" P  "]=tmp; CGU1[" P  "][0]=   4.464; CGU1[" P  "][1]=  -7.339; CGU1[" P  "][2]=-294.607;
    CGU1[" C4'"]=tmp; CGU1[" C4'"][0]=   1.325; CGU1[" C4'"][1]=  -9.531; CGU1[" C4'"][2]=-293.730;
    CGU1[" C1'"]=tmp; CGU1[" C1'"][0]=  -0.263; CGU1[" C1'"][1]=  -8.496; CGU1[" C1'"][2]=-292.349;
    ideal_rna["  C  G  U1"]=CGU1;
    map<string,vector<float> >().swap(CGU1);

    map<string,vector<float> >CGU2;
    CGU2[" P  "]=tmp; CGU2[" P  "][0]=  -0.209; CGU2[" P  "][1]=  -8.587; CGU2[" P  "][2]=-297.155;
    CGU2[" C4'"]=tmp; CGU2[" C4'"][0]=  -4.034; CGU2[" C4'"][1]=  -8.737; CGU2[" C4'"][2]=-296.278;
    CGU2[" C1'"]=tmp; CGU2[" C1'"][0]=  -4.811; CGU2[" C1'"][1]=  -7.007; CGU2[" C1'"][2]=-294.897;
    ideal_rna["  C  G  U2"]=CGU2;
    map<string,vector<float> >().swap(CGU2);

    map<string,vector<float> >CUA0;
    CUA0[" P  "]=tmp; CUA0[" P  "][0]=  -4.815; CUA0[" P  "][1]=  -7.114; CUA0[" P  "][2]=-299.703;
    CUA0[" C4'"]=tmp; CUA0[" C4'"][0]=  -8.115; CUA0[" C4'"][1]=  -5.173; CUA0[" C4'"][2]=-298.826;
    CUA0[" C1'"]=tmp; CUA0[" C1'"][0]=  -7.834; CUA0[" C1'"][1]=  -3.298; CUA0[" C1'"][2]=-297.445;
    ideal_rna["  C  U  A0"]=CUA0;
    map<string,vector<float> >().swap(CUA0);

    map<string,vector<float> >CUA1;
    CUA1[" P  "]=tmp; CUA1[" P  "][0]=  -7.895; CUA1[" P  "][1]=  -3.385; CUA1[" P  "][2]=-302.251;
    CUA1[" C4'"]=tmp; CUA1[" C4'"][0]=  -9.623; CUA1[" C4'"][1]=   0.031; CUA1[" C4'"][2]=-301.374;
    CUA1[" C1'"]=tmp; CUA1[" C1'"][0]=  -8.374; CUA1[" C1'"][1]=   1.458; CUA1[" C1'"][2]=-299.993;
    ideal_rna["  C  U  A1"]=CUA1;
    map<string,vector<float> >().swap(CUA1);

    map<string,vector<float> >CUA2;
    CUA2[" P  "]=tmp; CUA2[" P  "][0]=  -8.472; CUA2[" P  "][1]=   1.416; CUA2[" P  "][2]=-304.799;
    CUA2[" C4'"]=tmp; CUA2[" C4'"][0]=  -8.081; CUA2[" C4'"][1]=   5.225; CUA2[" C4'"][2]=-303.922;
    CUA2[" C1'"]=tmp; CUA2[" C1'"][0]=  -6.259; CUA2[" C1'"][1]=   5.751; CUA2[" C1'"][2]=-302.541;
    ideal_rna["  C  U  A2"]=CUA2;
    map<string,vector<float> >().swap(CUA2);

    map<string,vector<float> >CUC0;
    CUC0[" P  "]=tmp; CUC0[" P  "][0]=  -6.364; CUC0[" P  "][1]=   5.769; CUC0[" P  "][2]=-307.347;
    CUC0[" C4'"]=tmp; CUC0[" C4'"][0]=  -3.978; CUC0[" C4'"][1]=   8.762; CUC0[" C4'"][2]=-306.470;
    CUC0[" C1'"]=tmp; CUC0[" C1'"][0]=  -2.161; CUC0[" C1'"][1]=   8.221; CUC0[" C1'"][2]=-305.089;
    ideal_rna["  C  U  C0"]=CUC0;
    map<string,vector<float> >().swap(CUC0);

    map<string,vector<float> >CUC1;
    CUC1[" P  "]=tmp; CUC1[" P  "][0]=  -2.239; CUC1[" P  "][1]=   8.293; CUC1[" P  "][2]=-309.895;
    CUC1[" C4'"]=tmp; CUC1[" C4'"][0]=   1.386; CUC1[" C4'"][1]=   9.523; CUC1[" C4'"][2]=-309.018;
    CUC1[" C1'"]=tmp; CUC1[" C1'"][0]=   2.623; CUC1[" C1'"][1]=   8.085; CUC1[" C1'"][2]=-307.637;
    ideal_rna["  C  U  C1"]=CUC1;
    map<string,vector<float> >().swap(CUC1);

    map<string,vector<float> >CUC2;
    CUC2[" P  "]=tmp; CUC2[" P  "][0]=   2.596; CUC2[" P  "][1]=   8.188; CUC2[" P  "][2]=-312.443;
    CUC2[" C4'"]=tmp; CUC2[" C4'"][0]=   6.311; CUC2[" C4'"][1]=   7.264; CUC2[" C4'"][2]=-311.566;
    CUC2[" C1'"]=tmp; CUC2[" C1'"][0]=   6.575; CUC2[" C1'"][1]=   5.387; CUC2[" C1'"][2]=-310.185;
    ideal_rna["  C  U  C2"]=CUC2;
    map<string,vector<float> >().swap(CUC2);

    map<string,vector<float> >CUG0;
    CUG0[" P  "]=tmp; CUG0[" P  "][0]=   6.608; CUG0[" P  "][1]=   5.488; CUG0[" P  "][2]=-314.991;
    CUG0[" C4'"]=tmp; CUG0[" C4'"][0]=   9.235; CUG0[" C4'"][1]=   2.703; CUG0[" C4'"][2]=-314.114;
    CUG0[" C1'"]=tmp; CUG0[" C1'"][0]=   8.443; CUG0[" C1'"][1]=   0.982; CUG0[" C1'"][2]=-312.733;
    ideal_rna["  C  U  G0"]=CUG0;
    map<string,vector<float> >().swap(CUG0);

    map<string,vector<float> >CUG1;
    CUG1[" P  "]=tmp; CUG1[" P  "][0]=   8.526; CUG1[" P  "][1]=   1.048; CUG1[" P  "][2]=-317.539;
    CUG1[" C4'"]=tmp; CUG1[" C4'"][0]=   9.232; CUG1[" C4'"][1]=  -2.714; CUG1[" C4'"][2]=-316.662;
    CUG1[" C1'"]=tmp; CUG1[" C1'"][0]=   7.635; CUG1[" C1'"][1]=  -3.736; CUG1[" C1'"][2]=-315.281;
    ideal_rna["  C  U  G1"]=CUG1;
    map<string,vector<float> >().swap(CUG1);

    map<string,vector<float> >CUG2;
    CUG2[" P  "]=tmp; CUG2[" P  "][0]=   7.741; CUG2[" P  "][1]=  -3.724; CUG2[" P  "][2]=-320.087;
    CUG2[" C4'"]=tmp; CUG2[" C4'"][0]=   6.303; CUG2[" C4'"][1]=  -7.272; CUG2[" C4'"][2]=-319.210;
    CUG2[" C1'"]=tmp; CUG2[" C1'"][0]=   4.407; CUG2[" C1'"][1]=  -7.269; CUG2[" C1'"][2]=-317.829;
    ideal_rna["  C  U  G2"]=CUG2;
    map<string,vector<float> >().swap(CUG2);

    map<string,vector<float> >CUU0;
    CUU0[" P  "]=tmp; CUU0[" P  "][0]=   4.502; CUU0[" P  "][1]=  -7.316; CUU0[" P  "][2]=-322.635;
    CUU0[" C4'"]=tmp; CUU0[" C4'"][0]=   1.375; CUU0[" C4'"][1]=  -9.524; CUU0[" C4'"][2]=-321.758;
    CUU0[" C1'"]=tmp; CUU0[" C1'"][0]=  -0.218; CUU0[" C1'"][1]=  -8.497; CUU0[" C1'"][2]=-320.377;
    ideal_rna["  C  U  U0"]=CUU0;
    map<string,vector<float> >().swap(CUU0);

    map<string,vector<float> >CUU1;
    CUU1[" P  "]=tmp; CUU1[" P  "][0]=  -0.164; CUU1[" P  "][1]=  -8.588; CUU1[" P  "][2]=-325.183;
    CUU1[" C4'"]=tmp; CUU1[" C4'"][0]=  -3.988; CUU1[" C4'"][1]=  -8.758; CUU1[" C4'"][2]=-324.306;
    CUU1[" C1'"]=tmp; CUU1[" C1'"][0]=  -4.775; CUU1[" C1'"][1]=  -7.032; CUU1[" C1'"][2]=-322.925;
    ideal_rna["  C  U  U1"]=CUU1;
    map<string,vector<float> >().swap(CUU1);

    map<string,vector<float> >CUU2;
    CUU2[" P  "]=tmp; CUU2[" P  "][0]=  -4.777; CUU2[" P  "][1]=  -7.139; CUU2[" P  "][2]=-327.731;
    CUU2[" C4'"]=tmp; CUU2[" C4'"][0]=  -8.087; CUU2[" C4'"][1]=  -5.215; CUU2[" C4'"][2]=-326.854;
    CUU2[" C1'"]=tmp; CUU2[" C1'"][0]=  -7.817; CUU2[" C1'"][1]=  -3.338; CUU2[" C1'"][2]=-325.473;
    ideal_rna["  C  U  U2"]=CUU2;
    map<string,vector<float> >().swap(CUU2);

    map<string,vector<float> >GAA0;
    GAA0[" P  "]=tmp; GAA0[" P  "][0]=  -7.877; GAA0[" P  "][1]=  -3.426; GAA0[" P  "][2]=-330.279;
    GAA0[" C4'"]=tmp; GAA0[" C4'"][0]=  -9.623; GAA0[" C4'"][1]=  -0.019; GAA0[" C4'"][2]=-329.402;
    GAA0[" C1'"]=tmp; GAA0[" C1'"][0]=  -8.382; GAA0[" C1'"][1]=   1.414; GAA0[" C1'"][2]=-328.021;
    ideal_rna["  G  A  A0"]=GAA0;
    map<string,vector<float> >().swap(GAA0);

    map<string,vector<float> >GAA1;
    GAA1[" P  "]=tmp; GAA1[" P  "][0]=  -8.480; GAA1[" P  "][1]=   1.372; GAA1[" P  "][2]=-332.827;
    GAA1[" C4'"]=tmp; GAA1[" C4'"][0]=  -8.108; GAA1[" C4'"][1]=   5.182; GAA1[" C4'"][2]=-331.950;
    GAA1[" C1'"]=tmp; GAA1[" C1'"][0]=  -6.289; GAA1[" C1'"][1]=   5.718; GAA1[" C1'"][2]=-330.569;
    ideal_rna["  G  A  A1"]=GAA1;
    map<string,vector<float> >().swap(GAA1);

    map<string,vector<float> >GAA2;
    GAA2[" P  "]=tmp; GAA2[" P  "][0]=  -6.395; GAA2[" P  "][1]=   5.736; GAA2[" P  "][2]=-335.375;
    GAA2[" C4'"]=tmp; GAA2[" C4'"][0]=  -4.024; GAA2[" C4'"][1]=   8.741; GAA2[" C4'"][2]=-334.498;
    GAA2[" C1'"]=tmp; GAA2[" C1'"][0]=  -2.204; GAA2[" C1'"][1]=   8.209; GAA2[" C1'"][2]=-333.117;
    ideal_rna["  G  A  A2"]=GAA2;
    map<string,vector<float> >().swap(GAA2);

    map<string,vector<float> >GAC0;
    GAC0[" P  "]=tmp; GAC0[" P  "][0]=  -2.282; GAC0[" P  "][1]=   8.281; GAC0[" P  "][2]=-337.923;
    GAC0[" C4'"]=tmp; GAC0[" C4'"][0]=   1.337; GAC0[" C4'"][1]=   9.530; GAC0[" C4'"][2]=-337.046;
    GAC0[" C1'"]=tmp; GAC0[" C1'"][0]=   2.581; GAC0[" C1'"][1]=   8.099; GAC0[" C1'"][2]=-335.665;
    ideal_rna["  G  A  C0"]=GAC0;
    map<string,vector<float> >().swap(GAC0);

    map<string,vector<float> >GAC1;
    GAC1[" P  "]=tmp; GAC1[" P  "][0]=   2.553; GAC1[" P  "][1]=   8.202; GAC1[" P  "][2]=-340.471;
    GAC1[" C4'"]=tmp; GAC1[" C4'"][0]=   6.273; GAC1[" C4'"][1]=   7.297; GAC1[" C4'"][2]=-339.594;
    GAC1[" C1'"]=tmp; GAC1[" C1'"][0]=   6.547; GAC1[" C1'"][1]=   5.421; GAC1[" C1'"][2]=-338.213;
    ideal_rna["  G  A  C1"]=GAC1;
    map<string,vector<float> >().swap(GAC1);

    map<string,vector<float> >GAC2;
    GAC2[" P  "]=tmp; GAC2[" P  "][0]=   6.579; GAC2[" P  "][1]=   5.523; GAC2[" P  "][2]=-343.019;
    GAC2[" C4'"]=tmp; GAC2[" C4'"][0]=   9.221; GAC2[" C4'"][1]=   2.752; GAC2[" C4'"][2]=-342.142;
    GAC2[" C1'"]=tmp; GAC2[" C1'"][0]=   8.438; GAC2[" C1'"][1]=   1.026; GAC2[" C1'"][2]=-340.761;
    ideal_rna["  G  A  C2"]=GAC2;
    map<string,vector<float> >().swap(GAC2);

    map<string,vector<float> >GAG0;
    GAG0[" P  "]=tmp; GAG0[" P  "][0]=   8.520; GAG0[" P  "][1]=   1.093; GAG0[" P  "][2]=-345.567;
    GAG0[" C4'"]=tmp; GAG0[" C4'"][0]=   9.246; GAG0[" C4'"][1]=  -2.666; GAG0[" C4'"][2]=-344.690;
    GAG0[" C1'"]=tmp; GAG0[" C1'"][0]=   7.654; GAG0[" C1'"][1]=  -3.696; GAG0[" C1'"][2]=-343.309;
    ideal_rna["  G  A  G0"]=GAG0;
    map<string,vector<float> >().swap(GAG0);

    map<string,vector<float> >GAG1;
    GAG1[" P  "]=tmp; GAG1[" P  "][0]=   7.760; GAG1[" P  "][1]=  -3.683; GAG1[" P  "][2]=-348.115;
    GAG1[" C4'"]=tmp; GAG1[" C4'"][0]=   6.341; GAG1[" C4'"][1]=  -7.239; GAG1[" C4'"][2]=-347.238;
    GAG1[" C1'"]=tmp; GAG1[" C1'"][0]=   4.445; GAG1[" C1'"][1]=  -7.245; GAG1[" C1'"][2]=-345.857;
    ideal_rna["  G  A  G1"]=GAG1;
    map<string,vector<float> >().swap(GAG1);

    map<string,vector<float> >GAG2;
    GAG2[" P  "]=tmp; GAG2[" P  "][0]=   4.540; GAG2[" P  "][1]=  -7.292; GAG2[" P  "][2]=-350.663;
    GAG2[" C4'"]=tmp; GAG2[" C4'"][0]=   1.425; GAG2[" C4'"][1]=  -9.517; GAG2[" C4'"][2]=-349.786;
    GAG2[" C1'"]=tmp; GAG2[" C1'"][0]=  -0.174; GAG2[" C1'"][1]=  -8.498; GAG2[" C1'"][2]=-348.405;
    ideal_rna["  G  A  G2"]=GAG2;
    map<string,vector<float> >().swap(GAG2);

    map<string,vector<float> >GAU0;
    GAU0[" P  "]=tmp; GAU0[" P  "][0]=  -0.119; GAU0[" P  "][1]=  -8.589; GAU0[" P  "][2]=-353.211;
    GAU0[" C4'"]=tmp; GAU0[" C4'"][0]=  -3.942; GAU0[" C4'"][1]=  -8.778; GAU0[" C4'"][2]=-352.334;
    GAU0[" C1'"]=tmp; GAU0[" C1'"][0]=  -4.738; GAU0[" C1'"][1]=  -7.057; GAU0[" C1'"][2]=-350.953;
    ideal_rna["  G  A  U0"]=GAU0;
    map<string,vector<float> >().swap(GAU0);

    map<string,vector<float> >GAU1;
    GAU1[" P  "]=tmp; GAU1[" P  "][0]=  -4.740; GAU1[" P  "][1]=  -7.164; GAU1[" P  "][2]=-355.759;
    GAU1[" C4'"]=tmp; GAU1[" C4'"][0]=  -8.060; GAU1[" C4'"][1]=  -5.257; GAU1[" C4'"][2]=-354.882;
    GAU1[" C1'"]=tmp; GAU1[" C1'"][0]=  -7.799; GAU1[" C1'"][1]=  -3.379; GAU1[" C1'"][2]=-353.501;
    ideal_rna["  G  A  U1"]=GAU1;
    map<string,vector<float> >().swap(GAU1);

    map<string,vector<float> >GAU2;
    GAU2[" P  "]=tmp; GAU2[" P  "][0]=  -7.859; GAU2[" P  "][1]=  -3.468; GAU2[" P  "][2]=-358.307;
    GAU2[" C4'"]=tmp; GAU2[" C4'"][0]=  -9.623; GAU2[" C4'"][1]=  -0.070; GAU2[" C4'"][2]=-357.430;
    GAU2[" C1'"]=tmp; GAU2[" C1'"][0]=  -8.389; GAU2[" C1'"][1]=   1.370; GAU2[" C1'"][2]=-356.049;
    ideal_rna["  G  A  U2"]=GAU2;
    map<string,vector<float> >().swap(GAU2);

    map<string,vector<float> >GCA0;
    GCA0[" P  "]=tmp; GCA0[" P  "][0]=  -8.487; GCA0[" P  "][1]=   1.328; GCA0[" P  "][2]=-360.855;
    GCA0[" C4'"]=tmp; GCA0[" C4'"][0]=  -8.135; GCA0[" C4'"][1]=   5.140; GCA0[" C4'"][2]=-359.978;
    GCA0[" C1'"]=tmp; GCA0[" C1'"][0]=  -6.319; GCA0[" C1'"][1]=   5.685; GCA0[" C1'"][2]=-358.597;
    ideal_rna["  G  C  A0"]=GCA0;
    map<string,vector<float> >().swap(GCA0);

    map<string,vector<float> >GCA1;
    GCA1[" P  "]=tmp; GCA1[" P  "][0]=  -6.424; GCA1[" P  "][1]=   5.702; GCA1[" P  "][2]=-363.403;
    GCA1[" C4'"]=tmp; GCA1[" C4'"][0]=  -4.069; GCA1[" C4'"][1]=   8.720; GCA1[" C4'"][2]=-362.526;
    GCA1[" C1'"]=tmp; GCA1[" C1'"][0]=  -2.247; GCA1[" C1'"][1]=   8.198; GCA1[" C1'"][2]=-361.145;
    ideal_rna["  G  C  A1"]=GCA1;
    map<string,vector<float> >().swap(GCA1);

    map<string,vector<float> >GCA2;
    GCA2[" P  "]=tmp; GCA2[" P  "][0]=  -2.326; GCA2[" P  "][1]=   8.269; GCA2[" P  "][2]=-365.951;
    GCA2[" C4'"]=tmp; GCA2[" C4'"][0]=   1.287; GCA2[" C4'"][1]=   9.537; GCA2[" C4'"][2]=-365.074;
    GCA2[" C1'"]=tmp; GCA2[" C1'"][0]=   2.538; GCA2[" C1'"][1]=   8.112; GCA2[" C1'"][2]=-363.693;
    ideal_rna["  G  C  A2"]=GCA2;
    map<string,vector<float> >().swap(GCA2);

    map<string,vector<float> >GCC0;
    GCC0[" P  "]=tmp; GCC0[" P  "][0]=   2.510; GCC0[" P  "][1]=   8.215; GCC0[" P  "][2]=-368.499;
    GCC0[" C4'"]=tmp; GCC0[" C4'"][0]=   6.235; GCC0[" C4'"][1]=   7.330; GCC0[" C4'"][2]=-367.622;
    GCC0[" C1'"]=tmp; GCC0[" C1'"][0]=   6.518; GCC0[" C1'"][1]=   5.455; GCC0[" C1'"][2]=-366.241;
    ideal_rna["  G  C  C0"]=GCC0;
    map<string,vector<float> >().swap(GCC0);

    map<string,vector<float> >GCC1;
    GCC1[" P  "]=tmp; GCC1[" P  "][0]=   6.550; GCC1[" P  "][1]=   5.557; GCC1[" P  "][2]=-371.047;
    GCC1[" C4'"]=tmp; GCC1[" C4'"][0]=   9.207; GCC1[" C4'"][1]=   2.800; GCC1[" C4'"][2]=-370.170;
    GCC1[" C1'"]=tmp; GCC1[" C1'"][0]=   8.432; GCC1[" C1'"][1]=   1.070; GCC1[" C1'"][2]=-368.789;
    ideal_rna["  G  C  C1"]=GCC1;
    map<string,vector<float> >().swap(GCC1);

    map<string,vector<float> >GCC2;
    GCC2[" P  "]=tmp; GCC2[" P  "][0]=   8.514; GCC2[" P  "][1]=   1.137; GCC2[" P  "][2]=-373.595;
    GCC2[" C4'"]=tmp; GCC2[" C4'"][0]=   9.260; GCC2[" C4'"][1]=  -2.618; GCC2[" C4'"][2]=-372.718;
    GCC2[" C1'"]=tmp; GCC2[" C1'"][0]=   7.674; GCC2[" C1'"][1]=  -3.655; GCC2[" C1'"][2]=-371.337;
    ideal_rna["  G  C  C2"]=GCC2;
    map<string,vector<float> >().swap(GCC2);

    map<string,vector<float> >GCG0;
    GCG0[" P  "]=tmp; GCG0[" P  "][0]=   7.779; GCG0[" P  "][1]=  -3.643; GCG0[" P  "][2]=-376.143;
    GCG0[" C4'"]=tmp; GCG0[" C4'"][0]=   6.378; GCG0[" C4'"][1]=  -7.205; GCG0[" C4'"][2]=-375.266;
    GCG0[" C1'"]=tmp; GCG0[" C1'"][0]=   4.482; GCG0[" C1'"][1]=  -7.222; GCG0[" C1'"][2]=-373.885;
    ideal_rna["  G  C  G0"]=GCG0;
    map<string,vector<float> >().swap(GCG0);

    map<string,vector<float> >GCG1;
    GCG1[" P  "]=tmp; GCG1[" P  "][0]=   4.579; GCG1[" P  "][1]=  -7.268; GCG1[" P  "][2]=-378.691;
    GCG1[" C4'"]=tmp; GCG1[" C4'"][0]=   1.475; GCG1[" C4'"][1]=  -9.509; GCG1[" C4'"][2]=-377.814;
    GCG1[" C1'"]=tmp; GCG1[" C1'"][0]=  -0.129; GCG1[" C1'"][1]=  -8.499; GCG1[" C1'"][2]=-376.433;
    ideal_rna["  G  C  G1"]=GCG1;
    map<string,vector<float> >().swap(GCG1);

    map<string,vector<float> >GCG2;
    GCG2[" P  "]=tmp; GCG2[" P  "][0]=  -0.074; GCG2[" P  "][1]=  -8.590; GCG2[" P  "][2]=-381.239;
    GCG2[" C4'"]=tmp; GCG2[" C4'"][0]=  -3.896; GCG2[" C4'"][1]=  -8.799; GCG2[" C4'"][2]=-380.362;
    GCG2[" C1'"]=tmp; GCG2[" C1'"][0]=  -4.701; GCG2[" C1'"][1]=  -7.082; GCG2[" C1'"][2]=-378.981;
    ideal_rna["  G  C  G2"]=GCG2;
    map<string,vector<float> >().swap(GCG2);

    map<string,vector<float> >GCU0;
    GCU0[" P  "]=tmp; GCU0[" P  "][0]=  -4.702; GCU0[" P  "][1]=  -7.189; GCU0[" P  "][2]=-383.787;
    GCU0[" C4'"]=tmp; GCU0[" C4'"][0]=  -8.032; GCU0[" C4'"][1]=  -5.300; GCU0[" C4'"][2]=-382.910;
    GCU0[" C1'"]=tmp; GCU0[" C1'"][0]=  -7.782; GCU0[" C1'"][1]=  -3.420; GCU0[" C1'"][2]=-381.529;
    ideal_rna["  G  C  U0"]=GCU0;
    map<string,vector<float> >().swap(GCU0);

    map<string,vector<float> >GCU1;
    GCU1[" P  "]=tmp; GCU1[" P  "][0]=  -7.841; GCU1[" P  "][1]=  -3.509; GCU1[" P  "][2]=-386.335;
    GCU1[" C4'"]=tmp; GCU1[" C4'"][0]=  -9.622; GCU1[" C4'"][1]=  -0.120; GCU1[" C4'"][2]=-385.458;
    GCU1[" C1'"]=tmp; GCU1[" C1'"][0]=  -8.396; GCU1[" C1'"][1]=   1.325; GCU1[" C1'"][2]=-384.077;
    ideal_rna["  G  C  U1"]=GCU1;
    map<string,vector<float> >().swap(GCU1);

    map<string,vector<float> >GCU2;
    GCU2[" P  "]=tmp; GCU2[" P  "][0]=  -8.494; GCU2[" P  "][1]=   1.283; GCU2[" P  "][2]=-388.883;
    GCU2[" C4'"]=tmp; GCU2[" C4'"][0]=  -8.162; GCU2[" C4'"][1]=   5.097; GCU2[" C4'"][2]=-388.006;
    GCU2[" C1'"]=tmp; GCU2[" C1'"][0]=  -6.349; GCU2[" C1'"][1]=   5.652; GCU2[" C1'"][2]=-386.625;
    ideal_rna["  G  C  U2"]=GCU2;
    map<string,vector<float> >().swap(GCU2);

    map<string,vector<float> >GGA0;
    GGA0[" P  "]=tmp; GGA0[" P  "][0]=  -6.454; GGA0[" P  "][1]=   5.668; GGA0[" P  "][2]=-391.431;
    GGA0[" C4'"]=tmp; GGA0[" C4'"][0]=  -4.115; GGA0[" C4'"][1]=   8.699; GGA0[" C4'"][2]=-390.554;
    GGA0[" C1'"]=tmp; GGA0[" C1'"][0]=  -2.290; GGA0[" C1'"][1]=   8.186; GGA0[" C1'"][2]=-389.173;
    ideal_rna["  G  G  A0"]=GGA0;
    map<string,vector<float> >().swap(GGA0);

    map<string,vector<float> >GGA1;
    GGA1[" P  "]=tmp; GGA1[" P  "][0]=  -2.369; GGA1[" P  "][1]=   8.257; GGA1[" P  "][2]=-393.979;
    GGA1[" C4'"]=tmp; GGA1[" C4'"][0]=   1.237; GGA1[" C4'"][1]=   9.543; GGA1[" C4'"][2]=-393.102;
    GGA1[" C1'"]=tmp; GGA1[" C1'"][0]=   2.496; GGA1[" C1'"][1]=   8.125; GGA1[" C1'"][2]=-391.721;
    ideal_rna["  G  G  A1"]=GGA1;
    map<string,vector<float> >().swap(GGA1);

    map<string,vector<float> >GGA2;
    GGA2[" P  "]=tmp; GGA2[" P  "][0]=   2.467; GGA2[" P  "][1]=   8.228; GGA2[" P  "][2]=-396.527;
    GGA2[" C4'"]=tmp; GGA2[" C4'"][0]=   6.196; GGA2[" C4'"][1]=   7.363; GGA2[" C4'"][2]=-395.650;
    GGA2[" C1'"]=tmp; GGA2[" C1'"][0]=   6.490; GGA2[" C1'"][1]=   5.489; GGA2[" C1'"][2]=-394.269;
    ideal_rna["  G  G  A2"]=GGA2;
    map<string,vector<float> >().swap(GGA2);

    map<string,vector<float> >GGC0;
    GGC0[" P  "]=tmp; GGC0[" P  "][0]=   6.521; GGC0[" P  "][1]=   5.591; GGC0[" P  "][2]=-399.075;
    GGC0[" C4'"]=tmp; GGC0[" C4'"][0]=   9.192; GGC0[" C4'"][1]=   2.848; GGC0[" C4'"][2]=-398.198;
    GGC0[" C1'"]=tmp; GGC0[" C1'"][0]=   8.427; GGC0[" C1'"][1]=   1.113; GGC0[" C1'"][2]=-396.817;
    ideal_rna["  G  G  C0"]=GGC0;
    map<string,vector<float> >().swap(GGC0);

    map<string,vector<float> >GGC1;
    GGC1[" P  "]=tmp; GGC1[" P  "][0]=   8.508; GGC1[" P  "][1]=   1.182; GGC1[" P  "][2]=-401.623;
    GGC1[" C4'"]=tmp; GGC1[" C4'"][0]=   9.274; GGC1[" C4'"][1]=  -2.569; GGC1[" C4'"][2]=-400.746;
    GGC1[" C1'"]=tmp; GGC1[" C1'"][0]=   7.693; GGC1[" C1'"][1]=  -3.616; GGC1[" C1'"][2]=-399.365;
    ideal_rna["  G  G  C1"]=GGC1;
    map<string,vector<float> >().swap(GGC1);

    map<string,vector<float> >GGC2;
    GGC2[" P  "]=tmp; GGC2[" P  "][0]=   7.798; GGC2[" P  "][1]=  -3.602; GGC2[" P  "][2]=-404.171;
    GGC2[" C4'"]=tmp; GGC2[" C4'"][0]=   6.416; GGC2[" C4'"][1]=  -7.172; GGC2[" C4'"][2]=-403.294;
    GGC2[" C1'"]=tmp; GGC2[" C1'"][0]=   4.521; GGC2[" C1'"][1]=  -7.198; GGC2[" C1'"][2]=-401.913;
    ideal_rna["  G  G  C2"]=GGC2;
    map<string,vector<float> >().swap(GGC2);

    map<string,vector<float> >GGG0;
    GGG0[" P  "]=tmp; GGG0[" P  "][0]=   4.617; GGG0[" P  "][1]=  -7.244; GGG0[" P  "][2]=-406.719;
    GGG0[" C4'"]=tmp; GGG0[" C4'"][0]=   1.525; GGG0[" C4'"][1]=  -9.501; GGG0[" C4'"][2]=-405.842;
    GGG0[" C1'"]=tmp; GGG0[" C1'"][0]=  -0.085; GGG0[" C1'"][1]=  -8.500; GGG0[" C1'"][2]=-404.461;
    ideal_rna["  G  G  G0"]=GGG0;
    map<string,vector<float> >().swap(GGG0);

    map<string,vector<float> >GGG1;
    GGG1[" P  "]=tmp; GGG1[" P  "][0]=  -0.029; GGG1[" P  "][1]=  -8.590; GGG1[" P  "][2]=-409.267;
    GGG1[" C4'"]=tmp; GGG1[" C4'"][0]=  -3.850; GGG1[" C4'"][1]=  -8.819; GGG1[" C4'"][2]=-408.390;
    GGG1[" C1'"]=tmp; GGG1[" C1'"][0]=  -4.663; GGG1[" C1'"][1]=  -7.106; GGG1[" C1'"][2]=-407.009;
    ideal_rna["  G  G  G1"]=GGG1;
    map<string,vector<float> >().swap(GGG1);

    map<string,vector<float> >GGG2;
    GGG2[" P  "]=tmp; GGG2[" P  "][0]=  -4.665; GGG2[" P  "][1]=  -7.213; GGG2[" P  "][2]=-411.815;
    GGG2[" C4'"]=tmp; GGG2[" C4'"][0]=  -8.004; GGG2[" C4'"][1]=  -5.342; GGG2[" C4'"][2]=-410.938;
    GGG2[" C1'"]=tmp; GGG2[" C1'"][0]=  -7.764; GGG2[" C1'"][1]=  -3.461; GGG2[" C1'"][2]=-409.557;
    ideal_rna["  G  G  G2"]=GGG2;
    map<string,vector<float> >().swap(GGG2);

    map<string,vector<float> >GGU0;
    GGU0[" P  "]=tmp; GGU0[" P  "][0]=  -7.822; GGU0[" P  "][1]=  -3.550; GGU0[" P  "][2]=-414.363;
    GGU0[" C4'"]=tmp; GGU0[" C4'"][0]=  -9.621; GGU0[" C4'"][1]=  -0.171; GGU0[" C4'"][2]=-413.486;
    GGU0[" C1'"]=tmp; GGU0[" C1'"][0]=  -8.403; GGU0[" C1'"][1]=   1.282; GGU0[" C1'"][2]=-412.105;
    ideal_rna["  G  G  U0"]=GGU0;
    map<string,vector<float> >().swap(GGU0);

    map<string,vector<float> >GGU1;
    GGU1[" P  "]=tmp; GGU1[" P  "][0]=  -8.500; GGU1[" P  "][1]=   1.239; GGU1[" P  "][2]=-416.911;
    GGU1[" C4'"]=tmp; GGU1[" C4'"][0]=  -8.189; GGU1[" C4'"][1]=   5.054; GGU1[" C4'"][2]=-416.034;
    GGU1[" C1'"]=tmp; GGU1[" C1'"][0]=  -6.378; GGU1[" C1'"][1]=   5.618; GGU1[" C1'"][2]=-414.653;
    ideal_rna["  G  G  U1"]=GGU1;
    map<string,vector<float> >().swap(GGU1);

    map<string,vector<float> >GGU2;
    GGU2[" P  "]=tmp; GGU2[" P  "][0]=  -6.484; GGU2[" P  "][1]=   5.635; GGU2[" P  "][2]=-419.459;
    GGU2[" C4'"]=tmp; GGU2[" C4'"][0]=  -4.160; GGU2[" C4'"][1]=   8.677; GGU2[" C4'"][2]=-418.582;
    GGU2[" C1'"]=tmp; GGU2[" C1'"][0]=  -2.332; GGU2[" C1'"][1]=   8.174; GGU2[" C1'"][2]=-417.201;
    ideal_rna["  G  G  U2"]=GGU2;
    map<string,vector<float> >().swap(GGU2);

    map<string,vector<float> >GUA0;
    GUA0[" P  "]=tmp; GUA0[" P  "][0]=  -2.412; GUA0[" P  "][1]=   8.244; GUA0[" P  "][2]=-422.007;
    GUA0[" C4'"]=tmp; GUA0[" C4'"][0]=   1.187; GUA0[" C4'"][1]=   9.550; GUA0[" C4'"][2]=-421.130;
    GUA0[" C1'"]=tmp; GUA0[" C1'"][0]=   2.453; GUA0[" C1'"][1]=   8.138; GUA0[" C1'"][2]=-419.749;
    ideal_rna["  G  U  A0"]=GUA0;
    map<string,vector<float> >().swap(GUA0);

    map<string,vector<float> >GUA1;
    GUA1[" P  "]=tmp; GUA1[" P  "][0]=   2.424; GUA1[" P  "][1]=   8.241; GUA1[" P  "][2]=-424.555;
    GUA1[" C4'"]=tmp; GUA1[" C4'"][0]=   6.158; GUA1[" C4'"][1]=   7.395; GUA1[" C4'"][2]=-423.678;
    GUA1[" C1'"]=tmp; GUA1[" C1'"][0]=   6.461; GUA1[" C1'"][1]=   5.523; GUA1[" C1'"][2]=-422.297;
    ideal_rna["  G  U  A1"]=GUA1;
    map<string,vector<float> >().swap(GUA1);

    map<string,vector<float> >GUA2;
    GUA2[" P  "]=tmp; GUA2[" P  "][0]=   6.492; GUA2[" P  "][1]=   5.625; GUA2[" P  "][2]=-427.103;
    GUA2[" C4'"]=tmp; GUA2[" C4'"][0]=   9.177; GUA2[" C4'"][1]=   2.896; GUA2[" C4'"][2]=-426.226;
    GUA2[" C1'"]=tmp; GUA2[" C1'"][0]=   8.421; GUA2[" C1'"][1]=   1.157; GUA2[" C1'"][2]=-424.845;
    ideal_rna["  G  U  A2"]=GUA2;
    map<string,vector<float> >().swap(GUA2);

    map<string,vector<float> >GUC0;
    GUC0[" P  "]=tmp; GUC0[" P  "][0]=   8.502; GUC0[" P  "][1]=   1.227; GUC0[" P  "][2]=-429.651;
    GUC0[" C4'"]=tmp; GUC0[" C4'"][0]=   9.287; GUC0[" C4'"][1]=  -2.520; GUC0[" C4'"][2]=-428.774;
    GUC0[" C1'"]=tmp; GUC0[" C1'"][0]=   7.711; GUC0[" C1'"][1]=  -3.575; GUC0[" C1'"][2]=-427.393;
    ideal_rna["  G  U  C0"]=GUC0;
    map<string,vector<float> >().swap(GUC0);

    map<string,vector<float> >GUC1;
    GUC1[" P  "]=tmp; GUC1[" P  "][0]=   7.817; GUC1[" P  "][1]=  -3.561; GUC1[" P  "][2]=-432.199;
    GUC1[" C4'"]=tmp; GUC1[" C4'"][0]=   6.454; GUC1[" C4'"][1]=  -7.138; GUC1[" C4'"][2]=-431.322;
    GUC1[" C1'"]=tmp; GUC1[" C1'"][0]=   4.558; GUC1[" C1'"][1]=  -7.175; GUC1[" C1'"][2]=-429.941;
    ideal_rna["  G  U  C1"]=GUC1;
    map<string,vector<float> >().swap(GUC1);

    map<string,vector<float> >GUC2;
    GUC2[" P  "]=tmp; GUC2[" P  "][0]=   4.654; GUC2[" P  "][1]=  -7.220; GUC2[" P  "][2]=-434.747;
    GUC2[" C4'"]=tmp; GUC2[" C4'"][0]=   1.574; GUC2[" C4'"][1]=  -9.493; GUC2[" C4'"][2]=-433.870;
    GUC2[" C1'"]=tmp; GUC2[" C1'"][0]=  -0.040; GUC2[" C1'"][1]=  -8.500; GUC2[" C1'"][2]=-432.489;
    ideal_rna["  G  U  C2"]=GUC2;
    map<string,vector<float> >().swap(GUC2);

    map<string,vector<float> >GUG0;
    GUG0[" P  "]=tmp; GUG0[" P  "][0]=   0.016; GUG0[" P  "][1]=  -8.590; GUG0[" P  "][2]=-437.295;
    GUG0[" C4'"]=tmp; GUG0[" C4'"][0]=  -3.804; GUG0[" C4'"][1]=  -8.839; GUG0[" C4'"][2]=-436.418;
    GUG0[" C1'"]=tmp; GUG0[" C1'"][0]=  -4.626; GUG0[" C1'"][1]=  -7.131; GUG0[" C1'"][2]=-435.037;
    ideal_rna["  G  U  G0"]=GUG0;
    map<string,vector<float> >().swap(GUG0);

    map<string,vector<float> >GUG1;
    GUG1[" P  "]=tmp; GUG1[" P  "][0]=  -4.627; GUG1[" P  "][1]=  -7.237; GUG1[" P  "][2]=-439.843;
    GUG1[" C4'"]=tmp; GUG1[" C4'"][0]=  -7.976; GUG1[" C4'"][1]=  -5.383; GUG1[" C4'"][2]=-438.966;
    GUG1[" C1'"]=tmp; GUG1[" C1'"][0]=  -7.745; GUG1[" C1'"][1]=  -3.501; GUG1[" C1'"][2]=-437.585;
    ideal_rna["  G  U  G1"]=GUG1;
    map<string,vector<float> >().swap(GUG1);

    map<string,vector<float> >GUG2;
    GUG2[" P  "]=tmp; GUG2[" P  "][0]=  -7.804; GUG2[" P  "][1]=  -3.591; GUG2[" P  "][2]=-442.391;
    GUG2[" C4'"]=tmp; GUG2[" C4'"][0]=  -9.620; GUG2[" C4'"][1]=  -0.221; GUG2[" C4'"][2]=-441.514;
    GUG2[" C1'"]=tmp; GUG2[" C1'"][0]=  -8.409; GUG2[" C1'"][1]=   1.238; GUG2[" C1'"][2]=-440.133;
    ideal_rna["  G  U  G2"]=GUG2;
    map<string,vector<float> >().swap(GUG2);

    map<string,vector<float> >GUU0;
    GUU0[" P  "]=tmp; GUU0[" P  "][0]=  -8.507; GUU0[" P  "][1]=   1.194; GUU0[" P  "][2]=-444.939;
    GUU0[" C4'"]=tmp; GUU0[" C4'"][0]=  -8.215; GUU0[" C4'"][1]=   5.011; GUU0[" C4'"][2]=-444.062;
    GUU0[" C1'"]=tmp; GUU0[" C1'"][0]=  -6.408; GUU0[" C1'"][1]=   5.585; GUU0[" C1'"][2]=-442.681;
    ideal_rna["  G  U  U0"]=GUU0;
    map<string,vector<float> >().swap(GUU0);

    map<string,vector<float> >GUU1;
    GUU1[" P  "]=tmp; GUU1[" P  "][0]=  -6.513; GUU1[" P  "][1]=   5.601; GUU1[" P  "][2]=-447.487;
    GUU1[" C4'"]=tmp; GUU1[" C4'"][0]=  -4.206; GUU1[" C4'"][1]=   8.655; GUU1[" C4'"][2]=-446.610;
    GUU1[" C1'"]=tmp; GUU1[" C1'"][0]=  -2.375; GUU1[" C1'"][1]=   8.161; GUU1[" C1'"][2]=-445.229;
    ideal_rna["  G  U  U1"]=GUU1;
    map<string,vector<float> >().swap(GUU1);

    map<string,vector<float> >GUU2;
    GUU2[" P  "]=tmp; GUU2[" P  "][0]=  -2.455; GUU2[" P  "][1]=   8.232; GUU2[" P  "][2]=-450.035;
    GUU2[" C4'"]=tmp; GUU2[" C4'"][0]=   1.137; GUU2[" C4'"][1]=   9.556; GUU2[" C4'"][2]=-449.158;
    GUU2[" C1'"]=tmp; GUU2[" C1'"][0]=   2.410; GUU2[" C1'"][1]=   8.151; GUU2[" C1'"][2]=-447.777;
    ideal_rna["  G  U  U2"]=GUU2;
    map<string,vector<float> >().swap(GUU2);

    map<string,vector<float> >UAA0;
    UAA0[" P  "]=tmp; UAA0[" P  "][0]=   2.381; UAA0[" P  "][1]=   8.253; UAA0[" P  "][2]=-452.583;
    UAA0[" C4'"]=tmp; UAA0[" C4'"][0]=   6.119; UAA0[" C4'"][1]=   7.427; UAA0[" C4'"][2]=-451.706;
    UAA0[" C1'"]=tmp; UAA0[" C1'"][0]=   6.432; UAA0[" C1'"][1]=   5.557; UAA0[" C1'"][2]=-450.325;
    ideal_rna["  U  A  A0"]=UAA0;
    map<string,vector<float> >().swap(UAA0);

    map<string,vector<float> >UAA1;
    UAA1[" P  "]=tmp; UAA1[" P  "][0]=   6.462; UAA1[" P  "][1]=   5.659; UAA1[" P  "][2]=-455.131;
    UAA1[" C4'"]=tmp; UAA1[" C4'"][0]=   9.162; UAA1[" C4'"][1]=   2.944; UAA1[" C4'"][2]=-454.254;
    UAA1[" C1'"]=tmp; UAA1[" C1'"][0]=   8.415; UAA1[" C1'"][1]=   1.201; UAA1[" C1'"][2]=-452.873;
    ideal_rna["  U  A  A1"]=UAA1;
    map<string,vector<float> >().swap(UAA1);

    map<string,vector<float> >UAA2;
    UAA2[" P  "]=tmp; UAA2[" P  "][0]=   8.495; UAA2[" P  "][1]=   1.271; UAA2[" P  "][2]=-457.679;
    UAA2[" C4'"]=tmp; UAA2[" C4'"][0]=   9.300; UAA2[" C4'"][1]=  -2.472; UAA2[" C4'"][2]=-456.802;
    UAA2[" C1'"]=tmp; UAA2[" C1'"][0]=   7.730; UAA2[" C1'"][1]=  -3.535; UAA2[" C1'"][2]=-455.421;
    ideal_rna["  U  A  A2"]=UAA2;
    map<string,vector<float> >().swap(UAA2);

    map<string,vector<float> >UAC0;
    UAC0[" P  "]=tmp; UAC0[" P  "][0]=   7.836; UAC0[" P  "][1]=  -3.520; UAC0[" P  "][2]=-460.227;
    UAC0[" C4'"]=tmp; UAC0[" C4'"][0]=   6.491; UAC0[" C4'"][1]=  -7.104; UAC0[" C4'"][2]=-459.350;
    UAC0[" C1'"]=tmp; UAC0[" C1'"][0]=   4.595; UAC0[" C1'"][1]=  -7.151; UAC0[" C1'"][2]=-457.969;
    ideal_rna["  U  A  C0"]=UAC0;
    map<string,vector<float> >().swap(UAC0);

    map<string,vector<float> >UAC1;
    UAC1[" P  "]=tmp; UAC1[" P  "][0]=   4.692; UAC1[" P  "][1]=  -7.195; UAC1[" P  "][2]=-462.775;
    UAC1[" C4'"]=tmp; UAC1[" C4'"][0]=   1.624; UAC1[" C4'"][1]=  -9.485; UAC1[" C4'"][2]=-461.898;
    UAC1[" C1'"]=tmp; UAC1[" C1'"][0]=   0.004; UAC1[" C1'"][1]=  -8.500; UAC1[" C1'"][2]=-460.517;
    ideal_rna["  U  A  C1"]=UAC1;
    map<string,vector<float> >().swap(UAC1);

    map<string,vector<float> >UAC2;
    UAC2[" P  "]=tmp; UAC2[" P  "][0]=   0.061; UAC2[" P  "][1]=  -8.590; UAC2[" P  "][2]=-465.323;
    UAC2[" C4'"]=tmp; UAC2[" C4'"][0]=  -3.758; UAC2[" C4'"][1]=  -8.859; UAC2[" C4'"][2]=-464.446;
    UAC2[" C1'"]=tmp; UAC2[" C1'"][0]=  -4.588; UAC2[" C1'"][1]=  -7.155; UAC2[" C1'"][2]=-463.065;
    ideal_rna["  U  A  C2"]=UAC2;
    map<string,vector<float> >().swap(UAC2);

    map<string,vector<float> >UAG0;
    UAG0[" P  "]=tmp; UAG0[" P  "][0]=  -4.589; UAG0[" P  "][1]=  -7.262; UAG0[" P  "][2]=-467.871;
    UAG0[" C4'"]=tmp; UAG0[" C4'"][0]=  -7.948; UAG0[" C4'"][1]=  -5.425; UAG0[" C4'"][2]=-466.994;
    UAG0[" C1'"]=tmp; UAG0[" C1'"][0]=  -7.727; UAG0[" C1'"][1]=  -3.542; UAG0[" C1'"][2]=-465.613;
    ideal_rna["  U  A  G0"]=UAG0;
    map<string,vector<float> >().swap(UAG0);

    map<string,vector<float> >UAG1;
    UAG1[" P  "]=tmp; UAG1[" P  "][0]=  -7.785; UAG1[" P  "][1]=  -3.632; UAG1[" P  "][2]=-470.419;
    UAG1[" C4'"]=tmp; UAG1[" C4'"][0]=  -9.619; UAG1[" C4'"][1]=  -0.271; UAG1[" C4'"][2]=-469.542;
    UAG1[" C1'"]=tmp; UAG1[" C1'"][0]=  -8.416; UAG1[" C1'"][1]=   1.194; UAG1[" C1'"][2]=-468.161;
    ideal_rna["  U  A  G1"]=UAG1;
    map<string,vector<float> >().swap(UAG1);

    map<string,vector<float> >UAG2;
    UAG2[" P  "]=tmp; UAG2[" P  "][0]=  -8.513; UAG2[" P  "][1]=   1.150; UAG2[" P  "][2]=-472.967;
    UAG2[" C4'"]=tmp; UAG2[" C4'"][0]=  -8.241; UAG2[" C4'"][1]=   4.968; UAG2[" C4'"][2]=-472.090;
    UAG2[" C1'"]=tmp; UAG2[" C1'"][0]=  -6.437; UAG2[" C1'"][1]=   5.551; UAG2[" C1'"][2]=-470.709;
    ideal_rna["  U  A  G2"]=UAG2;
    map<string,vector<float> >().swap(UAG2);

    map<string,vector<float> >UAU0;
    UAU0[" P  "]=tmp; UAU0[" P  "][0]=  -6.542; UAU0[" P  "][1]=   5.566; UAU0[" P  "][2]=-475.515;
    UAU0[" C4'"]=tmp; UAU0[" C4'"][0]=  -4.251; UAU0[" C4'"][1]=   8.633; UAU0[" C4'"][2]=-474.638;
    UAU0[" C1'"]=tmp; UAU0[" C1'"][0]=  -2.418; UAU0[" C1'"][1]=   8.149; UAU0[" C1'"][2]=-473.257;
    ideal_rna["  U  A  U0"]=UAU0;
    map<string,vector<float> >().swap(UAU0);

    map<string,vector<float> >UAU1;
    UAU1[" P  "]=tmp; UAU1[" P  "][0]=  -2.498; UAU1[" P  "][1]=   8.219; UAU1[" P  "][2]=-478.063;
    UAU1[" C4'"]=tmp; UAU1[" C4'"][0]=   1.087; UAU1[" C4'"][1]=   9.561; UAU1[" C4'"][2]=-477.186;
    UAU1[" C1'"]=tmp; UAU1[" C1'"][0]=   2.368; UAU1[" C1'"][1]=   8.164; UAU1[" C1'"][2]=-475.805;
    ideal_rna["  U  A  U1"]=UAU1;
    map<string,vector<float> >().swap(UAU1);

    map<string,vector<float> >UAU2;
    UAU2[" P  "]=tmp; UAU2[" P  "][0]=   2.338; UAU2[" P  "][1]=   8.266; UAU2[" P  "][2]=-480.611;
    UAU2[" C4'"]=tmp; UAU2[" C4'"][0]=   6.080; UAU2[" C4'"][1]=   7.459; UAU2[" C4'"][2]=-479.734;
    UAU2[" C1'"]=tmp; UAU2[" C1'"][0]=   6.403; UAU2[" C1'"][1]=   5.591; UAU2[" C1'"][2]=-478.353;
    ideal_rna["  U  A  U2"]=UAU2;
    map<string,vector<float> >().swap(UAU2);

    map<string,vector<float> >UCA0;
    UCA0[" P  "]=tmp; UCA0[" P  "][0]=   6.433; UCA0[" P  "][1]=   5.693; UCA0[" P  "][2]=-483.159;
    UCA0[" C4'"]=tmp; UCA0[" C4'"][0]=   9.146; UCA0[" C4'"][1]=   2.992; UCA0[" C4'"][2]=-482.282;
    UCA0[" C1'"]=tmp; UCA0[" C1'"][0]=   8.408; UCA0[" C1'"][1]=   1.246; UCA0[" C1'"][2]=-480.901;
    ideal_rna["  U  C  A0"]=UCA0;
    map<string,vector<float> >().swap(UCA0);

    map<string,vector<float> >UCA1;
    UCA1[" P  "]=tmp; UCA1[" P  "][0]=   8.489; UCA1[" P  "][1]=   1.315; UCA1[" P  "][2]=-485.707;
    UCA1[" C4'"]=tmp; UCA1[" C4'"][0]=   9.313; UCA1[" C4'"][1]=  -2.423; UCA1[" C4'"][2]=-484.830;
    UCA1[" C1'"]=tmp; UCA1[" C1'"][0]=   7.749; UCA1[" C1'"][1]=  -3.494; UCA1[" C1'"][2]=-483.449;
    ideal_rna["  U  C  A1"]=UCA1;
    map<string,vector<float> >().swap(UCA1);

    map<string,vector<float> >UCA2;
    UCA2[" P  "]=tmp; UCA2[" P  "][0]=   7.854; UCA2[" P  "][1]=  -3.479; UCA2[" P  "][2]=-488.255;
    UCA2[" C4'"]=tmp; UCA2[" C4'"][0]=   6.528; UCA2[" C4'"][1]=  -7.070; UCA2[" C4'"][2]=-487.378;
    UCA2[" C1'"]=tmp; UCA2[" C1'"][0]=   4.633; UCA2[" C1'"][1]=  -7.127; UCA2[" C1'"][2]=-485.997;
    ideal_rna["  U  C  A2"]=UCA2;
    map<string,vector<float> >().swap(UCA2);

    map<string,vector<float> >UCC0;
    UCC0[" P  "]=tmp; UCC0[" P  "][0]=   4.730; UCC0[" P  "][1]=  -7.171; UCC0[" P  "][2]=-490.803;
    UCC0[" C4'"]=tmp; UCC0[" C4'"][0]=   1.674; UCC0[" C4'"][1]=  -9.476; UCC0[" C4'"][2]=-489.926;
    UCC0[" C1'"]=tmp; UCC0[" C1'"][0]=   0.048; UCC0[" C1'"][1]=  -8.500; UCC0[" C1'"][2]=-488.545;
    ideal_rna["  U  C  C0"]=UCC0;
    map<string,vector<float> >().swap(UCC0);

    map<string,vector<float> >UCC1;
    UCC1[" P  "]=tmp; UCC1[" P  "][0]=   0.106; UCC1[" P  "][1]=  -8.589; UCC1[" P  "][2]=-493.351;
    UCC1[" C4'"]=tmp; UCC1[" C4'"][0]=  -3.711; UCC1[" C4'"][1]=  -8.879; UCC1[" C4'"][2]=-492.474;
    UCC1[" C1'"]=tmp; UCC1[" C1'"][0]=  -4.551; UCC1[" C1'"][1]=  -7.179; UCC1[" C1'"][2]=-491.093;
    ideal_rna["  U  C  C1"]=UCC1;
    map<string,vector<float> >().swap(UCC1);

    map<string,vector<float> >UCC2;
    UCC2[" P  "]=tmp; UCC2[" P  "][0]=  -4.551; UCC2[" P  "][1]=  -7.285; UCC2[" P  "][2]=-495.899;
    UCC2[" C4'"]=tmp; UCC2[" C4'"][0]=  -7.920; UCC2[" C4'"][1]=  -5.467; UCC2[" C4'"][2]=-495.022;
    UCC2[" C1'"]=tmp; UCC2[" C1'"][0]=  -7.708; UCC2[" C1'"][1]=  -3.583; UCC2[" C1'"][2]=-493.641;
    ideal_rna["  U  C  C2"]=UCC2;
    map<string,vector<float> >().swap(UCC2);

    map<string,vector<float> >UCG0;
    UCG0[" P  "]=tmp; UCG0[" P  "][0]=  -7.765; UCG0[" P  "][1]=  -3.672; UCG0[" P  "][2]=-498.447;
    UCG0[" C4'"]=tmp; UCG0[" C4'"][0]=  -9.618; UCG0[" C4'"][1]=  -0.322; UCG0[" C4'"][2]=-497.570;
    UCG0[" C1'"]=tmp; UCG0[" C1'"][0]=  -8.422; UCG0[" C1'"][1]=   1.150; UCG0[" C1'"][2]=-496.189;
    ideal_rna["  U  C  G0"]=UCG0;
    map<string,vector<float> >().swap(UCG0);

    map<string,vector<float> >UCG1;
    UCG1[" P  "]=tmp; UCG1[" P  "][0]=  -8.519; UCG1[" P  "][1]=   1.105; UCG1[" P  "][2]=-500.995;
    UCG1[" C4'"]=tmp; UCG1[" C4'"][0]=  -8.267; UCG1[" C4'"][1]=   4.925; UCG1[" C4'"][2]=-500.118;
    UCG1[" C1'"]=tmp; UCG1[" C1'"][0]=  -6.466; UCG1[" C1'"][1]=   5.517; UCG1[" C1'"][2]=-498.737;
    ideal_rna["  U  C  G1"]=UCG1;
    map<string,vector<float> >().swap(UCG1);

    map<string,vector<float> >UCG2;
    UCG2[" P  "]=tmp; UCG2[" P  "][0]=  -6.572; UCG2[" P  "][1]=   5.532; UCG2[" P  "][2]=-503.543;
    UCG2[" C4'"]=tmp; UCG2[" C4'"][0]=  -4.296; UCG2[" C4'"][1]=   8.611; UCG2[" C4'"][2]=-502.666;
    UCG2[" C1'"]=tmp; UCG2[" C1'"][0]=  -2.460; UCG2[" C1'"][1]=   8.136; UCG2[" C1'"][2]=-501.285;
    ideal_rna["  U  C  G2"]=UCG2;
    map<string,vector<float> >().swap(UCG2);

    map<string,vector<float> >UCU0;
    UCU0[" P  "]=tmp; UCU0[" P  "][0]=  -2.541; UCU0[" P  "][1]=   8.205; UCU0[" P  "][2]=-506.091;
    UCU0[" C4'"]=tmp; UCU0[" C4'"][0]=   1.037; UCU0[" C4'"][1]=   9.567; UCU0[" C4'"][2]=-505.214;
    UCU0[" C1'"]=tmp; UCU0[" C1'"][0]=   2.325; UCU0[" C1'"][1]=   8.176; UCU0[" C1'"][2]=-503.833;
    ideal_rna["  U  C  U0"]=UCU0;
    map<string,vector<float> >().swap(UCU0);

    map<string,vector<float> >UCU1;
    UCU1[" P  "]=tmp; UCU1[" P  "][0]=   2.294; UCU1[" P  "][1]=   8.278; UCU1[" P  "][2]=-508.639;
    UCU1[" C4'"]=tmp; UCU1[" C4'"][0]=   6.041; UCU1[" C4'"][1]=   7.491; UCU1[" C4'"][2]=-507.762;
    UCU1[" C1'"]=tmp; UCU1[" C1'"][0]=   6.373; UCU1[" C1'"][1]=   5.625; UCU1[" C1'"][2]=-506.381;
    ideal_rna["  U  C  U1"]=UCU1;
    map<string,vector<float> >().swap(UCU1);

    map<string,vector<float> >UCU2;
    UCU2[" P  "]=tmp; UCU2[" P  "][0]=   6.403; UCU2[" P  "][1]=   5.727; UCU2[" P  "][2]=-511.187;
    UCU2[" C4'"]=tmp; UCU2[" C4'"][0]=   9.130; UCU2[" C4'"][1]=   3.040; UCU2[" C4'"][2]=-510.310;
    UCU2[" C1'"]=tmp; UCU2[" C1'"][0]=   8.402; UCU2[" C1'"][1]=   1.290; UCU2[" C1'"][2]=-508.929;
    ideal_rna["  U  C  U2"]=UCU2;
    map<string,vector<float> >().swap(UCU2);

    map<string,vector<float> >UGA0;
    UGA0[" P  "]=tmp; UGA0[" P  "][0]=   8.482; UGA0[" P  "][1]=   1.360; UGA0[" P  "][2]=-513.735;
    UGA0[" C4'"]=tmp; UGA0[" C4'"][0]=   9.326; UGA0[" C4'"][1]=  -2.374; UGA0[" C4'"][2]=-512.858;
    UGA0[" C1'"]=tmp; UGA0[" C1'"][0]=   7.767; UGA0[" C1'"][1]=  -3.454; UGA0[" C1'"][2]=-511.477;
    ideal_rna["  U  G  A0"]=UGA0;
    map<string,vector<float> >().swap(UGA0);

    map<string,vector<float> >UGA1;
    UGA1[" P  "]=tmp; UGA1[" P  "][0]=   7.872; UGA1[" P  "][1]=  -3.438; UGA1[" P  "][2]=-516.283;
    UGA1[" C4'"]=tmp; UGA1[" C4'"][0]=   6.565; UGA1[" C4'"][1]=  -7.036; UGA1[" C4'"][2]=-515.406;
    UGA1[" C1'"]=tmp; UGA1[" C1'"][0]=   4.670; UGA1[" C1'"][1]=  -7.102; UGA1[" C1'"][2]=-514.025;
    ideal_rna["  U  G  A1"]=UGA1;
    map<string,vector<float> >().swap(UGA1);

    map<string,vector<float> >UGA2;
    UGA2[" P  "]=tmp; UGA2[" P  "][0]=   4.767; UGA2[" P  "][1]=  -7.146; UGA2[" P  "][2]=-518.831;
    UGA2[" C4'"]=tmp; UGA2[" C4'"][0]=   1.723; UGA2[" C4'"][1]=  -9.467; UGA2[" C4'"][2]=-517.954;
    UGA2[" C1'"]=tmp; UGA2[" C1'"][0]=   0.093; UGA2[" C1'"][1]=  -8.499; UGA2[" C1'"][2]=-516.573;
    ideal_rna["  U  G  A2"]=UGA2;
    map<string,vector<float> >().swap(UGA2);

    map<string,vector<float> >UGC0;
    UGC0[" P  "]=tmp; UGC0[" P  "][0]=   0.151; UGC0[" P  "][1]=  -8.589; UGC0[" P  "][2]=-521.379;
    UGC0[" C4'"]=tmp; UGC0[" C4'"][0]=  -3.665; UGC0[" C4'"][1]=  -8.898; UGC0[" C4'"][2]=-520.502;
    UGC0[" C1'"]=tmp; UGC0[" C1'"][0]=  -4.514; UGC0[" C1'"][1]=  -7.203; UGC0[" C1'"][2]=-519.121;
    ideal_rna["  U  G  C0"]=UGC0;
    map<string,vector<float> >().swap(UGC0);

    map<string,vector<float> >UGC1;
    UGC1[" P  "]=tmp; UGC1[" P  "][0]=  -4.513; UGC1[" P  "][1]=  -7.309; UGC1[" P  "][2]=-523.927;
    UGC1[" C4'"]=tmp; UGC1[" C4'"][0]=  -7.891; UGC1[" C4'"][1]=  -5.508; UGC1[" C4'"][2]=-523.050;
    UGC1[" C1'"]=tmp; UGC1[" C1'"][0]=  -7.689; UGC1[" C1'"][1]=  -3.623; UGC1[" C1'"][2]=-521.669;
    ideal_rna["  U  G  C1"]=UGC1;
    map<string,vector<float> >().swap(UGC1);

    map<string,vector<float> >UGC2;
    UGC2[" P  "]=tmp; UGC2[" P  "][0]=  -7.746; UGC2[" P  "][1]=  -3.713; UGC2[" P  "][2]=-526.475;
    UGC2[" C4'"]=tmp; UGC2[" C4'"][0]=  -9.616; UGC2[" C4'"][1]=  -0.372; UGC2[" C4'"][2]=-525.598;
    UGC2[" C1'"]=tmp; UGC2[" C1'"][0]=  -8.428; UGC2[" C1'"][1]=   1.105; UGC2[" C1'"][2]=-524.217;
    ideal_rna["  U  G  C2"]=UGC2;
    map<string,vector<float> >().swap(UGC2);

    map<string,vector<float> >UGG0;
    UGG0[" P  "]=tmp; UGG0[" P  "][0]=  -8.524; UGG0[" P  "][1]=   1.060; UGG0[" P  "][2]=-529.023;
    UGG0[" C4'"]=tmp; UGG0[" C4'"][0]=  -8.293; UGG0[" C4'"][1]=   4.882; UGG0[" C4'"][2]=-528.146;
    UGG0[" C1'"]=tmp; UGG0[" C1'"][0]=  -6.495; UGG0[" C1'"][1]=   5.483; UGG0[" C1'"][2]=-526.765;
    ideal_rna["  U  G  G0"]=UGG0;
    map<string,vector<float> >().swap(UGG0);

    map<string,vector<float> >UGG1;
    UGG1[" P  "]=tmp; UGG1[" P  "][0]=  -6.600; UGG1[" P  "][1]=   5.498; UGG1[" P  "][2]=-531.571;
    UGG1[" C4'"]=tmp; UGG1[" C4'"][0]=  -4.341; UGG1[" C4'"][1]=   8.588; UGG1[" C4'"][2]=-530.694;
    UGG1[" C1'"]=tmp; UGG1[" C1'"][0]=  -2.503; UGG1[" C1'"][1]=   8.123; UGG1[" C1'"][2]=-529.313;
    ideal_rna["  U  G  G1"]=UGG1;
    map<string,vector<float> >().swap(UGG1);

    map<string,vector<float> >UGG2;
    UGG2[" P  "]=tmp; UGG2[" P  "][0]=  -2.584; UGG2[" P  "][1]=   8.192; UGG2[" P  "][2]=-534.119;
    UGG2[" C4'"]=tmp; UGG2[" C4'"][0]=   0.986; UGG2[" C4'"][1]=   9.572; UGG2[" C4'"][2]=-533.242;
    UGG2[" C1'"]=tmp; UGG2[" C1'"][0]=   2.282; UGG2[" C1'"][1]=   8.188; UGG2[" C1'"][2]=-531.861;
    ideal_rna["  U  G  G2"]=UGG2;
    map<string,vector<float> >().swap(UGG2);

    map<string,vector<float> >UGU0;
    UGU0[" P  "]=tmp; UGU0[" P  "][0]=   2.251; UGU0[" P  "][1]=   8.290; UGU0[" P  "][2]=-536.667;
    UGU0[" C4'"]=tmp; UGU0[" C4'"][0]=   6.001; UGU0[" C4'"][1]=   7.522; UGU0[" C4'"][2]=-535.790;
    UGU0[" C1'"]=tmp; UGU0[" C1'"][0]=   6.344; UGU0[" C1'"][1]=   5.657; UGU0[" C1'"][2]=-534.409;
    ideal_rna["  U  G  U0"]=UGU0;
    map<string,vector<float> >().swap(UGU0);

    map<string,vector<float> >UGU1;
    UGU1[" P  "]=tmp; UGU1[" P  "][0]=   6.373; UGU1[" P  "][1]=   5.760; UGU1[" P  "][2]=-539.215;
    UGU1[" C4'"]=tmp; UGU1[" C4'"][0]=   9.114; UGU1[" C4'"][1]=   3.088; UGU1[" C4'"][2]=-538.338;
    UGU1[" C1'"]=tmp; UGU1[" C1'"][0]=   8.395; UGU1[" C1'"][1]=   1.334; UGU1[" C1'"][2]=-536.957;
    ideal_rna["  U  G  U1"]=UGU1;
    map<string,vector<float> >().swap(UGU1);

    map<string,vector<float> >UGU2;
    UGU2[" P  "]=tmp; UGU2[" P  "][0]=   8.474; UGU2[" P  "][1]=   1.404; UGU2[" P  "][2]=-541.763;
    UGU2[" C4'"]=tmp; UGU2[" C4'"][0]=   9.338; UGU2[" C4'"][1]=  -2.325; UGU2[" C4'"][2]=-540.886;
    UGU2[" C1'"]=tmp; UGU2[" C1'"][0]=   7.785; UGU2[" C1'"][1]=  -3.413; UGU2[" C1'"][2]=-539.505;
    ideal_rna["  U  G  U2"]=UGU2;
    map<string,vector<float> >().swap(UGU2);

    map<string,vector<float> >UUA0;
    UUA0[" P  "]=tmp; UUA0[" P  "][0]=   7.890; UUA0[" P  "][1]=  -3.396; UUA0[" P  "][2]=-544.311;
    UUA0[" C4'"]=tmp; UUA0[" C4'"][0]=   6.602; UUA0[" C4'"][1]=  -7.002; UUA0[" C4'"][2]=-543.434;
    UUA0[" C1'"]=tmp; UUA0[" C1'"][0]=   4.707; UUA0[" C1'"][1]=  -7.078; UUA0[" C1'"][2]=-542.053;
    ideal_rna["  U  U  A0"]=UUA0;
    map<string,vector<float> >().swap(UUA0);

    map<string,vector<float> >UUA1;
    UUA1[" P  "]=tmp; UUA1[" P  "][0]=   4.805; UUA1[" P  "][1]=  -7.121; UUA1[" P  "][2]=-546.859;
    UUA1[" C4'"]=tmp; UUA1[" C4'"][0]=   1.773; UUA1[" C4'"][1]=  -9.458; UUA1[" C4'"][2]=-545.982;
    UUA1[" C1'"]=tmp; UUA1[" C1'"][0]=   0.137; UUA1[" C1'"][1]=  -8.499; UUA1[" C1'"][2]=-544.601;
    ideal_rna["  U  U  A1"]=UUA1;
    map<string,vector<float> >().swap(UUA1);

    map<string,vector<float> >UUA2;
    UUA2[" P  "]=tmp; UUA2[" P  "][0]=   0.196; UUA2[" P  "][1]=  -8.588; UUA2[" P  "][2]=-549.407;
    UUA2[" C4'"]=tmp; UUA2[" C4'"][0]=  -3.618; UUA2[" C4'"][1]=  -8.917; UUA2[" C4'"][2]=-548.530;
    UUA2[" C1'"]=tmp; UUA2[" C1'"][0]=  -4.476; UUA2[" C1'"][1]=  -7.226; UUA2[" C1'"][2]=-547.149;
    ideal_rna["  U  U  A2"]=UUA2;
    map<string,vector<float> >().swap(UUA2);

    map<string,vector<float> >UUC0;
    UUC0[" P  "]=tmp; UUC0[" P  "][0]=  -4.474; UUC0[" P  "][1]=  -7.333; UUC0[" P  "][2]=-551.955;
    UUC0[" C4'"]=tmp; UUC0[" C4'"][0]=  -7.862; UUC0[" C4'"][1]=  -5.549; UUC0[" C4'"][2]=-551.078;
    UUC0[" C1'"]=tmp; UUC0[" C1'"][0]=  -7.670; UUC0[" C1'"][1]=  -3.663; UUC0[" C1'"][2]=-549.697;
    ideal_rna["  U  U  C0"]=UUC0;
    map<string,vector<float> >().swap(UUC0);

    map<string,vector<float> >UUC1;
    UUC1[" P  "]=tmp; UUC1[" P  "][0]=  -7.727; UUC1[" P  "][1]=  -3.753; UUC1[" P  "][2]=-554.503;
    UUC1[" C4'"]=tmp; UUC1[" C4'"][0]=  -9.614; UUC1[" C4'"][1]=  -0.422; UUC1[" C4'"][2]=-553.626;
    UUC1[" C1'"]=tmp; UUC1[" C1'"][0]=  -8.433; UUC1[" C1'"][1]=   1.062; UUC1[" C1'"][2]=-552.245;
    ideal_rna["  U  U  C1"]=UUC1;
    map<string,vector<float> >().swap(UUC1);

    map<string,vector<float> >UUC2;
    UUC2[" P  "]=tmp; UUC2[" P  "][0]=  -8.530; UUC2[" P  "][1]=   1.016; UUC2[" P  "][2]=-557.051;
    UUC2[" C4'"]=tmp; UUC2[" C4'"][0]=  -8.318; UUC2[" C4'"][1]=   4.838; UUC2[" C4'"][2]=-556.174;
    UUC2[" C1'"]=tmp; UUC2[" C1'"][0]=  -6.524; UUC2[" C1'"][1]=   5.449; UUC2[" C1'"][2]=-554.793;
    ideal_rna["  U  U  C2"]=UUC2;
    map<string,vector<float> >().swap(UUC2);

    map<string,vector<float> >UUG0;
    UUG0[" P  "]=tmp; UUG0[" P  "][0]=  -6.629; UUG0[" P  "][1]=   5.463; UUG0[" P  "][2]=-559.599;
    UUG0[" C4'"]=tmp; UUG0[" C4'"][0]=  -4.386; UUG0[" C4'"][1]=   8.565; UUG0[" C4'"][2]=-558.722;
    UUG0[" C1'"]=tmp; UUG0[" C1'"][0]=  -2.546; UUG0[" C1'"][1]=   8.110; UUG0[" C1'"][2]=-557.341;
    ideal_rna["  U  U  G0"]=UUG0;
    map<string,vector<float> >().swap(UUG0);

    map<string,vector<float> >UUG1;
    UUG1[" P  "]=tmp; UUG1[" P  "][0]=  -2.627; UUG1[" P  "][1]=   8.178; UUG1[" P  "][2]=-562.147;
    UUG1[" C4'"]=tmp; UUG1[" C4'"][0]=   0.936; UUG1[" C4'"][1]=   9.577; UUG1[" C4'"][2]=-561.270;
    UUG1[" C1'"]=tmp; UUG1[" C1'"][0]=   2.239; UUG1[" C1'"][1]=   8.200; UUG1[" C1'"][2]=-559.889;
    ideal_rna["  U  U  G1"]=UUG1;
    map<string,vector<float> >().swap(UUG1);

    map<string,vector<float> >UUG2;
    UUG2[" P  "]=tmp; UUG2[" P  "][0]=   2.207; UUG2[" P  "][1]=   8.302; UUG2[" P  "][2]=-564.695;
    UUG2[" C4'"]=tmp; UUG2[" C4'"][0]=   5.962; UUG2[" C4'"][1]=   7.554; UUG2[" C4'"][2]=-563.818;
    UUG2[" C1'"]=tmp; UUG2[" C1'"][0]=   6.314; UUG2[" C1'"][1]=   5.690; UUG2[" C1'"][2]=-562.437;
    ideal_rna["  U  U  G2"]=UUG2;
    map<string,vector<float> >().swap(UUG2);

    map<string,vector<float> >UUU0;
    UUU0[" P  "]=tmp; UUU0[" P  "][0]=   6.342; UUU0[" P  "][1]=   5.793; UUU0[" P  "][2]=-567.243;
    UUU0[" C4'"]=tmp; UUU0[" C4'"][0]=   9.098; UUU0[" C4'"][1]=   3.135; UUU0[" C4'"][2]=-566.366;
    UUU0[" C1'"]=tmp; UUU0[" C1'"][0]=   8.388; UUU0[" C1'"][1]=   1.377; UUU0[" C1'"][2]=-564.985;
    ideal_rna["  U  U  U0"]=UUU0;
    map<string,vector<float> >().swap(UUU0);

    map<string,vector<float> >UUU1;
    UUU1[" P  "]=tmp; UUU1[" P  "][0]=   8.467; UUU1[" P  "][1]=   1.449; UUU1[" P  "][2]=-569.791;
    UUU1[" C4'"]=tmp; UUU1[" C4'"][0]=   9.350; UUU1[" C4'"][1]=  -2.276; UUU1[" C4'"][2]=-568.914;
    UUU1[" C1'"]=tmp; UUU1[" C1'"][0]=   7.802; UUU1[" C1'"][1]=  -3.372; UUU1[" C1'"][2]=-567.533;
    ideal_rna["  U  U  U1"]=UUU1;
    map<string,vector<float> >().swap(UUU1);

    map<string,vector<float> >UUU2;
    UUU2[" P  "]=tmp; UUU2[" P  "][0]=   7.908; UUU2[" P  "][1]=  -3.355; UUU2[" P  "][2]=-572.339;
    UUU2[" C4'"]=tmp; UUU2[" C4'"][0]=   6.638; UUU2[" C4'"][1]=  -6.967; UUU2[" C4'"][2]=-571.462;
    UUU2[" C1'"]=tmp; UUU2[" C1'"][0]=   4.744; UUU2[" C1'"][1]=  -7.053; UUU2[" C1'"][2]=-570.081;
    ideal_rna["  U  U  U2"]=UUU2;
    map<string,vector<float> >().swap(UUU2);


    vector<vector<float> > xyz_list1(3,tmp);
    vector<vector<float> > xyz_list2(3,tmp);
    vector<vector<float> > RotMatix;  // U
    vector<float> TranVect;  // t
    for (map<string, map<string,vector<float> > >::iterator iter = ideal_rna.begin();
        iter != ideal_rna.end(); iter++)
    {
        string key =  iter->first;
        if (key.size()==3) continue;
        int idx=(char)(key[key.size()-1])-'0';
        string nt=key.substr(idx*3,3);
    
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
