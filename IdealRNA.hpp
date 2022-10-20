/* IdealRNA.hpp - Idealized RNA nucleotide conformations */
#include <string>
#include <map>
#include <vector>

using namespace std;

map<string, map<string,vector<double> > >parse_ideal_rna()
{
    vector<double> tmp(3,0);
    map<string, map<string,vector<double> > >ideal_rna;

    map<string,vector<double> >P;
    P[" O3'"]=tmp; P[" O3'"][0]=   7.397; P[" O3'"][1]=   5.908; P[" O3'"][2]=  -5.390;
    P[" P  "]=tmp; P[" P  "][0]=   6.913; P[" P  "][1]=   5.099; P[" P  "][2]=  -6.683;
    P[" OP1"]=tmp; P[" OP1"][0]=   7.496; P[" OP1"][1]=   5.711; P[" OP1"][2]=  -7.898;
    P[" OP2"]=tmp; P[" OP2"][0]=   5.439; P[" OP2"][1]=   4.971; P[" OP2"][2]=  -6.666;
    P[" O5'"]=tmp; P[" O5'"][0]=   7.579; P[" O5'"][1]=   3.667; P[" O5'"][2]=  -6.429;
    ideal_rna["P"]=P;
    map<string,vector<double> >().swap(P);

    map<string,vector<double> >A;
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
    map<string,vector<double> >().swap(A);

    map<string,vector<double> >C;
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
    map<string,vector<double> >().swap(C);

    map<string,vector<double> >G;
    G[" P  "]=tmp; G[" P  "][0]=  -1.758; G[" P  "][1]=  -8.408; G[" P  "][2]=   1.587;
    G[" OP1"]=tmp; G[" OP1"][0]=  -2.072; G[" OP1"][1]=  -9.193; G[" OP1"][2]=   2.802;
    G[" OP2"]=tmp; G[" OP2"][0]=  -2.256; G[" OP2"][1]=  -7.014; G[" OP2"][2]=   1.570;
    G[" O5'"]=tmp; G[" O5'"][0]=  -0.179; G[" O5'"][1]=  -8.418; G[" O5'"][2]=   1.333;
    G[" C5'"]=tmp; G[" C5'"][0]=   0.473; G[" C5'"][1]=  -9.668; G[" C5'"][2]=   1.039;
    G[" C4'"]=tmp; G[" C4'"][0]=   1.932; G[" C4'"][1]=  -9.427; G[" C4'"][2]=   0.710;
    G[" O4'"]=tmp; G[" O4'"][0]=   2.020; G[" O4'"][1]=  -8.859; G[" O4'"][2]=  -0.629;
    G[" C3'"]=tmp; G[" C3'"][0]=   2.653; G[" C3'"][1]=  -8.408; G[" C3'"][2]=   1.592;
    G[" O3'"]=tmp; G[" O3'"][0]=   3.034; G[" O3'"][1]=  -8.969; G[" O3'"][2]=   2.842;
    G[" C2'"]=tmp; G[" C2'"][0]=   3.815; G[" C2'"][1]=  -8.065; G[" C2'"][2]=   0.660;
    G[" O2'"]=tmp; G[" O2'"][0]=   4.887; G[" O2'"][1]=  -8.973; G[" O2'"][2]=   0.520;
    G[" C1'"]=tmp; G[" C1'"][0]=   3.084; G[" C1'"][1]=  -7.921; G[" C1'"][2]=  -0.671;
    G[" N9 "]=tmp; G[" N9 "][0]=   2.522; G[" N9 "][1]=  -6.560; G[" N9 "][2]=  -0.901;
    G[" C8 "]=tmp; G[" C8 "][0]=   1.251; G[" C8 "][1]=  -6.088; G[" C8 "][2]=  -0.649;
    G[" N7 "]=tmp; G[" N7 "][0]=   1.079; G[" N7 "][1]=  -4.828; G[" N7 "][2]=  -0.969;
    G[" C5 "]=tmp; G[" C5 "][0]=   2.320; G[" C5 "][1]=  -4.440; G[" C5 "][2]=  -1.469;
    G[" C6 "]=tmp; G[" C6 "][0]=   2.748; G[" C6 "][1]=  -3.184; G[" C6 "][2]=  -1.974;
    G[" O6 "]=tmp; G[" O6 "][0]=   2.110; G[" O6 "][1]=  -2.141; G[" O6 "][2]=  -2.086;
    G[" N1 "]=tmp; G[" N1 "][0]=   4.091; G[" N1 "][1]=  -3.226; G[" N1 "][2]=  -2.376;
    G[" C2 "]=tmp; G[" C2 "][0]=   4.910; G[" C2 "][1]=  -4.335; G[" C2 "][2]=  -2.300;
    G[" N2 "]=tmp; G[" N2 "][0]=   6.162; G[" N2 "][1]=  -4.170; G[" N2 "][2]=  -2.736;
    G[" N3 "]=tmp; G[" N3 "][0]=   4.507; G[" N3 "][1]=  -5.512; G[" N3 "][2]=  -1.826;
    G[" C4 "]=tmp; G[" C4 "][0]=   3.207; G[" C4 "][1]=  -5.489; G[" C4 "][2]=  -1.431;
    ideal_rna["  G"]=G;
    map<string,vector<double> >().swap(G);

    map<string,vector<double> >U;
    U[" P  "]=tmp; U[" P  "][0]=   3.063; U[" P  "][1]=  -8.025; U[" P  "][2]=   4.135;
    U[" OP1"]=tmp; U[" OP1"][0]=   3.223; U[" OP1"][1]=  -8.856; U[" OP1"][2]=   5.350;
    U[" OP2"]=tmp; U[" OP2"][0]=   1.891; U[" OP2"][1]=  -7.121; U[" OP2"][2]=   4.118;
    U[" O5'"]=tmp; U[" O5'"][0]=   4.397; U[" O5'"][1]=  -7.181; U[" O5'"][2]=   3.881;
    U[" C5'"]=tmp; U[" C5'"][0]=   5.621; U[" C5'"][1]=  -7.881; U[" C5'"][2]=   3.587;
    U[" C4'"]=tmp; U[" C4'"][0]=   6.719; U[" C4'"][1]=  -6.889; U[" C4'"][2]=   3.258;
    U[" O4'"]=tmp; U[" O4'"][0]=   6.486; U[" O4'"][1]=  -6.364; U[" O4'"][2]=   1.919;
    U[" C3'"]=tmp; U[" C3'"][0]=   6.776; U[" C3'"][1]=  -5.642; U[" C3'"][2]=   4.140;
    U[" O3'"]=tmp; U[" O3'"][0]=   7.397; U[" O3'"][1]=  -5.908; U[" O3'"][2]=   5.390;
    U[" C2'"]=tmp; U[" C2'"][0]=   7.567; U[" C2'"][1]=  -4.725; U[" C2'"][2]=   3.208;
    U[" O2'"]=tmp; U[" O2'"][0]=   8.962; U[" O2'"][1]=  -4.910; U[" O2'"][2]=   3.068;
    U[" C1'"]=tmp; U[" C1'"][0]=   6.874; U[" C1'"][1]=  -4.999; U[" C1'"][2]=   1.877;
    U[" N1 "]=tmp; U[" N1 "][0]=   5.666; U[" N1 "][1]=  -4.158; U[" N1 "][2]=   1.647;
    U[" C2 "]=tmp; U[" C2 "][0]=   5.867; U[" C2 "][1]=  -2.911; U[" C2 "][2]=   1.107;
    U[" O2 "]=tmp; U[" O2 "][0]=   6.971; U[" O2 "][1]=  -2.482; U[" O2 "][2]=   0.819;
    U[" N3 "]=tmp; U[" N3 "][0]=   4.725; U[" N3 "][1]=  -2.160; U[" N3 "][2]=   0.908;
    U[" C4 "]=tmp; U[" C4 "][0]=   3.431; U[" C4 "][1]=  -2.543; U[" C4 "][2]=   1.197;
    U[" O4 "]=tmp; U[" O4 "][0]=   2.487; U[" O4 "][1]=  -1.783; U[" O4 "][2]=   0.974;
    U[" C5 "]=tmp; U[" C5 "][0]=   3.322; U[" C5 "][1]=  -3.869; U[" C5 "][2]=   1.760;
    U[" C6 "]=tmp; U[" C6 "][0]=   4.416; U[" C6 "][1]=  -4.621; U[" C6 "][2]=   1.964;
    ideal_rna["  U"]=U;
    map<string,vector<double> >().swap(U);

    map<string,vector<double> >A0;
    A0[" P  "]=tmp; A0[" P  "][0]=   3.063; A0[" P  "][1]=   8.025; A0[" P  "][2]=  -4.135;
    A0[" C4'"]=tmp; A0[" C4'"][0]=   6.719; A0[" C4'"][1]=   6.889; A0[" C4'"][2]=  -3.258;
    A0[" C1'"]=tmp; A0[" C1'"][0]=   6.874; A0[" C1'"][1]=   4.999; A0[" C1'"][2]=  -1.877;
    map<string,vector<double> >A1;
    A1[" P  "]=tmp; A1[" P  "][0]=   6.913; A1[" P  "][1]=   5.099; A1[" P  "][2]=  -6.683;
    A1[" C4'"]=tmp; A1[" C4'"][0]=   9.376; A1[" C4'"][1]=   2.167; A1[" C4'"][2]=  -5.806;
    A1[" C1'"]=tmp; A1[" C1'"][0]=   8.486; A1[" C1'"][1]=   0.494; A1[" C1'"][2]=  -4.425;
    map<string,vector<double> >A2;
    A2[" P  "]=tmp; A2[" P  "][0]=   8.572; A2[" P  "][1]=   0.556; A2[" P  "][2]=  -9.231;
    A2[" C4'"]=tmp; A2[" C4'"][0]=   9.061; A2[" C4'"][1]=  -3.241; A2[" C4'"][2]=  -8.354;
    A2[" C1'"]=tmp; A2[" C1'"][0]=   7.407; A2[" C1'"][1]=  -4.169; A2[" C1'"][2]=  -6.973;
    map<string,vector<double> >B0;
    B0[" P  "]=tmp; B0[" P  "][0]=  -6.022; B0[" P  "][1]=  -6.126; B0[" P  "][2]=  -0.961;
    B0[" C4'"]=tmp; B0[" C4'"][0]=  -3.467; B0[" C4'"][1]=  -8.977; B0[" C4'"][2]=  -1.838;
    B0[" C1'"]=tmp; B0[" C1'"][0]=  -1.685; B0[" C1'"][1]=  -8.331; B0[" C1'"][2]=  -3.219;
    map<string,vector<double> >B1;
    B1[" P  "]=tmp; B1[" P  "][0]=  -1.758; B1[" P  "][1]=  -8.408; B1[" P  "][2]=   1.587;
    B1[" C4'"]=tmp; B1[" C4'"][0]=   1.932; B1[" C4'"][1]=  -9.427; B1[" C4'"][2]=   0.710;
    B1[" C1'"]=tmp; B1[" C1'"][0]=   3.084; B1[" C1'"][1]=  -7.921; B1[" C1'"][2]=  -0.671;
    map<string,vector<double> >B2;
    B2[" P  "]=tmp; B2[" P  "][0]=   3.063; B2[" P  "][1]=  -8.025; B2[" P  "][2]=   4.135;
    B2[" C4'"]=tmp; B2[" C4'"][0]=   6.719; B2[" C4'"][1]=  -6.889; B2[" C4'"][2]=   3.258;
    B2[" C1'"]=tmp; B2[" C1'"][0]=   6.874; B2[" C1'"][1]=  -4.999; B2[" C1'"][2]=   1.877;

    ideal_rna["  A  A0"]=A0; ideal_rna["  A  A1"]=A1; ideal_rna["  A  A2"]=B1; ideal_rna["  A  A3"]=B2;
    ideal_rna["  A  C0"]=A0; ideal_rna["  A  C1"]=A1; ideal_rna["  A  C2"]=B1; ideal_rna["  A  C3"]=B2;
    ideal_rna["  A  G0"]=A0; ideal_rna["  A  G1"]=A1; ideal_rna["  A  G2"]=B1; ideal_rna["  A  G3"]=B2;
    ideal_rna["  A  U0"]=A0; ideal_rna["  A  U1"]=A1; ideal_rna["  A  U2"]=B1; ideal_rna["  A  U3"]=B2;
    ideal_rna["  C  A0"]=A0; ideal_rna["  C  A1"]=A1; ideal_rna["  C  A2"]=B1; ideal_rna["  C  A3"]=B2;
    ideal_rna["  C  C0"]=A0; ideal_rna["  C  C1"]=A1; ideal_rna["  C  C2"]=B1; ideal_rna["  C  C3"]=B2;
    ideal_rna["  C  G0"]=A0; ideal_rna["  C  G1"]=A1; ideal_rna["  C  G2"]=B1; ideal_rna["  C  G3"]=B2;
    ideal_rna["  C  U0"]=A0; ideal_rna["  C  U1"]=A1; ideal_rna["  C  U2"]=B1; ideal_rna["  C  U3"]=B2;
    ideal_rna["  G  A0"]=A0; ideal_rna["  G  A1"]=A1; ideal_rna["  G  A2"]=B1; ideal_rna["  G  A3"]=B2;
    ideal_rna["  G  C0"]=A0; ideal_rna["  G  C1"]=A1; ideal_rna["  G  C2"]=B1; ideal_rna["  G  C3"]=B2;
    ideal_rna["  G  G0"]=A0; ideal_rna["  G  G1"]=A1; ideal_rna["  G  G2"]=B1; ideal_rna["  G  G3"]=B2;
    ideal_rna["  G  U0"]=A0; ideal_rna["  G  U1"]=A1; ideal_rna["  G  U2"]=B1; ideal_rna["  G  U3"]=B2;
    ideal_rna["  U  A0"]=A0; ideal_rna["  U  A1"]=A1; ideal_rna["  U  A2"]=B1; ideal_rna["  U  A3"]=B2;
    ideal_rna["  U  C0"]=A0; ideal_rna["  U  C1"]=A1; ideal_rna["  U  C2"]=B1; ideal_rna["  U  C3"]=B2;
    ideal_rna["  U  G0"]=A0; ideal_rna["  U  G1"]=A1; ideal_rna["  U  G2"]=B1; ideal_rna["  U  G3"]=B2;
    ideal_rna["  U  U0"]=A0; ideal_rna["  U  U1"]=A1; ideal_rna["  U  U2"]=B1; ideal_rna["  U  U3"]=B2;
    ideal_rna["  A  A  A0"]=A0; ideal_rna["  A  A  A1"]=A1; ideal_rna["  A  A  A2"]=A2;  ideal_rna["  A  A  A3"]=B0; ideal_rna["  A  A  A4"]=B1; ideal_rna["  A  A  A5"]=B2;
    ideal_rna["  A  A  C0"]=A0; ideal_rna["  A  A  C1"]=A1; ideal_rna["  A  A  C2"]=A2;  ideal_rna["  A  A  C3"]=B0; ideal_rna["  A  A  C4"]=B1; ideal_rna["  A  A  C5"]=B2;
    ideal_rna["  A  A  G0"]=A0; ideal_rna["  A  A  G1"]=A1; ideal_rna["  A  A  G2"]=A2;  ideal_rna["  A  A  G3"]=B0; ideal_rna["  A  A  G4"]=B1; ideal_rna["  A  A  G5"]=B2;
    ideal_rna["  A  A  U0"]=A0; ideal_rna["  A  A  U1"]=A1; ideal_rna["  A  A  U2"]=A2;  ideal_rna["  A  A  U3"]=B0; ideal_rna["  A  A  U4"]=B1; ideal_rna["  A  A  U5"]=B2;
    ideal_rna["  A  C  A0"]=A0; ideal_rna["  A  C  A1"]=A1; ideal_rna["  A  C  A2"]=A2;  ideal_rna["  A  C  A3"]=B0; ideal_rna["  A  C  A4"]=B1; ideal_rna["  A  C  A5"]=B2;
    ideal_rna["  A  C  C0"]=A0; ideal_rna["  A  C  C1"]=A1; ideal_rna["  A  C  C2"]=A2;  ideal_rna["  A  C  C3"]=B0; ideal_rna["  A  C  C4"]=B1; ideal_rna["  A  C  C5"]=B2;
    ideal_rna["  A  C  G0"]=A0; ideal_rna["  A  C  G1"]=A1; ideal_rna["  A  C  G2"]=A2;  ideal_rna["  A  C  G3"]=B0; ideal_rna["  A  C  G4"]=B1; ideal_rna["  A  C  G5"]=B2;
    ideal_rna["  A  C  U0"]=A0; ideal_rna["  A  C  U1"]=A1; ideal_rna["  A  C  U2"]=A2;  ideal_rna["  A  C  U3"]=B0; ideal_rna["  A  C  U4"]=B1; ideal_rna["  A  C  U5"]=B2;
    ideal_rna["  A  G  A0"]=A0; ideal_rna["  A  G  A1"]=A1; ideal_rna["  A  G  A2"]=A2;  ideal_rna["  A  G  A3"]=B0; ideal_rna["  A  G  A4"]=B1; ideal_rna["  A  G  A5"]=B2;
    ideal_rna["  A  G  C0"]=A0; ideal_rna["  A  G  C1"]=A1; ideal_rna["  A  G  C2"]=A2;  ideal_rna["  A  G  C3"]=B0; ideal_rna["  A  G  C4"]=B1; ideal_rna["  A  G  C5"]=B2;
    ideal_rna["  A  G  G0"]=A0; ideal_rna["  A  G  G1"]=A1; ideal_rna["  A  G  G2"]=A2;  ideal_rna["  A  G  G3"]=B0; ideal_rna["  A  G  G4"]=B1; ideal_rna["  A  G  G5"]=B2;
    ideal_rna["  A  G  U0"]=A0; ideal_rna["  A  G  U1"]=A1; ideal_rna["  A  G  U2"]=A2;  ideal_rna["  A  G  U3"]=B0; ideal_rna["  A  G  U4"]=B1; ideal_rna["  A  G  U5"]=B2;
    ideal_rna["  A  U  A0"]=A0; ideal_rna["  A  U  A1"]=A1; ideal_rna["  A  U  A2"]=A2;  ideal_rna["  A  U  A3"]=B0; ideal_rna["  A  U  A4"]=B1; ideal_rna["  A  U  A5"]=B2;
    ideal_rna["  A  U  C0"]=A0; ideal_rna["  A  U  C1"]=A1; ideal_rna["  A  U  C2"]=A2;  ideal_rna["  A  U  C3"]=B0; ideal_rna["  A  U  C4"]=B1; ideal_rna["  A  U  C5"]=B2;
    ideal_rna["  A  U  G0"]=A0; ideal_rna["  A  U  G1"]=A1; ideal_rna["  A  U  G2"]=A2;  ideal_rna["  A  U  G3"]=B0; ideal_rna["  A  U  G4"]=B1; ideal_rna["  A  U  G5"]=B2;
    ideal_rna["  A  U  U0"]=A0; ideal_rna["  A  U  U1"]=A1; ideal_rna["  A  U  U2"]=A2;  ideal_rna["  A  U  U3"]=B0; ideal_rna["  A  U  U4"]=B1; ideal_rna["  A  U  U5"]=B2;
    ideal_rna["  C  A  A0"]=A0; ideal_rna["  C  A  A1"]=A1; ideal_rna["  C  A  A2"]=A2;  ideal_rna["  C  A  A3"]=B0; ideal_rna["  C  A  A4"]=B1; ideal_rna["  C  A  A5"]=B2;
    ideal_rna["  C  A  C0"]=A0; ideal_rna["  C  A  C1"]=A1; ideal_rna["  C  A  C2"]=A2;  ideal_rna["  C  A  C3"]=B0; ideal_rna["  C  A  C4"]=B1; ideal_rna["  C  A  C5"]=B2;
    ideal_rna["  C  A  G0"]=A0; ideal_rna["  C  A  G1"]=A1; ideal_rna["  C  A  G2"]=A2;  ideal_rna["  C  A  G3"]=B0; ideal_rna["  C  A  G4"]=B1; ideal_rna["  C  A  G5"]=B2;
    ideal_rna["  C  A  U0"]=A0; ideal_rna["  C  A  U1"]=A1; ideal_rna["  C  A  U2"]=A2;  ideal_rna["  C  A  U3"]=B0; ideal_rna["  C  A  U4"]=B1; ideal_rna["  C  A  U5"]=B2;
    ideal_rna["  C  C  A0"]=A0; ideal_rna["  C  C  A1"]=A1; ideal_rna["  C  C  A2"]=A2;  ideal_rna["  C  C  A3"]=B0; ideal_rna["  C  C  A4"]=B1; ideal_rna["  C  C  A5"]=B2;
    ideal_rna["  C  C  C0"]=A0; ideal_rna["  C  C  C1"]=A1; ideal_rna["  C  C  C2"]=A2;  ideal_rna["  C  C  C3"]=B0; ideal_rna["  C  C  C4"]=B1; ideal_rna["  C  C  C5"]=B2;
    ideal_rna["  C  C  G0"]=A0; ideal_rna["  C  C  G1"]=A1; ideal_rna["  C  C  G2"]=A2;  ideal_rna["  C  C  G3"]=B0; ideal_rna["  C  C  G4"]=B1; ideal_rna["  C  C  G5"]=B2;
    ideal_rna["  C  C  U0"]=A0; ideal_rna["  C  C  U1"]=A1; ideal_rna["  C  C  U2"]=A2;  ideal_rna["  C  C  U3"]=B0; ideal_rna["  C  C  U4"]=B1; ideal_rna["  C  C  U5"]=B2;
    ideal_rna["  C  G  A0"]=A0; ideal_rna["  C  G  A1"]=A1; ideal_rna["  C  G  A2"]=A2;  ideal_rna["  C  G  A3"]=B0; ideal_rna["  C  G  A4"]=B1; ideal_rna["  C  G  A5"]=B2;
    ideal_rna["  C  G  C0"]=A0; ideal_rna["  C  G  C1"]=A1; ideal_rna["  C  G  C2"]=A2;  ideal_rna["  C  G  C3"]=B0; ideal_rna["  C  G  C4"]=B1; ideal_rna["  C  G  C5"]=B2;
    ideal_rna["  C  G  G0"]=A0; ideal_rna["  C  G  G1"]=A1; ideal_rna["  C  G  G2"]=A2;  ideal_rna["  C  G  G3"]=B0; ideal_rna["  C  G  G4"]=B1; ideal_rna["  C  G  G5"]=B2;
    ideal_rna["  C  G  U0"]=A0; ideal_rna["  C  G  U1"]=A1; ideal_rna["  C  G  U2"]=A2;  ideal_rna["  C  G  U3"]=B0; ideal_rna["  C  G  U4"]=B1; ideal_rna["  C  G  U5"]=B2;
    ideal_rna["  C  U  A0"]=A0; ideal_rna["  C  U  A1"]=A1; ideal_rna["  C  U  A2"]=A2;  ideal_rna["  C  U  A3"]=B0; ideal_rna["  C  U  A4"]=B1; ideal_rna["  C  U  A5"]=B2;
    ideal_rna["  C  U  C0"]=A0; ideal_rna["  C  U  C1"]=A1; ideal_rna["  C  U  C2"]=A2;  ideal_rna["  C  U  C3"]=B0; ideal_rna["  C  U  C4"]=B1; ideal_rna["  C  U  C5"]=B2;
    ideal_rna["  C  U  G0"]=A0; ideal_rna["  C  U  G1"]=A1; ideal_rna["  C  U  G2"]=A2;  ideal_rna["  C  U  G3"]=B0; ideal_rna["  C  U  G4"]=B1; ideal_rna["  C  U  G5"]=B2;
    ideal_rna["  C  U  U0"]=A0; ideal_rna["  C  U  U1"]=A1; ideal_rna["  C  U  U2"]=A2;  ideal_rna["  C  U  U3"]=B0; ideal_rna["  C  U  U4"]=B1; ideal_rna["  C  U  U5"]=B2;
    ideal_rna["  G  A  A0"]=A0; ideal_rna["  G  A  A1"]=A1; ideal_rna["  G  A  A2"]=A2;  ideal_rna["  G  A  A3"]=B0; ideal_rna["  G  A  A4"]=B1; ideal_rna["  G  A  A5"]=B2;
    ideal_rna["  G  A  C0"]=A0; ideal_rna["  G  A  C1"]=A1; ideal_rna["  G  A  C2"]=A2;  ideal_rna["  G  A  C3"]=B0; ideal_rna["  G  A  C4"]=B1; ideal_rna["  G  A  C5"]=B2;
    ideal_rna["  G  A  G0"]=A0; ideal_rna["  G  A  G1"]=A1; ideal_rna["  G  A  G2"]=A2;  ideal_rna["  G  A  G3"]=B0; ideal_rna["  G  A  G4"]=B1; ideal_rna["  G  A  G5"]=B2;
    ideal_rna["  G  A  U0"]=A0; ideal_rna["  G  A  U1"]=A1; ideal_rna["  G  A  U2"]=A2;  ideal_rna["  G  A  U3"]=B0; ideal_rna["  G  A  U4"]=B1; ideal_rna["  G  A  U5"]=B2;
    ideal_rna["  G  C  A0"]=A0; ideal_rna["  G  C  A1"]=A1; ideal_rna["  G  C  A2"]=A2;  ideal_rna["  G  C  A3"]=B0; ideal_rna["  G  C  A4"]=B1; ideal_rna["  G  C  A5"]=B2;
    ideal_rna["  G  C  C0"]=A0; ideal_rna["  G  C  C1"]=A1; ideal_rna["  G  C  C2"]=A2;  ideal_rna["  G  C  C3"]=B0; ideal_rna["  G  C  C4"]=B1; ideal_rna["  G  C  C5"]=B2;
    ideal_rna["  G  C  G0"]=A0; ideal_rna["  G  C  G1"]=A1; ideal_rna["  G  C  G2"]=A2;  ideal_rna["  G  C  G3"]=B0; ideal_rna["  G  C  G4"]=B1; ideal_rna["  G  C  G5"]=B2;
    ideal_rna["  G  C  U0"]=A0; ideal_rna["  G  C  U1"]=A1; ideal_rna["  G  C  U2"]=A2;  ideal_rna["  G  C  U3"]=B0; ideal_rna["  G  C  U4"]=B1; ideal_rna["  G  C  U5"]=B2;
    ideal_rna["  G  G  A0"]=A0; ideal_rna["  G  G  A1"]=A1; ideal_rna["  G  G  A2"]=A2;  ideal_rna["  G  G  A3"]=B0; ideal_rna["  G  G  A4"]=B1; ideal_rna["  G  G  A5"]=B2;
    ideal_rna["  G  G  C0"]=A0; ideal_rna["  G  G  C1"]=A1; ideal_rna["  G  G  C2"]=A2;  ideal_rna["  G  G  C3"]=B0; ideal_rna["  G  G  C4"]=B1; ideal_rna["  G  G  C5"]=B2;
    ideal_rna["  G  G  G0"]=A0; ideal_rna["  G  G  G1"]=A1; ideal_rna["  G  G  G2"]=A2;  ideal_rna["  G  G  G3"]=B0; ideal_rna["  G  G  G4"]=B1; ideal_rna["  G  G  G5"]=B2;
    ideal_rna["  G  G  U0"]=A0; ideal_rna["  G  G  U1"]=A1; ideal_rna["  G  G  U2"]=A2;  ideal_rna["  G  G  U3"]=B0; ideal_rna["  G  G  U4"]=B1; ideal_rna["  G  G  U5"]=B2;
    ideal_rna["  G  U  A0"]=A0; ideal_rna["  G  U  A1"]=A1; ideal_rna["  G  U  A2"]=A2;  ideal_rna["  G  U  A3"]=B0; ideal_rna["  G  U  A4"]=B1; ideal_rna["  G  U  A5"]=B2;
    ideal_rna["  G  U  C0"]=A0; ideal_rna["  G  U  C1"]=A1; ideal_rna["  G  U  C2"]=A2;  ideal_rna["  G  U  C3"]=B0; ideal_rna["  G  U  C4"]=B1; ideal_rna["  G  U  C5"]=B2;
    ideal_rna["  G  U  G0"]=A0; ideal_rna["  G  U  G1"]=A1; ideal_rna["  G  U  G2"]=A2;  ideal_rna["  G  U  G3"]=B0; ideal_rna["  G  U  G4"]=B1; ideal_rna["  G  U  G5"]=B2;
    ideal_rna["  G  U  U0"]=A0; ideal_rna["  G  U  U1"]=A1; ideal_rna["  G  U  U2"]=A2;  ideal_rna["  G  U  U3"]=B0; ideal_rna["  G  U  U4"]=B1; ideal_rna["  G  U  U5"]=B2;
    ideal_rna["  U  A  A0"]=A0; ideal_rna["  U  A  A1"]=A1; ideal_rna["  U  A  A2"]=A2;  ideal_rna["  U  A  A3"]=B0; ideal_rna["  U  A  A4"]=B1; ideal_rna["  U  A  A5"]=B2;
    ideal_rna["  U  A  C0"]=A0; ideal_rna["  U  A  C1"]=A1; ideal_rna["  U  A  C2"]=A2;  ideal_rna["  U  A  C3"]=B0; ideal_rna["  U  A  C4"]=B1; ideal_rna["  U  A  C5"]=B2;
    ideal_rna["  U  A  G0"]=A0; ideal_rna["  U  A  G1"]=A1; ideal_rna["  U  A  G2"]=A2;  ideal_rna["  U  A  G3"]=B0; ideal_rna["  U  A  G4"]=B1; ideal_rna["  U  A  G5"]=B2;
    ideal_rna["  U  A  U0"]=A0; ideal_rna["  U  A  U1"]=A1; ideal_rna["  U  A  U2"]=A2;  ideal_rna["  U  A  U3"]=B0; ideal_rna["  U  A  U4"]=B1; ideal_rna["  U  A  U5"]=B2;
    ideal_rna["  U  C  A0"]=A0; ideal_rna["  U  C  A1"]=A1; ideal_rna["  U  C  A2"]=A2;  ideal_rna["  U  C  A3"]=B0; ideal_rna["  U  C  A4"]=B1; ideal_rna["  U  C  A5"]=B2;
    ideal_rna["  U  C  C0"]=A0; ideal_rna["  U  C  C1"]=A1; ideal_rna["  U  C  C2"]=A2;  ideal_rna["  U  C  C3"]=B0; ideal_rna["  U  C  C4"]=B1; ideal_rna["  U  C  C5"]=B2;
    ideal_rna["  U  C  G0"]=A0; ideal_rna["  U  C  G1"]=A1; ideal_rna["  U  C  G2"]=A2;  ideal_rna["  U  C  G3"]=B0; ideal_rna["  U  C  G4"]=B1; ideal_rna["  U  C  G5"]=B2;
    ideal_rna["  U  C  U0"]=A0; ideal_rna["  U  C  U1"]=A1; ideal_rna["  U  C  U2"]=A2;  ideal_rna["  U  C  U3"]=B0; ideal_rna["  U  C  U4"]=B1; ideal_rna["  U  C  U5"]=B2;
    ideal_rna["  U  G  A0"]=A0; ideal_rna["  U  G  A1"]=A1; ideal_rna["  U  G  A2"]=A2;  ideal_rna["  U  G  A3"]=B0; ideal_rna["  U  G  A4"]=B1; ideal_rna["  U  G  A5"]=B2;
    ideal_rna["  U  G  C0"]=A0; ideal_rna["  U  G  C1"]=A1; ideal_rna["  U  G  C2"]=A2;  ideal_rna["  U  G  C3"]=B0; ideal_rna["  U  G  C4"]=B1; ideal_rna["  U  G  C5"]=B2;
    ideal_rna["  U  G  G0"]=A0; ideal_rna["  U  G  G1"]=A1; ideal_rna["  U  G  G2"]=A2;  ideal_rna["  U  G  G3"]=B0; ideal_rna["  U  G  G4"]=B1; ideal_rna["  U  G  G5"]=B2;
    ideal_rna["  U  G  U0"]=A0; ideal_rna["  U  G  U1"]=A1; ideal_rna["  U  G  U2"]=A2;  ideal_rna["  U  G  U3"]=B0; ideal_rna["  U  G  U4"]=B1; ideal_rna["  U  G  U5"]=B2;
    ideal_rna["  U  U  A0"]=A0; ideal_rna["  U  U  A1"]=A1; ideal_rna["  U  U  A2"]=A2;  ideal_rna["  U  U  A3"]=B0; ideal_rna["  U  U  A4"]=B1; ideal_rna["  U  U  A5"]=B2;
    ideal_rna["  U  U  C0"]=A0; ideal_rna["  U  U  C1"]=A1; ideal_rna["  U  U  C2"]=A2;  ideal_rna["  U  U  C3"]=B0; ideal_rna["  U  U  C4"]=B1; ideal_rna["  U  U  C5"]=B2;
    ideal_rna["  U  U  G0"]=A0; ideal_rna["  U  U  G1"]=A1; ideal_rna["  U  U  G2"]=A2;  ideal_rna["  U  U  G3"]=B0; ideal_rna["  U  U  G4"]=B1; ideal_rna["  U  U  G5"]=B2;
    ideal_rna["  U  U  U0"]=A0; ideal_rna["  U  U  U1"]=A1; ideal_rna["  U  U  U2"]=A2;  ideal_rna["  U  U  U3"]=B0; ideal_rna["  U  U  U4"]=B1; ideal_rna["  U  U  U5"]=B2;

    map<string,vector<double> >().swap(A0);
    map<string,vector<double> >().swap(B0);
    map<string,vector<double> >().swap(A1);
    map<string,vector<double> >().swap(B1);
    map<string,vector<double> >().swap(A2);
    map<string,vector<double> >().swap(B2);

    vector<vector<double> > xyz_list1(3,tmp);
    vector<vector<double> > xyz_list2(3,tmp);
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    for (map<string, map<string,vector<double> > >::iterator iter = ideal_rna.begin();
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
        
        for (map<string, vector<double> >::iterator i = ideal_rna[nt].begin();
            i != ideal_rna[nt].end(); i++)
        {
            string k=i-> first;
            if (ideal_rna[key].count(k)) continue;
            ideal_rna[key][k]=tmp;
            ChangeCoor(ideal_rna[nt][k], RotMatix, TranVect, ideal_rna[key][k]);
        }
    }

    vector<vector<double> >().swap(xyz_list1);
    vector<vector<double> >().swap(xyz_list2);
    vector<vector<double> >().swap(RotMatix);  // U
    vector<double> ().swap(TranVect);  // t
    vector<double> ().swap(tmp);
    return ideal_rna;
}
