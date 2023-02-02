/* Purpose: Refine bond lengths and angles
Ideal bond and pseudobond lengths are taken from the Amber forcefield parameters.
*/

#ifndef BondLengths_HPP
#define BondLengths_HPP 1

#include <math.h>
#include <stdbool.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include "GeometryTools.hpp"

using namespace std;

// A
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4
double A_bonds[22][22] = {
	{0, 1.48, 1.48, 1.598, 2.637, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 2.410, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 2.410, 0, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 2.532, 2.532, 0, 1.438, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{2.637, 0, 0, 1.438, 0, 1.521, 2.406, 2.425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 2.406, 1.521, 0, 1.438, 1.521, 2.406, 2.425, 0, 2.343, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.406, 1.438, 0, 2.406, 0, 2.406, 0, 1.438, 2.371, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.425, 1.521, 2.406, 0, 1.438, 1.521, 2.410, 2.425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.406, 0, 1.438, 0, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.425, 2.406, 1.521, 2.406, 0, 1.418, 1.521, 2.484, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 2.410, 0, 1.418, 0, 2.410, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.343, 1.438, 2.425, 0, 1.521, 2.410, 0, 1.467, 2.551, 0, 0, 0, 0, 0, 0, 0, 2.522},
	{0, 0, 0, 0, 0, 0, 2.371, 0, 0, 2.484, 0, 1.467, 0, 1.371, 2.223, 2.187, 0, 0, 0, 0, 2.411, 1.359},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.551, 1.371, 0, 1.279, 2.092, 0, 0, 0, 0, 0, 2.174},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.223, 1.279, 0, 1.368, 2.548, 0, 0, 0, 0, 2.250},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.187, 2.092, 1.368, 0, 1.403, 2.433, 2.316, 0, 2.433, 1.373},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.548, 1.403, 0, 1.363, 1.296, 2.256, 0, 2.353},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.433, 1.363, 0, 2.272, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.316, 1.296, 2.272, 0, 1.319, 2.375, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.256, 0, 1.319, 0, 1.319, 2.201},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.411, 0, 0, 2.433, 0, 0, 2.375, 1.319, 0, 1.342},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.522, 1.359, 2.174, 2.250, 1.373, 2.353, 0, 0, 2.201, 1.342, 0}
};

// C
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6
double C_bonds[20][20] = {
	{0, 1.48, 1.48, 1.598, 2.637, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 2.410, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 2.410, 0, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 2.532, 2.532, 0, 1.438, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{2.637, 0, 0, 1.438, 0, 1.521, 2.406, 2.425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 2.406, 1.521, 0, 1.438, 1.521, 2.406, 2.425, 0, 2.343, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.406, 1.438, 0, 2.406, 0, 2.406, 0, 1.438, 2.371, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.425, 1.521, 2.406, 0, 1.438, 1.521, 2.410, 2.425, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.406, 0, 1.438, 0, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.425, 2.406, 1.521, 2.406, 0, 1.418, 1.521, 2.484, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 2.410, 0, 1.418, 0, 2.410, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.343, 1.438, 2.425, 0, 1.521, 2.410, 0, 1.467, 2.434, 0, 0, 0, 0, 0, 2.474},
	{0, 0, 0, 0, 0, 0, 2.371, 0, 0, 2.484, 0, 1.467, 0, 1.386, 2.284, 2.346, 0, 0, 2.391, 1.372},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.434, 1.386, 0, 1.229, 1.358, 2.310, 0, 0, 2.403},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.284, 1.229, 0, 2.290, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.346, 1.358, 2.290, 0, 1.296, 2.272, 2.409, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.310, 0, 1.296, 0, 1.363, 1.442, 2.352},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.272, 1.363, 0, 2.424, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.391, 0, 0, 2.409, 1.442, 2.424, 0, 1.344},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.474, 1.372, 2.403, 0, 0, 2.352, 0, 1.344, 0}
};

// G
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4
double G_bonds[23][23] = {
	{0, 1.48, 1.48, 1.598, 2.637, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 2.410, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 2.410, 0, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 2.532, 2.532, 0, 1.438, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{2.637, 0, 0, 1.438, 0, 1.521, 2.406, 2.425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 2.406, 1.521, 0, 1.438, 1.521, 2.406, 2.425, 0, 2.343, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.406, 1.438, 0, 2.406, 0, 2.406, 0, 1.438, 2.371, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.425, 1.521, 2.406, 0, 1.438, 1.521, 2.410, 2.425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.406, 0, 1.438, 0, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.425, 2.406, 1.521, 2.406, 0, 1.418, 1.521, 2.484, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 2.410, 0, 1.418, 0, 2.410, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.343, 1.438, 2.425, 0, 1.521, 2.410, 0, 1.467, 2.551, 0, 0, 0, 0, 0, 0, 0, 0, 2.522},
	{0, 0, 0, 0, 0, 0, 2.371, 0, 0, 2.484, 0, 1.467, 0, 1.371, 2.223, 2.187, 0, 0, 0, 0, 0, 2.411, 1.359},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.551, 1.371, 0, 1.279, 2.092, 0, 0, 0, 0, 0, 0, 2.174},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.223, 1.279, 0, 1.368, 2.539, 0, 0, 0, 0, 0, 2.250},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.187, 2.092, 1.368, 0, 1.439, 2.413, 2.319, 0, 0, 2.433, 1.373},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.539, 1.439, 0, 1.229, 1.389, 2.475, 0, 0, 2.432},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.413, 1.229, 0, 2.276, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.319, 1.389, 2.276, 0, 1.361, 2.317, 2.342, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.475, 0, 1.361, 0, 1.363, 1.296, 2.204},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.317, 1.363, 0, 2.272, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.411, 0, 0, 2.433, 0, 0, 2.342, 1.296, 2.272, 0, 1.342},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.522, 1.359, 2.174, 2.250, 1.373, 2.432, 0, 0, 2.204, 0, 1.342, 0}
};

// U
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6
double U_bonds[20][20] = {
	{0, 1.48, 1.48, 1.598, 2.637, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 2.410, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 2.410, 0, 2.532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 2.532, 2.532, 0, 1.438, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{2.637, 0, 0, 1.438, 0, 1.521, 2.406, 2.425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 2.406, 1.521, 0, 1.438, 1.521, 2.406, 2.425, 0, 2.343, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.406, 1.438, 0, 2.406, 0, 2.406, 0, 1.438, 2.371, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 2.425, 1.521, 2.406, 0, 1.438, 1.521, 2.410, 2.425, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.406, 0, 1.438, 0, 2.406, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.425, 2.406, 1.521, 2.406, 0, 1.418, 1.521, 2.484, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 2.410, 0, 1.418, 0, 2.410, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 2.343, 1.438, 2.425, 0, 1.521, 2.410, 0, 1.467, 2.434, 0, 0, 0, 0, 0, 2.474},
	{0, 0, 0, 0, 0, 0, 2.371, 0, 0, 2.484, 0, 1.467, 0, 1.386, 2.284, 2.335, 0, 0, 2.391, 1.372},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.434, 1.386, 0, 1.229, 1.389, 2.497, 0, 0, 2.403},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.284, 1.229, 0, 2.276, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.335, 1.389, 2.276, 0, 1.389, 2.276, 2.392, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.497, 0, 1.389, 0, 1.229, 1.466, 2.428},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.276, 1.229, 0, 2.399, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.391, 0, 0, 2.392, 1.466, 2.399, 0, 1.344},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.474, 1.372, 2.403, 0, 0, 2.428, 0, 1.344, 0}
};

// OSP
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3'
double OSP_bonds[9][4] = {
	{0, 0, 0, 0},
	{0, 0, 0, 0},
	{0, 0, 0, 0},
	{0, 0, 0, 0},
	{0, 0, 0, 0},
	{0, 0, 0, 0},
	{0, 0, 0, 0},
	{2.648, 0, 0, 0},
	{1.610, 2.505, 2.505, 2.504}
};

// Bond between O3' and P of the next residue
double OSP = 1.610;

/* Check if residue r and residue r+1 are connected */
bool isConnected(const ModelUnit &pep, const int c, const int r)
{
    for (int a=0;a<=12;a++) if (Points2Distance(
        pep.chains[c].residues[r].atoms[a].xyz,
        pep.chains[c].residues[r+1].atoms[a].xyz)<=8) return true;
    return false;
}

int fix_bond_lengths (ModelUnit &pep, double tolerance){
	
	// Initialize variables
	double b_ideal = 0;
	double b = 0;
	double deltaB = 0;
	double diff = 0;
	int a2 = 0;

	double a1x = 0;
	double a1y = 0;
	double a1z = 0;
	double a2x = 0;
	double a2y = 0;
	double a2z = 0;

	int moved = 0;

	// Iterate through chains
	for (int c=0; c<pep.chains.size(); c++){
		// Iterate through residues
		for (int r=0; r<pep.chains[c].residues.size(); r++){
			string resname = pep.chains[c].residues[r].resn;
			int resnumber = pep.chains[c].residues[r].resi;
			
			/* Special case to check the O3'-P bonds */
			if (r != (pep.chains[c].residues.size()-1) && isConnected(pep,c,r)) {
                for (int a1=7;a1<=8;a1++) {
			        // Get (x,y,z) coordinates
			        a1x = pep.chains[c].residues[r].atoms[a1].xyz[0];
			        a1y = pep.chains[c].residues[r].atoms[a1].xyz[1];
			        a1z = pep.chains[c].residues[r].atoms[a1].xyz[2];
                    for (a2=0;a2<=3;a2++) {
                        b_ideal = OSP_bonds[a1][a2];
                        if (b_ideal<=0) continue;
					    a2x = pep.chains[c].residues[r+1].atoms[a2].xyz[0];
					    a2y = pep.chains[c].residues[r+1].atoms[a2].xyz[1];
					    a2z = pep.chains[c].residues[r+1].atoms[a2].xyz[2];
					    // Calculate the current bond length
					    b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					    // Calculate fractional difference between the current and ideal bond lengths
					    diff = fabs(b - b_ideal) / b_ideal;
					    // Check if difference is greater than or equal to the allowed tolerance
					    if (diff >= tolerance && b>0.0001){
						    //cout << resname << resnumber << atom1 << atom2 << b << ' ' << b_ideal << ' ' << diff <<endl;
						    deltaB = b - b_ideal; // difference in bond length from ideal
						    // Check which atoms are movable and assign new coordinates to the movable atom
						    if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==0){
							    //cout << "a1 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							    moved++;
							    pep.chains[c].residues[r].atoms[a1].xyz[0] = (deltaB/b)*(a2x-a1x) + a1x;
							    pep.chains[c].residues[r].atoms[a1].xyz[1] = (deltaB/b)*(a2y-a1y) + a1y;
							    pep.chains[c].residues[r].atoms[a1].xyz[2] = (deltaB/b)*(a2z-a1z) + a1z;
						    } else if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==1){
							    //cout << "a2 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							    moved++;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[0] = (deltaB/b)*(a1x-a2x) + a2x;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[1] = (deltaB/b)*(a1y-a2y) + a2y;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[2] = (deltaB/b)*(a1z-a2z) + a2z;
						    } else if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==1){
							    //cout << "both movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							    moved++;
							    // move each atom half the required distance
							    pep.chains[c].residues[r+1].atoms[a2].xyz[0] = 0.5*((deltaB/b)*(a1x-a2x)) + a2x;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[1] = 0.5*((deltaB/b)*(a1y-a2y)) + a2y;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[2] = 0.5*((deltaB/b)*(a1z-a2z)) + a2z;
							    pep.chains[c].residues[r].atoms[a1].xyz[0] = 0.5*((deltaB/b)*(a2x-a1x)) + a1x;
							    pep.chains[c].residues[r].atoms[a1].xyz[1] = 0.5*((deltaB/b)*(a2y-a1y)) + a1y;
							    pep.chains[c].residues[r].atoms[a1].xyz[2] = 0.5*((deltaB/b)*(a2z-a1z)) + a1z;
						    }
					    }	
					}
				}
			}

			// Iterate through all pairs of backbone atoms
			for (int a1=0; a1<pep.chains[c].residues[r].atoms.size(); a1++){
				string atom1 = pep.chains[c].residues[r].atoms[a1].name;
				for (a2=a1+1; a2<pep.chains[c].residues[r].atoms.size(); a2++){
					string atom2 = pep.chains[c].residues[r].atoms[a2].name;
					bool pass = false;
					// if neither atom is movable, immediately continue
					if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==0){
						//cout << "none movable " << resname << resnumber << atom1 << atom2 << endl;
						continue;
					}
					//cout << atom1 << pep.chains[c].residues[r].atoms[a1].movable << endl;
					//cout << "find atom in array" << endl;
					
					// Get ideal bond lengths
					if (resname == "  A"){
						if (A_bonds[a1][a2]>0){
							pass = true;
							b_ideal = A_bonds[a1][a2];
						}
					} else if (resname == "  C"){
						if (C_bonds[a1][a2]>0){
							pass = true;
							b_ideal = C_bonds[a1][a2];
						}
					} else if (resname == "  G"){
						if (G_bonds[a1][a2]>0){
							pass = true;
							b_ideal = G_bonds[a1][a2];
						}
					} else if (resname == "  U"){
						if (U_bonds[a1][a2]>0){
							pass = true;
							b_ideal = U_bonds[a1][a2];
						}
					}

					if (! pass) continue;
					// Get (x,y,z) coordinates
					a1x = pep.chains[c].residues[r].atoms[a1].xyz[0];
					a1y = pep.chains[c].residues[r].atoms[a1].xyz[1];
					a1z = pep.chains[c].residues[r].atoms[a1].xyz[2];
					a2x = pep.chains[c].residues[r].atoms[a2].xyz[0];
					a2y = pep.chains[c].residues[r].atoms[a2].xyz[1];
					a2z = pep.chains[c].residues[r].atoms[a2].xyz[2];
					// Calculate the current bond length
					b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					// Calculate fractional difference between the current and ideal bond lengths
					diff = fabs(b - b_ideal) / b_ideal;
					// Check if difference is greater than or equal to the allowed tolerance
					if (diff >= tolerance && b>0.0001){
						//cout << resname << resnumber << atom1 << atom2 << b << ' ' << b_ideal << ' ' << diff <<endl;
						deltaB = b - b_ideal; // difference in bond length from ideal
						// Check which atoms are movable and assign new coordinates to the movable atom
						if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==0){
							//cout << "a1 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							moved++;
							pep.chains[c].residues[r].atoms[a1].xyz[0] = (deltaB/b)*(a2x-a1x) + a1x;
							pep.chains[c].residues[r].atoms[a1].xyz[1] = (deltaB/b)*(a2y-a1y) + a1y;
							pep.chains[c].residues[r].atoms[a1].xyz[2] = (deltaB/b)*(a2z-a1z) + a1z;
						} else if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==1){
							//cout << "a2 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							moved++;
							pep.chains[c].residues[r].atoms[a2].xyz[0] = (deltaB/b)*(a1x-a2x) + a2x;
							pep.chains[c].residues[r].atoms[a2].xyz[1] = (deltaB/b)*(a1y-a2y) + a2y;
							pep.chains[c].residues[r].atoms[a2].xyz[2] = (deltaB/b)*(a1z-a2z) + a2z;
						} else if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==1){
							//cout << "both movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							moved++;
							// move each atom half the required distance
							pep.chains[c].residues[r].atoms[a2].xyz[0] = 0.5*((deltaB/b)*(a1x-a2x)) + a2x;
							pep.chains[c].residues[r].atoms[a2].xyz[1] = 0.5*((deltaB/b)*(a1y-a2y)) + a2y;
							pep.chains[c].residues[r].atoms[a2].xyz[2] = 0.5*((deltaB/b)*(a1z-a2z)) + a2z;
							pep.chains[c].residues[r].atoms[a1].xyz[0] = 0.5*((deltaB/b)*(a2x-a1x)) + a1x;
							pep.chains[c].residues[r].atoms[a1].xyz[1] = 0.5*((deltaB/b)*(a2y-a1y)) + a1y;
							pep.chains[c].residues[r].atoms[a1].xyz[2] = 0.5*((deltaB/b)*(a2z-a1z)) + a1z;
						}
					}	
				}
			}
		}
	}

	return moved;
}
#endif
