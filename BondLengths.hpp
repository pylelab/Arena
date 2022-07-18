/* Purpose: Optimize bond lengths
Ideal bond lengths are taken from the Amber forcefield parameters.
The order of each array along both rows and columns is the standard PDB order, with hydrogen atoms omitted.
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
float A_bonds[22][22] = {
	{0, 1.48, 1.48, 1.598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 1.438, 0, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 1.521, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.521, 0, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.521, 0, 0, 1.418, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1.418, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 1.438, 0, 0, 1.521, 0, 0, 1.467, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.467, 0, 1.371, 0, 0, 0, 0, 0, 0, 0, 1.359},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.371, 0, 1.279, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.279, 0, 1.368, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.368, 0, 1.403, 0, 0, 0, 0, 1.373},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.403, 0, 1.363, 1.296, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.363, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.296, 0, 0, 1.319, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.319, 0, 1.319, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.319, 0, 1.342},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.359, 0, 0, 1.373, 0, 0, 0, 0, 1.342, 0}
};

// C
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6
float C_bonds[20][20] = {
	{0, 1.48, 1.48, 1.598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 1.438, 0, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 1.521, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.521, 0, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.521, 0, 0, 1.418, 1.521, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1.418, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 1.438, 0, 0, 1.521, 0, 0, 1.467, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.467, 0, 1.386, 0, 0, 0, 0, 0, 1.372},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.386, 0, 1.229, 1.358, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.229, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.358, 0, 0, 1.296, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.296, 0, 1.363, 1.442, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.363, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.442, 0, 0, 1.344},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.372, 0, 0, 0, 0, 0, 1.344, 0}
};

// G
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4
float G_bonds[23][23] = {
	{0, 1.48, 1.48, 1.598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 1.438, 0, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 1.521, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.521, 0, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.521, 0, 0, 1.418, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1.418, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 1.438, 0, 0, 1.521, 0, 0, 1.467, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.467, 0, 1.371, 0, 0, 0, 0, 0, 0, 0, 0, 1.359},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.371, 0, 1.279, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.279, 0, 1.368, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.368, 0, 1.439, 0, 0, 0, 0, 0, 1.373},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.439, 0, 1.229, 1.389, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.229, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.389, 0, 0, 1.361, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.361, 0, 1.363, 1.296, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.363, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.296, 0, 0, 1.342},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.359, 0, 0, 1.373, 0, 0, 0, 0, 0, 1.342, 0}
};

// U
//Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6
float U_bonds[20][20] = {
	{0, 1.48, 1.48, 1.598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1.598, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 1.438, 0, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 1.521, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 1.521, 0, 0, 1.438, 1.521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 1.521, 0, 0, 1.418, 1.521, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1.418, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 1.438, 0, 0, 1.521, 0, 0, 1.467, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.467, 0, 1.386, 0, 0, 0, 0, 0, 1.372},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.386, 0, 1.229, 1.389, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.229, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.389, 0, 0, 1.389, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.389, 0, 1.229, 1.466, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.229, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.466, 0, 0, 1.344},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.372, 0, 0, 0, 0, 0, 1.344, 0}
};

// Bond between O3' and P of the next residue
float OSP = 1.610;

/* check is residue r and residue r+1 is connected */
bool isConnected(const ModelUnit &pep, const int c, const int r)
{
    for (int a=0;a<=12;a++) if (Points2Distance(
        pep.chains[c].residues[r].atoms[0].xyz,
        pep.chains[c].residues[r+1].atoms[0].xyz)<=8) return true;
    return false;
}

int fix_bond_lengths (ModelUnit &pep, float tolerance){
	
	// Initialize variables
	float b_ideal = 0;
	float b = 0;
	float deltaB = 0;
	float diff = 0;
	int a2 = 0;

	float a1x = 0;
	float a1y = 0;
	float a1z = 0;
	float a2x = 0;
	float a2y = 0;
	float a2z = 0;

	int moved = 0;

	// Iterate through chains
	for (int c=0; c<pep.chains.size(); c++){
		// Iterate through residues
		for (int r=0; r<pep.chains[c].residues.size(); r++){
			string resname = pep.chains[c].residues[r].resn;
			int resnumber = pep.chains[c].residues[r].resi;
			
			/* Special case to check the O3'-P bonds */
			b_ideal = OSP;
			// Get (x,y,z) coordinates
			a1x = pep.chains[c].residues[r].atoms[8].xyz[0];
			a1y = pep.chains[c].residues[r].atoms[8].xyz[1];
			a1z = pep.chains[c].residues[r].atoms[8].xyz[2];
			if (r != (pep.chains[c].residues.size()-1) && isConnected(pep,c,r)){
				//cout << resname << resnumber << pep.chains[c].residues[r].atoms[8].name << pep.chains[c].residues[r+1].atoms[0].name << endl;
				a2x = pep.chains[c].residues[r+1].atoms[0].xyz[0];
				a2y = pep.chains[c].residues[r+1].atoms[0].xyz[1];
				a2z = pep.chains[c].residues[r+1].atoms[0].xyz[2];
				// Calculate the current bond length
				b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
				// Calculate fractional difference between the current and ideal bond lengths
				diff = abs(b - b_ideal) / b_ideal;
				// Check if difference is greater than or equal to the allowed tolerance
				if (diff >= tolerance){
					deltaB = b - b_ideal; // difference in bond length from ideal
					// Check which atoms are movable and assign new coordinates to the movable atom
					if (pep.chains[c].residues[r].atoms[8].movable==1 && pep.chains[c].residues[r+1].atoms[0].movable==0){
						//cout << "a1 movable " << resname << resnumber << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
						moved++;
						pep.chains[c].residues[r].atoms[8].xyz[0] = (deltaB/b)*(a2x-a1x) + a1x;
						pep.chains[c].residues[r].atoms[8].xyz[1] = (deltaB/b)*(a2y-a1y) + a1y;
						pep.chains[c].residues[r].atoms[8].xyz[2] = (deltaB/b)*(a2z-a1z) + a1z;
					} else if (pep.chains[c].residues[r].atoms[8].movable==0 && pep.chains[c].residues[r+1].atoms[0].movable==1){
						//cout << "a2 movable " << resname << resnumber << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
						moved++;
						pep.chains[c].residues[r+1].atoms[0].xyz[0] = (deltaB/b)*(a1x-a2x) + a2x;
						pep.chains[c].residues[r+1].atoms[0].xyz[1] = (deltaB/b)*(a1y-a2y) + a2y;
						pep.chains[c].residues[r+1].atoms[0].xyz[2] = (deltaB/b)*(a1z-a2z) + a2z;
					} else if (pep.chains[c].residues[r].atoms[8].movable==1 && pep.chains[c].residues[r+1].atoms[0].movable==1){
						//cout << "both movable " << resname << resnumber << "b: " << b << "b_ideal " << b_ideal << "diff: " << diff << endl;
						moved++;
						// move each atom half the required distance
						pep.chains[c].residues[r+1].atoms[0].xyz[0] = 0.5*((deltaB/b)*(a1x-a2x)) + a2x;
						pep.chains[c].residues[r+1].atoms[0].xyz[1] = 0.5*((deltaB/b)*(a1y-a2y)) + a2y;
						pep.chains[c].residues[r+1].atoms[0].xyz[2] = 0.5*((deltaB/b)*(a1z-a2z)) + a2z;
						pep.chains[c].residues[r].atoms[8].xyz[0] = 0.5*((deltaB/b)*(a2x-a1x)) + a1x;
						pep.chains[c].residues[r].atoms[8].xyz[1] = 0.5*((deltaB/b)*(a2y-a1y)) + a1y;
						pep.chains[c].residues[r].atoms[8].xyz[2] = 0.5*((deltaB/b)*(a2z-a1z)) + a1z;
					}
				}
			}
			/* */

			// bool check_OSP = false; // check O3'-P bond
			// Iterate through all pairs of backbone atoms
			for (int a1=0; a1<pep.chains[c].residues[r].atoms.size(); a1++){
				string atom1 = pep.chains[c].residues[r].atoms[a1].name;
				//cout << "atom1" << endl;
				// int r2 = r;
				// if (a1 == 8){
				// 	r2 = r + 1; // set r2 equal to the next residue for the O3'-P bond
				// }
				for (a2=a1+1; a2<pep.chains[c].residues[r].atoms.size(); a2++){
					//cout << "atom2" << endl;
					string atom2 = pep.chains[c].residues[r].atoms[a2].name;
					bool pass = false;
					// if neither atom is movable, immediately continue
					if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==0){
						//cout << "none movable " << resname << resnumber << atom1 << atom2 << endl;
						continue;
					}
					//cout << atom1 << pep.chains[c].residues[r].atoms[a1].movable << endl;
					// Check if atoms are connected and, if so, save their ideal bond length
					//cout << "find atom in array" << endl;
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
					// if (! pass && atom1==" O3'" && check_OSP==false){
					// 	pass = true;
					// 	check_OSP=true;
					// 	b_ideal = OSP;
						// cout << atom1 << endl;
					//cout<<"resname="<<resname<<". a1="<<a1<<".a2="<<a2<<".pass="<<pass<<endl;
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
					diff = abs(b - b_ideal) / b_ideal;
					// Check if difference is greater than or equal to the allowed tolerance
					if (diff >= tolerance){
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
							//if (b == 0){ // if atoms exactly overlapping
								// generate random direction
							//}
					}	
				}
			}
		}
	}

	return moved;
}
#endif
