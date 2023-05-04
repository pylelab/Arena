# Arena (Atomic Reconstruction of RNA)

Arena is a stand-alone program used for the full-atomic reconstruction of RNA structures from coarse-grained models. Arena can build in full-atomic detail for any residue with at least one heavy atom (of any identity). Arena takes an input in PDB format and first fills in the missing heavy atoms via superimposition. Arena iterates through three refinement steps to correct bond lengths and angles, correct base and base pair conformations, and remove atomic clashes. Structures generated by Arena are highly accurate in terms of RMSD (< 3.7 Å on average), contain almost no clashes, and run in < 3.5 sec on average. The runtime will depend on the length of your RNA; the longest RNA we tested was 692 nt and ran in < 35 sec. Arena can be downloaded from here and run on your local computer or you can use the webserver at [https://zhanggroup.org/Arena/](https://zhanggroup.org/Arena/). We have provided example PDB structures to run with Arena in the "Examples" folder.

## Installation

Download Arena from GitHub into your desired location, such as your Desktop:
```
git clone https://github.com/pylelab/Arena.git
```

In the new Arena folder, compile with Clang:
```
cd Arena
make Arena
```

## Usage

Run Arena using one of the provided examples:
```
./Arena Examples/5osg2.input.pdb Examples/5osg2.output.pdb 5
```

To see the help menu:
```
./Arena
```
Note: If the RNA is interacting with other molecules (proteins, ligands, etc), the user must manually remove the non-RNA chains for the Arena standalone program. Arena can, however, handle multichain RNA inputs without user intervention (see below). On the other hand, [the Arena webserver](https://zhanggroup.org/Arena/) automatically excludes non-RNA residues.

## Special Case: Multiple Models

Split into separate model files before running Arena:
```
python split_models.py Examples/2koc.pdb
```

## Benchmarking Scripts

Calculate RMSD:
```
python calculate_RMSD.py -r Examples/5osg2.pdb -i Examples/5osg2.output.pdb
```

Calculate number of clashes:
```
make Arena_counter
./Arena_counter Examples/5osg2.output.pdb Examples/5osg2.clashes.txt 5
```

## Reference
Zion R Perry, Anna Marie Pyle, Chengxin Zhang (2023)
Arena: rapid and accurate reconstruction of full atomic RNA structures from coarse-grained models.

