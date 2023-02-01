# Arena (Atomic Reconstruction of RNA)

## Instructions for running

Download Arena from GitHub into your desired location, such as your Desktop:
```
cd Desktop/
git clone git@github.com:pylelab/Arena.git
```

In the new Arena folder, compile with Clang:
```
cd Arena
make Arena
```

Run Arena using one of the provided examples:
```
./Arena Examples/5osg2.input.pdb Examples/5osg2.output.pdb 5
```

Calculate RMSD:
```
python calculate_RMSD.py -r Examples/5osg2.pdb -i Examples/5osg2.output.pdb
```

Calculate number of clashes:
```
make Arena_counter
./Arena_counter Examples/5osg2.output.pdb Examples/5osg2.clashes.txt 5
```
