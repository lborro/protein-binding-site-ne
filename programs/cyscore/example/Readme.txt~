CYSCORE v2.0.0


Cyscore is an empirical scoring function for accurate protein-ligand binding affinty prediction. It is compromised of hydrophobic free energy, van der Waals interaction energy and the ligand's entropy. To improve the prediction accuracy, a novel curvature weighted surface area model was developed for the hydrophobic free energy calculation. Our analyses show that the new model is superior to the conventional surface area model indeed, implying that surface shape is also important other than surface area for the prediction of hydrophobic free energy. Cyscore will be a useful tool for the accurate protein-ligand binding affinity prediction.

Cyscore is free for non-commercial users. Other users please contact cy_scu@yeah.net .
Redistribution is NOT allowed without explicit permission.
The program is distributed in the hope that they will be useful,but WITHOUT ANY WARRANTY.

Please cite "Improved protein-ligand binding affinity prediction by using a curvature dependent surface area model. Yang Cao and Lei Li,  Bioinformatics (2014) 30 (12): 1674-1680. " 

Cyscore is a command-line application under Linux-x86 (32-bit Intel/AMD) The usage is :

     Cyscore [Protein PDB file] [Ligand MOL2 or PDB file]

Please note: All hydrogens of proteins and ligands are needed for correct calculation.

Improvements:
1. Fixed a bug in assigning van der Waals radius to ligand atoms.
2. Using more efficient surface area model, which is 4 times faster compared to the old version.
3. Reduced the Entropy term which has 7 or more rotatable bonds.


Example:

$../bin/Cyscore 3nw9_protein.pdb 3nw9_ligand.mol2

-----------------------------------------------------
 Cyscore:  Protein-Ligand Binding Affinity Predictor
       Yang Cao  2016 All Rights Reserved.
                      V 2.0.0
           Sat Apr 23 23:54:24 2016
-----------------------------------------------------

Protein: 3nw9_protein.pdb    Ligand: 3nw9_ligand.mol2
Hydrophobic -2.3241 Vdw -3.9487 HBond 0.0000 Ent 0.2520
Cyscore = -6.0208


Note:  The lower the predicted value, the better the binding affinity. 




