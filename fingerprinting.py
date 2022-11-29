## Fingerprints
# (C) Michael Schmuker
# Substructure fingerprints for odorants
# please cite this git repo and also [1] when using it
#
# Shawn D Burton, Audrey Brown, Thomas P Eiting, Isaac A Youngstrom, 
# Thomas C Rust, Michael Schmuker, Matt Wachowiak (2022). Mapping 
# odorant sensitivities reveals a sparse but structured representation 
# of olfactory chemical space by sensory input to the mouse olfactory 
# bulb. 
# eLife 11:e80470 
# https://doi.org/10.7554/eLife.80470

# Works similarly to MACCS keys in RDKit. 

from rdkit import Chem
import numpy as np


odorant_patterns_elife22 = {
     0: ('*-C(=O)-[OH1]', 'carboxylic acid'),
     1: ('[CH1]=O', 'aldehyde'),
     2: ('C-C(=O)-[O]-C', 'ester'),
     3: ('C-C(=O)-[S]-C', 'thioester'),
     4: ('[!O&!S]-C(=O)-[!O&!S]', 'ketone'),
     5: ('[OX2H][CX4&!$(C([OX2H])[O,S,#7,#15]),c]', 'alcohol'),
     6: ('c1ccccc1', 'benzyl'),
     7: ('C~C(~C)~C~C~C~C(~C)~C', 'monoterpene'),
     8: ('[#8]1~[#6]~[#6]~[#6]~[#6]1', 'furanoid'),
     9: ('o1cccc1', 'furan'),
     10: ('[NH2][C]', 'primary amine'),
     11: ('[NH](C)C', 'secondary amine'),
     12: ('[NH0](C)(C)C', 'tertiary amine'),
     13: ('[N,n]1~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]1', 'pyridine'),
     14: ('[n,N]1~[C,c]~[C,c]~[C,c]~[C,c]1', 'pyrrole'),
     15: ('[N,n]1~[C,c]~[C,c]~[N,n]~[C,c]~[C,c]1', 'pyrazine'),
     16: ('[#16]1~[#6]~[#7]~[#6]~[#6]1', 'thiazoline'),
     17: ('[!#8]~C-S-C~[!#8]', 'thioether'),
     18: ('[$(C-S-S-C),$(C-S-S-S-C)]', 'sulfide'),
     19: ('[#6]-[SH]', 'thiol'),
     20: ('[#6]=[#6]', 'Alkene'),
     21: ('[#16]', 'sulfur'),
     22: ('[#7]', 'nitrogen'),
     23: ('[#8]', 'oxygen'),
     24: ('[R]', 'ring'),
     25: ('[CH3]-*-[CH2]-*', 'terminal 4-bond chain'),
     26: ('*!@*@*!@*', 'ortho-substituted ring'),
     27: ('*!@*@*@*!@*', 'meta-substituted rings'),
     28: ('*1(!@*)@*@*@*(!@*)@*@*@1',
      'para-substituted 6-ring\n(but not fused ring)'),
     29: ('C~C(~C)~[R1]1~[R1]~[R1]~[R1](~C)~[R1]~[R1]~1', 'menthane scaffold'),
     30: ('C~C(~C)~2~[R2]1~[R2]~2~[R1]~[R1](~C)~[R1]~[R1]~1', 'carene scaffold'),
     31: ('C~C(~C)~[R2]12~[R1]~[R2]~2~[R1](~C)~[R1]~[R1]~1', 'thujane scaffold'),
     32: ('C~C2(~C)~[R]1~[R]~[R]~2~[R](~C)~[R]~[R]~1', 'pinane scaffold'),
     33: ('[!H]~[!H]2(~[!H])~[R]1~[R]~[R]~[R](~[!H])~2~[R]~[R]~1',
      'camphane scaffold'),
     34: ('[!H]~[!H]2(~[!H])~[R]~[R](~[!H])1~[R]~[R]~2~[R]~[R]~1',
      'fenchane scaffold'),
     35: ('C(-C)(-C)(-C)-C', 'quadra C'),
     36: ('C-C-C-C-C-C', 'six C single bond'),
     37: ('C-C-C-C-C-C-C', 'seven C single bond'),
     38: ('C-C-C-C-C-C-C-C', 'eight C single bond'),
     39: ('C-C-C-C-C-C-C-C-C', 'nine C single bond'),
     40: ('C-C-C-C-C-C-C-C-C-C', 'ten C single bond'),
     41: ('C-C-C-C-C-C-C-C-C-C-C', 'eleven C single bond')
 }

## Todo: generate Chem.MolFromSmarts once upon first run of get_fingerprint,
##       e.g. by encapsulating the computation into an object. Should speed 
##       up the procedure considerably for very large molecule collections.


def get_fingerprint(mols, patterns=odorant_patterns_elife22):
    """
    compute substructure matching fingerprint using a dictionary of SMARTS.
    mols: iterable of rdkit molecules
    patterns: dictionary with SMARTS patterns. 
        Format: {<int>:('<SMARTS>', '<name>')}, e.g. {19:('[#6]-[SH]','thiol'), ...}
    """
    num_pats = len(patterns.keys())
    fp_mat = np.zeros((len(mols),num_pats),dtype=float)
    for i,m in enumerate(mols):
        if not m:
            a = np.ones(num_pats)
            a.fill(np.nan)
            fp_mat[i,:] = a
            continue 
        for j in range(num_pats):
            pattern = patterns[j][0]
            smarts = Chem.MolFromSmarts(pattern)
            if m.HasSubstructMatch(smarts):
                fp_mat[i,j]=1.

    return fp_mat


