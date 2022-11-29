# Elife '22 smarts fingerprint

This code computes the 42-element substructure fingerprint, as in \[1\].

## Requirements
1. Numpy
1. RDKit: `pip install rdkit` or `conda install -c rdkit rdkit`

## Quick start
```
import fingerprinting as fp
from rdkit import Chem

smiles  = [
	'CC(=O)OC',            #ethyl acetate
	'CCOC(=O)/C(=C/C)/C',  #ethyl tiglate
	'CC1=CCC2C(C1)C2(C)C'] #carene

mols = [Chem.MolFromSmiles(s)  for s in smiles]

res = fp.get_fingerprint(mols)

```

## Limitations of the feature set
This set of fingerprints captures well the chemical diversity of the 185 molecules used in [1]. Other datasets may require extending the dictionary of SMARTS patterns to achieve fill coverage of chemical space. Alternatively, consider MACCS keys, which worked well in our experience. MACCS keys can be computed with RDKIT \[2\].


## References
[1] Shawn D Burton, Audrey Brown, Thomas P Eiting, Isaac A Youngstrom, Thomas C Rust, Michael Schmuker, Matt Wachowiak (2022) Mapping odorant sensitivities reveals a sparse but structured representation of olfactory chemical space by sensory input to the mouse olfactory bulb eLife 11:e80470, [https://doi.org/10.7554/eLife.80470](https://doi.org/10.7554/eLife.80470)

[2] http://rdkit.org/docs/GettingStartedInPython.html?highlight=maccs#maccs-keys