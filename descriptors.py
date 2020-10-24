# A script to calculate molecular descriptors and fingerprints based on molecule SMILES
# Created by: svavil
# Created on: 2020-10-24
# Edited on: 2020-10-25

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

data_block_1 = pd.read_excel("bcl-xl_1.xlsx").drop(labels = ["Unnamed: 6", "Unnamed: 7", "Unnamed: 8"], axis = "columns")
#reference_fingerprint = Chem.RDKFingerprint(Chem.MolFromSmiles("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C(C4=C(C(=C23)O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O.Cl"))
doxrub_fingerprint = Chem.RDKFingerprint(Chem.MolFromSmiles("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O"))

#data_block_1["reference_similarity"] = data_block_1.SMILES.apply(lambda s: DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(Chem.MolFromSmiles(s)), reference_fingerprint))
data_block_1["doxrub_similarity"] = data_block_1.SMILES.apply(lambda s: DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(Chem.MolFromSmiles(s)), doxrub_fingerprint))
data_block_1["exact_mol_weight"] = data_block_1.SMILES.apply(lambda s: Descriptors.ExactMolWt(Chem.MolFromSmiles(s)))
data_block_1["ring_count"] = data_block_1.SMILES.apply(lambda s: Descriptors.RingCount(Chem.MolFromSmiles(s)))

data_block_1["fingerprint"] = data_block_1.SMILES.apply(lambda s: list(Chem.RDKFingerprint(Chem.MolFromSmiles(s)).GetOnBits()))
fingerprint_block = data_block_1[["fingerprint"]].explode("fingerprint").assign(value = 1).pivot(columns = "fingerprint")
fingerprint_block.columns = fingerprint_block.columns.get_level_values(1)
fingerprint_block = fingerprint_block.add_prefix("fp.").fillna(value = 0)
data_block_1 = pd.concat([data_block_1, fingerprint_block], axis = "columns").drop(labels = "fingerprint", axis = "columns")

data_block_1["pharmacophores"] = data_block_1.SMILES.apply(lambda s: list(Generate.Gen2DFingerprint(Chem.MolFromSmiles(s), Gobbi_Pharm2D.factory).GetOnBits()))
pharmacophores_block = data_block_1[["pharmacophores"]].explode("pharmacophores").assign(value = 1).pivot(columns = "pharmacophores")
pharmacophores_block.columns = pharmacophores_block.columns.get_level_values(1)
pharmacophores_block = pharmacophores_block.add_prefix("ph.").fillna(value = 0)
data_block_1 = pd.concat([data_block_1, pharmacophores_block], axis = "columns").drop(labels = "pharmacophores", axis = "columns")

data_block_1.to_excel("bcl-xl_1_descriptors.xlsx", index = False) 