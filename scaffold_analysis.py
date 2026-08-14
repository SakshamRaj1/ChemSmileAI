from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# Define target molecule using a SMILES string
smiles_string = "Cc1ccccc1C(=O)N2CCOCC2" 
molecule = Chem.MolFromSmiles(smiles_string)

# Extracting the basic Murcko Scaffold 
scaffold_molecule = MurckoScaffold.GetScaffoldForMol(molecule)
scaffold_smiles = Chem.MolToSmiles(scaffold_molecule)

print("Original Molecule SMILES:", smiles_string)
print("Extracted Core Scaffold SMILES:", scaffold_smiles)
