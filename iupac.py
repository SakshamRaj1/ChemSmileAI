import pandas as pd
import requests
import time

# Built a local CSV with IUPAC names for all known structures in your library.
# Now, IUPAC names can be retrieved instantly without any network calls.
# But it will only work for known structures. For novel structures, I still need to use the model-based approach.

# df = pd.read_csv("molecules_smiles.csv")
# iupac_names = []

# for idx, row in df.iterrows():
#     smiles = row['canonical_smiles']
#     try:
#         url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON"
#         res = requests.get(url, timeout=3).json()
#         name = res['PropertyTable']['Properties'][0]['IUPACName']
#         iupac_names.append(name)
#     except:
#         iupac_names.append("Unknown Compound Structure")
#     time.sleep(0.2) # Polite request spacing

# df['iupac_name'] = iupac_names
# df.to_csv("molecules_smiles_with_iupac_names.csv", index=False)

# Load at Flask startup
df_metadata = pd.read_csv('D:/ChemSmileAI/molecules_smiles_with_iupac_names.csv')

def get_local_iupac(chembl_id):
    # Match the ID instantly against your local text bank
    match = df_metadata[df_metadata['chembl_id'] == chembl_id]
    if not match.empty:
        return match['iupac_name'].values[0]
    return "Novel structure - Name not in local database"