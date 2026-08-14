import requests
import urllib.parse

def get_iupac_name_by_smiles(smiles):
    # Calls the PubChem API to retrieve the IUPAC name for a given SMILES string. 
    safe_smiles = urllib.parse.quote(smiles)
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{safe_smiles}/property/IUPACName/JSON"
    
    try:
        response = requests.get(pubchem_url, timeout=3)
        if response.status_code == 200:
            data = response.json()
            iupac = data['PropertyTable']['Properties'][0]['IUPACName']
            print(f"IUPAC name from PubCHem", iupac)
            return iupac
    except Exception:
        return "[Unknown Compound Structure]"
    
