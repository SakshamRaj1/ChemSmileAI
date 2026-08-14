# from chembl_webresource_client.new_client import new_client
# Imported in the app.py using try-except block to avoid app crash if CHEMBL API is down.
from chembl_webresource_client.http_errors import HttpApplicationError
import requests
import urllib.parse
import pubchempy as pcp
from predict_iupac_by_smiles import predict_iupac_name_by_smiles

# def get_chembl_client():
#     """Safely attempts to initialize the ChEMBL client on-demand."""
#     try:
#         from chembl_webresource_client.new_client import new_client
#         return new_client
#     except Exception as e:
#         print(f"ChEMBL API client initialization failed")
#         return None

def get_chembl_client():
    """Safely attempts to initialize the ChEMBL client on-demand."""
    try:
        from chembl_webresource_client.settings import Settings
        Settings.Instance().TIMEOUT = 0.5  # Set global HTTP timeout to 0.5s
        
        from chembl_webresource_client.new_client import new_client
        return new_client
    except Exception as e:
        print(f"ChEMBL API client initialization failed")
        return None
    

def get_iupac_name_by_smiles(smiles):
    safe_smiles = urllib.parse.quote(smiles)
    pubchem_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{safe_smiles}/property/IUPACName/JSON"
    )

    try:
        response = requests.get(pubchem_url, timeout=3)
        response.raise_for_status()

        data = response.json()
        return data["PropertyTable"]["Properties"][0]["IUPACName"]

    except Exception:
        return predict_iupac_name_by_smiles(smiles)

def get_molecule_details_by_id(chembl_id):
    """
    Given a ChEMBL ID, retrieves the preferred molecule name
    and its canonical SMILES string.
    """
    client = get_chembl_client()

    if client is None:
        # Fallback gracefully instead of crashing the app
        error = "ChEMBL database servers are down."
    else:
        # Running ChEMBL API code normally here
        pass

        try:
            molecule_client = client.molecule

            # Query using the exact ChEMBL ID
            # .only() limits the API response data to speed up the download
            res = molecule_client.filter(molecule_chembl_id=chembl_id).only(['pref_name', 'molecule_structures'])

            if res:
                record = res[0]

                # 1. Fetch Molecule Name
                # If the molecule doesn't have a common approved name, pref_name will be None
                name = record.get('pref_name') or "No preferred name listed"

                # 2. Fetch SMILES String
                # Structures are nested inside the 'molecule_structures' dictionary
                structures = record.get('molecule_structures')
                smiles = structures.get('canonical_smiles') if structures else None

                return {
                    'chembl_id': chembl_id,
                    'name': name,
                    'smiles': smiles if smiles else "No SMILES available"
                }
            else:
                print(f"Error: {chembl_id} not found in the ChEMBL database.")
                return None
            
        except HttpApplicationError as e:
                print(f"ChEMBL API is down.")
                return None

def get_chembl_id_by_smiles(smiles_string):
    smiles_string = smiles_string.replace("\r", "").replace("\n", "")

    client = get_chembl_client()

    if client is None:
        # Fallback gracefully instead of crashing the app
        error = "ChEMBL database servers are down."
        error
    else:
        # Run ChEMBL API code normally here
        pass

        molecule = client.molecule
        # Query using connectivity search for structural robustness
        res = molecule.filter(molecule_structures__canonical_smiles=smiles_string).only(['molecule_chembl_id'])

        try:
            if res and len(res) > 0:
                return res[0]['molecule_chembl_id']
            else:
                # Removed the print statement here
                return None
            
        except HttpApplicationError as e:
            print(f"ChEMBL API is down.")
            return None

        except Exception as e:
            return None

def get_name_by_smiles(smiles_string):
    client = get_chembl_client()
    
    if client is None:
        # Fallback gracefully instead of crashing the app
        print("ChEMBL database servers are down.")
  
    else:
        # Run ChEMBL API code normally here
        chembl_id = get_chembl_id_by_smiles(smiles_string)
        if chembl_id:
            details = get_molecule_details_by_id(chembl_id)
            print(f"Molecule details function said:, {details['name']}")
            return details['name'] if details else None    


def get_home_mol_name(smiles):
    # 1. Try fetching the ChEMBL common name first
    try:
        chembl_name = get_name_by_smiles(smiles)
    except Exception as e:
        print(f"Error fetching from ChEMBL: {e}")
        chembl_name = None

    invalid_chembl_responses = ["No chembl record found", "No preferred name listed", "Unknown Compound Structure", "None"]
    
    # Check if ChEMBL returned a valid name
    if chembl_name and chembl_name not in invalid_chembl_responses:
        return chembl_name.capitalize()

    # 2. If ChEMBL failed, try PubChemPy
    print(f"No valid ChEMBL name for {smiles}. Trying PubChemPy...")
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        
        if compounds:
            compound = compounds[0]
            if compound.synonyms:
                # Pick the first synonym that isn't empty or an error message
                compound_name = compound.synonyms[0]
                if compound_name and compound_name not in invalid_chembl_responses and not compound_name.lower().startswith("chembl"):
                    return compound_name.capitalize()
                    
        print(f"PubChem found a compound, but no clean name was listed.")
    except Exception as e:
        print(f"PubChem API failed: {e}")

    # 3. Final Fallback: IUPAC Name generation
    print("Fetching IUPAC name as a final fallback...")
    try:
        iupac_name = get_iupac_name_by_smiles(smiles)
        if not iupac_name:
            iupac_name = predict_iupac_name_by_smiles(smiles)
            
        if iupac_name:
            return iupac_name.capitalize()
    except Exception as e:
        print(f"IUPAC generation failed: {e}")

    return "Unknown Molecule"



# def get_home_mol_name(smiles):
#     # Try fetching the ChEMBL common name first
#     chembl_name = get_name_by_smiles(smiles)
#     if chembl_name != "No chembl record found" or chembl_name != "No preferred name listed" or chembl_name != None:
#         if chembl_name == "No chembl record found" or chembl_name == "No preferred name listed" or None:
#             print("Chembl returned", chembl_name, "for Molecule name Home page")
#             # ChEMBL might return None or a specific string depending on the API wrapper
#             print(f"No ChEMBL name for {smiles} in CHEMBL Database. Trying from PubChempy")
    
#             try:
#                 # Search PubChem by SMILES
#                 compounds = pcp.get_compounds(smiles, namespace='smiles')
                
#                 if not compounds:
#                     return "No matching molecule found in PubChem."
                    
#                 # Get the first matching compound object
#                 compound = compounds[0]
                
#                 # Return the first common name/synonym if available; otherwise, return IUPAC name
#                 if compound.synonyms:
#                     compound_name = compound.synonyms[0]
#                     print(f"Found molecule (CID: {compound.cid}), but no name is listed.")
#                     if not compound_name or compound_name == "No record found" or compound_name == "No preferred name listed":
#                         print("Chembl, Pubchem API failed to fetch Compound name, showing IUPAC name instead")
#                         iupac_name = get_iupac_name_by_smiles(smiles) or predict_iupac_name_by_smiles(smiles)
#                         return iupac_name.capitalize() if iupac_name else "Unknown Molecule"
#                     else:
#                         return compound_name.capitalize()
                    
#             except Exception as e:
#                 print("Chembl, Pubchem API failed to fetch Compound name, showing IUPAC name instead")
#                 print(f"str{e}")
#                 return None
                
#         # return chembl_name.capitalize()

#     else:
#         iupac_name = get_iupac_name_by_smiles(smiles)
#         return iupac_name.capitalize() if iupac_name else chembl_name.capitalize() or "Unknown Molecule"
    

# def get_name_from_pubchempy(mol):

#     # Query compound by its SMILES string
#     compounds = pcp.get_compounds('CC(=O)OC1=CC=CC=C1C(=O)O', 'smiles')

#     for comp in compounds:
#         print("CID:", comp.cid)
#         print("Name:", comp.synonyms[0] if comp.synonyms else "No name found")