# # from chembl_webresource_client.new_client import new_client
# Imported in the app.py using try-except block to avoid app crash if CHEMBL API is down.
from chembl_webresource_client.http_errors import HttpApplicationError

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

    
def get_chembl_id_by_smiles(smiles_string):
    smiles_string = smiles_string.replace("\r", "").replace("\n", "")

    client = get_chembl_client()

    if client is None:
        # Fallback gracefully instead of crashing the app
        error = "ChEMBL database servers are down. Only raw SMILES lookups are supported right now."
    else:
        # Run your ChEMBL API code normally here
        pass

        molecule = client.molecule
        # Query using connectivity search for structural robustness
        res = molecule.filter(molecule_structures__canonical_smiles=smiles_string).only(['molecule_chembl_id'])

        try:
            if res and len(res) > 0:
                return res[0]['molecule_chembl_id']
            else:
                return "No ChEMBL record found for SMILES"
            
        except HttpApplicationError as e:
            print(f"ChEMBL API is down: {e}")
            return None
        
        except Exception as e:
            # Catches server-side conversion issues or network dropouts
            # print(f"API Error: Could not compute structural identifiers on ChEMBL backend. Details: {e}"
            return f"No ChEMBL record found for SMILES: {smiles_string}"
