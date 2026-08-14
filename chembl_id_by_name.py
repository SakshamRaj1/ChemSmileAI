# from chembl_webresource_client.new_client import new_client
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

    
def get_chembl_id_by_name(molecule_name):
    # Query ChEMBL for the specific preferred name

    client = get_chembl_client()

    if client is None:
        # Fallback gracefully instead of crashing the app
        error = "ChEMBL database servers are down. Only raw SMILES lookups are supported right now."
    else:
        # Run your ChEMBL API code normally here
        pass

        try:
            molecule = client.molecule
            res = molecule.filter(pref_name__iexact=molecule_name).only(['molecule_chembl_id'])

            if res:
                # Return the first match's ID
                return res[0]['molecule_chembl_id']
            else:
                print(f"No ChEMBL record found for name: {molecule_name}")
                return None
            
        except HttpApplicationError as e:
            print(f"ChEMBL API is down: {e}")
            return None