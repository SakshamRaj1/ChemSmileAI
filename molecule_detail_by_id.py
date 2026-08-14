# from chembl_webresource_client.new_client import new_client
# Imported in the app.py using try-except block to avoid app crash if CHEMBL API is down.
from chembl_webresource_client.http_errors import HttpApplicationError

def get_chembl_client():
    """Safely attempts to initialize the ChEMBL client on-demand."""
    try:
        from chembl_webresource_client.new_client import new_client
        return new_client
    except Exception as e:
        print(f"ChEMBL API client initialization failed")
        return None


def get_molecule_details_by_id(chembl_id):
    """
    Given a ChEMBL ID, retrieves the preferred molecule name
    and its canonical SMILES string.
    """
    client = get_chembl_client()
    
    if client is None:
        # Fallback gracefully instead of crashing the app
        error = "ChEMBL database servers are down. Only raw SMILES lookups are supported right now."
    else:
        # Run your ChEMBL API code normally here
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