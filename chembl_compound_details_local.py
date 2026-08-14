import csv
import pandas as pd
from rdkit import Chem

df_metadata = pd.read_csv(r'chembl_smiles_prefname_synonyms.csv', dtype = {"pref_name": "string", "all_synonyms": "string"})

# df_metadata = df_metadata.where(pd.notna(df_metadata), None)

def get_chemblid_from_smiles_locally(smiles):
    # Check if the SMILES exists in the local DataFrame and retrieve the 
    # corresponding ChEMBL ID if exists else "Novel Query"

    row = df_metadata[df_metadata['canonical_smiles'] == smiles]

    if not row.empty:
        chembl_id = str(row.iloc[0]['chembl_id']) 
        # print(f"Found ChEMBL ID for SMILES {smiles}: {chembl_id}")
        return chembl_id
    else:
        return "NOVEL_QUERY"


def get_name_from_smiles_locally(smiles):
    # Check if the SMILES exists in the local DataFrame and retrieve the 
    # corresponding synonyms and preferred name if exists else "No preferred name found"

    row = df_metadata[df_metadata['canonical_smiles'] == smiles]

    if not row.empty:
        all_synonyms = row.iloc[0]['all_synonyms']
        pref_name = row.iloc[0]['pref_name']
        print("Preferred Name: ", pref_name)
        print("All Names: ", all_synonyms)
        names = {
            "pref_name": pref_name,
            "all_synonyms": all_synonyms            
        }
        return names
    else:
        return "No preferred name found"


def get_pref_name_from_smiles_locally(smiles):
    # Check if the SMILES exists in the local DataFrame and retrieve the 
    # corresponding synonyms and preferred name if exists else "No preferred name found"

    row = df_metadata[df_metadata['canonical_smiles'] == smiles]

    if not row.empty:
        pref_name = row.iloc[0]['pref_name']
        pref_name = "No preferred name found" if pd.isna(pref_name) else pref_name
        # handling pandas NA in dataframe
        print("Preferred Name: ", pref_name)
        names = {
            "pref_name": pref_name,          
        }
        return names
    else:
        return "No preferred name found"




def get_smiles_from_chemblid_locally(chembl_id):
    # Check if the ChEMBL ID exists in the local DataFrame and retrieve the
    # corresponding SMILES if exists else "No preferred name found"

    row = df_metadata[df_metadata['chembl_id'] == chembl_id]

    if not row.empty:
        smiles = row.iloc[0]['canonical_smiles']
        return smiles
    else:
        return "No preferred SMILES found"


def get_name_from_chemblid_locally(chembl_id):
    # Check if the ChEMBL ID exists in the local DataFrame and retrieve the
    # corresponding synonyms and preferred name if exists else "No preferred name found"

    row = df_metadata[df_metadata['chembl_id'] == chembl_id]

    if not row.empty:
        all_synonyms = row.iloc[0]['all_synonyms']
        pref_name = row.iloc[0]['pref_name']
        return all_synonyms, pref_name
    else:
        return "No preferred name found"




def get_smiles_from_name_locally(name):
    # Check if the name (synonym) exists in the local DataFrame and retrieve the
    # corresponding synonyms and preferred name if exists else "No preferred name found"

    row1 = df_metadata[df_metadata['pref_name'] == name]
    # row2 = df_metadata[df_metadata['all_synonyms'] == smiles]

    if not row1.empty:
        smiles = row1.iloc[0]['canonical_smiles']
        return smiles
    else:
        return "No preferred SMILES found"


def get_chemblid_from_name_locally(name):
    # Check if the name (synonym) exists in the local DataFrame and retrieve the
    # corresponding ChEMBL ID if exists else "Novel Query"

    row1 = df_metadata[df_metadata['pref_name'] ==  name]
    # row2 = df_metadata[df_metadata['all_synonyms'] == smiles]

    if not row1.empty:
        chembl_id = str(row1.iloc[0]['chembl_id']) 
        return chembl_id
    else:
        return "NOVEL_QUERY"
    


def get_all_details_from_chemblid_locally(clean_str):
    # Single function to get all details from chemblid
    chembl_id = clean_str.strip()
    row = df_metadata[df_metadata['chembl_id'] == chembl_id]

    if not row.empty:
        smiles = row.iloc[0]['canonical_smiles'] if smiles is not None else "None"
        pref_name = row.iloc[0]['pref_name'] if pref_name is not None else "None"
        all_synonyms = row.iloc[0]['all_synonyms'] if all_synonyms is not None else "None"

        detailsfromchemblid = {
            "chemblid": chembl_id,
            "smiles": smiles,
            "pref_name": pref_name,
            "all_synonyms": all_synonyms
        }
        return detailsfromchemblid
    else:
        return "Invalid CHEMBL ID"


def get_all_details_from_smiles_locally(clean_str):
    # Single function to get all details from SMILES
       
    smiles = clean_str.strip()
    if Chem.MolFromSmiles(smiles) is True:
        print("invalid smiles")
    row = df_metadata[df_metadata['canonical_smiles'] == smiles]
    
    if not row.empty:
        chembl_id = str(row.iloc[0]['chembl_id']) if not row.empty else "No CHEMBL ID exists"
        # chembl_id = "No CHEMBL ID exists" if pd.isna(chembl_id) is not None else chembl_id
        pref_name = row.iloc[0]['pref_name'] 
        pref_name = "No preferred name found" if pd.isna(pref_name) else pref_name
        all_synonyms = row.iloc[0]['all_synonyms']
        all_synonyms = "No preferred name found" if pd.isna(all_synonyms) else all_synonyms

        detailsfromsmiles = {
            "chemblid": chembl_id,
            "smiles": smiles,
            "pref_name": pref_name,
            "all_synonyms": all_synonyms
        }
        return detailsfromsmiles
    # If smiles does not exists, it may be a correct smiles but not present in CHEMBL DB
    else:
        detailsfromsmiles = {
            "chemblid": "No CHEMBL ID exists", 
            "smiles": smiles,
            "pref_name": "No preferred name found",             
            "all_synonyms": "No preferred name found"}
        return detailsfromsmiles


def get_all_details_from_name_locally(clean_str):
    # Single function to get all details from name
    name = clean_str.strip()
    row1 = df_metadata[df_metadata['pref_name'] == name]
    row2 = df_metadata[df_metadata['all_synonyms'] == name]
    
    if not row1.empty or not row2.empty:
        chembl_id = row1.iloc[0]['chembl_id'] if chembl_id is not None else "None"
        smiles = row1.iloc[0]['canonical_smiles'] if smiles is not None else "None"
        pref_name = row1.iloc[0]['pref_name'] if pref_name is not None else "None"
        all_synonyms = row1.iloc[0]['all_synonyms'] if all_synonyms is not None else "None"

        detailsfromname = {
            "chemblid": chembl_id,
            "smiles": smiles,
            "pref_name": pref_name,
            "all_synonyms": all_synonyms
        }
        return detailsfromname
    else:
        return "Invalid CHEMBL ID"








def one_function_to_identify_input(input_str):
    clean_str = input_str.strip()
    print("Input for 1F_CHEMBL search engine:", clean_str)
    check_str = clean_str.upper()
    # 1. Check for ChEMBL ID
    if clean_str[6:].isdigit() or check_str.startswith("CHEMBL"): 
        try:
            clean_str = clean_str.upper()
            clean_str = "".join(clean_str.split()) # allowing chembl 1000 to be accepted despite having a space in between
            print(f"1F_CHEMBL identified input as a format of CHEMBL ID:", clean_str)
            data = get_all_details_from_chemblid_locally(clean_str)
            return data
        
        except:
            # error = "This CHEMBL ID does not exists or Invalid CHEMBL ID"
            print(f"1F_CHEMBL search engine received invalid input:", clean_str)
            return "Invalid CHEMBL ID"

    # 2. Check for SMILES String
    # Common structural syntax characters and lowercase aromatic elements
    smiles_indicators = {'=', '#', '@', '(', ')', '[', ']', 'c', 'n', 'o', 's'}

    # If it contains SMILES symbols or contains numbers (for rings)
    has_smiles_chars = any(char in smiles_indicators for char in clean_str)
    has_digits = any(char.isdigit() for char in clean_str)

    if has_smiles_chars or has_digits:
        try:
            mol = Chem.MolFromSmiles(clean_str)
            if mol is None:
                error = "Invalid SMILES string structure"
                return False
            
            print(f"1F_CHEMBL identified input as a format of SMILES:", clean_str)
            data = get_all_details_from_smiles_locally(clean_str)
            return data

        except Exception as e:
            print(f"1F_CHEMBL search engine received invalid input:", {clean_str})

    else:
        if clean_str.isalpha() or " " in clean_str:
            try:
                print(f"1F_CHEMBL identified input as a (could be preferred name, synonym, IUPAC, or common name):", clean_str)
                data = get_all_details_from_name_locally(clean_str)
                return data

            except Exception as e:
                print(f"1F_CHEMBL search engine received invalid/unknown input type:", {clean_str})




