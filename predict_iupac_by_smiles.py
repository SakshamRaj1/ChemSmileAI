from openclatura import name_smiles
# from chemicalconverters import NamesConverter
import requests
import urllib.parse
from rdkit import Chem

def predict_iupac_name_by_smiles(smiles):
    """
    Predict IUPAC name from a SMILES string using Openclatura.

    Returns:
        str: Predicted IUPAC name, or
             "No preferred IUPAC name found" if unavailable.
    """
    try:
        # Ensure the Openclatura function is available
        if not callable(name_smiles):
            return "No preferred IUPAC name found"

        # Predict IUPAC name
        iupac_name = name_smiles(smiles)
        print(f"Predicted IUPAC name from Openclatura: {iupac_name}")

        # Handle invalid or missing predictions
        if (
            iupac_name is None
            or str(iupac_name).strip() == ""
            or str(iupac_name).strip() == "No preferred IUPAC name found"
        ):
            return "No preferred IUPAC name found"

        return iupac_name

    except Exception as e:
        print(f"Error predicting IUPAC name: {e}")
        return "No preferred IUPAC name found"

    
#  Removing any kind of dependency on PubChem API, so created a new function above
# def predict_iupac_name_by_smiles(smiles):
#     # -----------------------------
#     # Step 1: Use Openclature library to get IUPAC name from SMILES
#     # -----------------------------
#     if name_smiles:
#         iupac_name = name_smiles(smiles)
#         print(f"Predicted IUPAC name from Openclatura is: , {iupac_name}")
#         if iupac_name is None:
#             pass
#             # URL encode the SMILES
#             safe_smiles = urllib.parse.quote(smiles)
#             pubchem_url = (
#                 f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{safe_smiles}/property/IUPACName/JSON"
#             )

#             iupac_name1 = "[Unknown Compound Structure]"

#             # ----------------------------------------------------------
#             # Step 2: Try PubChem lookup if Openclatura fails
#             # ----------------------------------------------------------
#             try:
#                 response = requests.get(pubchem_url, timeout=3)

#                 if response.status_code == 200:
#                     data = response.json()
#                     iupac_name1 = data["PropertyTable"]["Properties"][0]["IUPACName"]

#             except Exception:
#                 # Ignore API/network errors and continue with fallback
#                 pass

#             # ----------------------------------------------------------
#             # Step 3: If PubChem found a name,
#             # return it immediately if it is not = "Unknown Compound Structure"
#             # ----------------------------------------------------------
#             print("Showing IUPAC name from PubChem after Openclatura failed")

#             if iupac_name1 != "[Unknown Compound Structure]" or iupac_name1 != "No preferred name listed" or iupac_name != None:
#                 print(f"IUPAC name from PubChem after Openclatura failed: , {iupac_name1}")
#                 return iupac_name1
        
#             # -----------------------------
#             # Step 4: In the case of problem in Pubchem or Openclature, return IUPAC name from Openclatura only
#             # -----------------------------
#         else:
#             return iupac_name
                
#             # -----------------------------
#             # Step 4: In the case of problem in Pubchem or Openclature, return IUPAC name from PubChem
#             # -----------------------------

#     else:
#         print(f"Predicted IUPAC name from PubChem also failed, Openclatura said: , {iupac_name}")
#         return iupac_name or "No preferred name predicted"
        




# Skipping the 'chemicalconverters' method as it was unable to predict the smiles having 
# bond directions like / or \ or @ or all these symbols


# converter = NamesConverter(model_name="knowledgator/SMILES2IUPAC-canonical-base")
# print("Hugging Face Model loaded: SMILES2IUPAC-canonical-base")

# def predict_iupac_name_by_smiles(smiles):
#     # URL encode the SMILES
#     safe_smiles = urllib.parse.quote(smiles)
#     pubchem_url = (
#         f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{safe_smiles}/property/IUPACName/JSON"
#     )

#     pub_iupac = "[Unknown Compound Structure]"

#     # -----------------------------
#     # Step 1: Try PubChem lookup
#     # -----------------------------
#     try:
#         response = requests.get(pubchem_url, timeout=3)

#         if response.status_code == 200:
#             data = response.json()
#             pub_iupac = data["PropertyTable"]["Properties"][0]["IUPACName"]

#     except Exception:
#         # Ignore API/network errors and continue with fallback
#         pass

#     # -----------------------------
#     # Step 2: If PubChem found a name,
#     # return it immediately
#     # -----------------------------
#     if pub_iupac != "[Unknown Compound Structure]":
#         return pub_iupac

#     # -----------------------------
#     # Step 3: PubChem didn't know it
#     # -----------------------------
#     if converter is None:
#         return pub_iupac

#     if len(smiles) > 400:
#         print("Structure too large for neural model translation")
#         return pub_iupac

#     # -----------------------------
#     # Step 4: Neural model fallback
#     # -----------------------------
#     try:
#         mol = Chem.MolFromSmiles(smiles)
#         canonical = Chem.MolToSmiles(mol, canonical=True)
#         prediction = converter.smiles_to_iupac(f"<BASE>{canonical}")
#         return prediction[0] if isinstance(prediction, list) else prediction

#     except Exception:
#         # If converter also fails, return the PubChem result
#         return pub_iupac