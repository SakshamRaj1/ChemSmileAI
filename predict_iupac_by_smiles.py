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

  
