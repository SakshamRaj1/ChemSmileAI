def identify_string_no_re(input_str):
    # Standardise input by stripping spaces
    clean_str = input_str.strip()
    
    # 1. Check for ChEMBL ID
    if clean_str.startswith("CHEMBL") and clean_str[6:].isdigit():
        return "ChEMBL ID"
    
    # 2. Check for SMILES String
    # Common structural syntax characters and lowercase aromatic elements
    smiles_indicators = {'=', '#', '@', '(', ')', '[', ']', 'c', 'n', 'o', 's'}
    
    # If it contains SMILES symbols or contains numbers (for rings)
    has_smiles_chars = any(char in smiles_indicators for char in clean_str)
    has_digits = any(char.isdigit() for char in clean_str)
    
    if has_smiles_chars or has_digits:
        return "SMILES String"
    
    # 3. Default to Chemical Name
    return "Chemical Name"

# Verification
print(identify_string_no_re("CHEMBL25"))               # ChEMBL ID
print(identify_string_no_re("CC(=O)Oc1ccccc1C(=O)O"))  # SMILES String
print(identify_string_no_re("Aspirin"))                # Chemical Name
