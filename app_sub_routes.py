from flask import Flask, render_template, request, jsonify, redirect, url_for
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit import RDLogger 

import os 
import pandas as pd
import numpy as np
from FPSim2 import FPSim2Engine
from chemicalconverters import NamesConverter
import requests
import urllib.parse

# from chembl_webresource_client.new_client import new_client

from chembl_id_by_name import get_chembl_id_by_name
from chembl_id_by_smiles import get_chembl_id_by_smiles
from molecule_detail_by_id import get_molecule_details_by_id
from home_mol_name import get_name_by_smiles
from iupac_by_smiles import get_iupac_name_by_smiles
from predict_iupac_by_smiles import predict_iupac_name_by_smiles
from home_mol_name import get_home_mol_name

# https://www.ebi.ac.uk/chembl/api/data/status.json 
RDLogger.DisableLog('rdApp.*') 
# This silences warnings and error logs from RDKit

converter = NamesConverter(model_name="knowledgator/SMILES2IUPAC-canonical-base")
print("Hugging Face Model loaded: SMILES2IUPAC-canonical-base")

app = Flask(__name__)
app.secret_key = "chemsmile_ai_secure_key"
app.json.sort_keys = False 
# This line preserves the order of keys in the JSON output, which is important for maintaining the order of results in the response.

def get_chembl_client():
    """Safely attempts to initialize the ChEMBL client on-demand."""
    try:
        from chembl_webresource_client.new_client import new_client
        print(f"ChEMBL API client working now.")
        return new_client
    except Exception as e:
        print(f"ChEMBL API client initialization failed: Error getting schema from url https://www.ebi.ac.uk/chembl/api/data/spore with status 500")
        return None

client = get_chembl_client()

if client is None:
    # Fallback gracefully instead of crashing the app
    error = "ChEMBL database servers are down. Only raw SMILES lookups are supported right now."
else:
    # Run your ChEMBL API code normally here
    pass

# --- 1. PRELOAD DATA AT STARTUP (<1 second) ---
print("Initializing Instant Fingerprint Search Engine...")
# Load the pre-compiled binary database file
fpe = FPSim2Engine('molecules_library.h5')

# Load the source CSV metadata file
# We preserve the natural integer row index matching the database matrix
df_metadata = pd.read_csv(r'D:\ChemSmileAI\molecules_smiles.csv')
print("Search engine is hot and ready.")

# using hardware-accelerated FPSim2 matrix screening

#
# < - - - - - - - ChemSmileAI Home Page- - - - - - - - - >
#
@app.route("/", methods=["GET", "POST"])
def index():
    data = None
    error = None
    
    if request.method == "POST":
        smiles = request.form.get("smiles", "").strip()
        
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                error = "Invalid connection schema. Could not parse structure."
            else:
                data = {
                    "svg": render_molecule_svg(smiles),
                    "smiles": Chem.MolToSmiles(mol, canonical=True),
                    "formula": rdMolDescriptors.CalcMolFormula(mol),
                    "mw": round(Descriptors.MolWt(mol), 2),
                    "logp": round(Crippen.MolLogP(mol), 2),
                    "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
                    "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                    "lipinski": compute_lipinski(mol)
                }
        else:
            error = "Please provide a valid molecular token identifier."

    return render_template("index.html", data=data, error=error)

#
# < - - - - - - - Molecular Analysis NavBar Route - - - - - - - - - >
#
@app.route("/molecule_analysis", methods=["GET", "POST"])
def molecule_analysis():
    data = None
    error = None
    
    if request.method == "POST":
        smiles = request.form.get("smiles1", "").strip()
        
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                error = "Invalid connection schema. Could not parse structure."
            else:
                data = {
                    # Fumction to display the Molecule's Name on the Molecula ANalysis.html
                    "home_mol_name" : get_home_mol_name(smiles),
                    # Calling function for accessing the sub-functions
                    # 1. Molecular Property
                    "molecular_property" : molecular_property(mol)
                    # 2. Lipinski Property of 5
                    # "lipinski5": lipinski_5(mol),
                    #3. Molecular Fingerprint
                    # "morganfp": morgan_fp(mol),
                    #4. Atom/Bond Analysis
                    # "atom_analysis": atom_analysis(mol),
                    # "bond_anaylsis": bond_analysis(mol),
                    
                    # SVG Image
                    # "svg": render_molecule_svg(smiles)
                }
        else:
            error = "Please submit a valid SMILES structure above to query properties."

    return render_template("molecule_analysis.html", data=data, error=error)


def molecular_property(mol):
    smiles = Chem.MolToSmiles(mol)
    return{
        "compound_name": get_name_by_smiles(smiles),
        "common_name": "",
        "iupac_name": get_iupac_name_by_smiles(smiles),
        "predicted_iupac_name": predict_iupac_name_by_smiles(smiles),
        "predicted_chembl_id": get_chembl_id_by_smiles(smiles),
        "smiles": Chem.MolToSmiles(mol, canonical=True),
        "svg": render_molecule_svg(smiles),
        "elemental_formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "lipinski": compute_lipinski(mol),
        "summary": mol_summary(mol),
    }

def compute_lipinski(mol):
    """Calculates molecular metrics and evaluates Lipinski compliance."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    
    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    
    return "Pass" if violations <= 1 else f"Fail ({violations} Violations)"

def mol_summary(mol):
    return 'This is a brief summary of the molecule.'


@app.route("/molecule_analysis/lipinski", methods=["POST"])
def lipinski_5(mol):
    smiles = request.form.get("smiles1", "").strip()
    mol = Chem.MolFromSmiles(smiles)
    
    if not mol:
        return redirect(url_for('molecule_analysis'))
    
    data = {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "lipinski": compute_lipinski(mol)
    }
    return render_template("molecule_analysis.html", data=data, active_tab="lipinski")


@app.route("/molecule_analysis/fingerprint", methods=["POST"])
def morgan_fp(mol):
    smiles = request.form.get("smiles1", "").strip()
    mol = Chem.MolFromSmiles(smiles)
    
    if not mol:
        return redirect(url_for('molecule_analysis'))
    
    mf = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    morganfp = list(mf)
    morganfp_2048 = morganfp[:2049]
    onbits = len(mf.GetOnBits())
    idxonbits = tuple(mf.GetOnBits())
    data = {
        "morganfingerprint": morganfp_2048,
        "lenmf": len(mf),
        "onbits": onbits,
        "idxonbits": idxonbits
    }
    return render_template("molecule_analysis.html", data=data, active_tab="fingerprint")


@app.route("/molecule_analysis/atom", methods=["POST"])
def atom_analysis(mol):
    smiles = request.form.get("smiles1", "").strip()
    mol = Chem.MolFromSmiles(smiles)
    
    if not mol:
        return redirect(url_for('molecule_analysis'))
    
    atom_data_list = []
    
    for atom in mol.GetAtoms():
        # Extracted metrics formatted as clean key-value dictionary items
        data = {
            "index": atom.GetIdx(),
            "symbol": atom.GetSymbol(),
            "atomic_number": atom.GetAtomicNum(),
            "atomic_mass": round(atom.GetMass(), 3),
            "hybridization": str(atom.GetHybridization()),
            "charge": atom.GetFormalCharge(),
            "degree": atom.GetDegree(),
            "total_hs": atom.GetTotalNumHs(includeNeighbors=True),
            "explicit_valence": atom.GetExplicitValence(), 
            "implicit_valence": atom.GetImplicitValence(),
            # "explicit_valence" : atom.GetValence(Chem.ValenceType.EXPLICIT),
            # "implicit_valence" : atom.GetValence(Chem.ValenceType.IMPLICIT),
            # This is a new method, will not give errors.
            "aromatic": "Yes" if atom.GetIsAromatic() else "No",
            "chiral": str(atom.GetChiralTag()),
            "isotope": atom.GetIsotope()
        }
        atom_data_list.append(data)
    return render_template("molecule_analysis.html", data=data, active_tab="atom")

@app.route("/molecule_analysis/bond", methods=["POST"])
def bond_analysis(mol):
    smiles = request.form.get("smiles1", "").strip()
    mol = Chem.MolFromSmiles(smiles)
    
    if not mol:
        return redirect(url_for('molecule_analysis'))
    
    bond_data_list = []
    
    for bond in mol.GetBonds():
        # Extracted metrics formatted as clean key-value dictionary items
        data = {
            "bond": bond.GetBondTypeAsDouble(),
            "bond-type": bond.GetBondType(),
            "begin-atom-index": bond.GetBeginAtomIdx(),
            "end-atom-index": bond.GetEndAtomIdx(),
            "is-ring": bond.IsInRing()
        }
        bond_data_list.append(data)
        
    return render_template("molecule_analysis.html", data=data, active_tab="bond")


@app.route("/contact", methods=["GET", "POST"])
def contact():
    return render_template("contact.html")


@app.route("/similarity", methods=["GET", "POST"])
def similarity():
    return render_template("similarity.html")


@app.route("/about", methods=["GET", "POST"])
def about():
    return render_template("about.html")

@app.route("/compare", methods=["GET", "POST"])
def compare():
    return render_template("compare.html")


# Molecular Similarity Search Engine
@app.route('/search', methods=['GET'])
def search_molecule():
    """
    Looks up a target molecule by ChEMBL ID, maps it to its internal SMILES string, and runs a high-speed Tanimoto coefficient matrix similarity query.
    """
    target_id = request.args.get('chembl_id1', '').strip().upper()
    if not target_id:
        return jsonify({"error": "No ChEMBL ID provided"}), 400

    # 1. Fetch the metadata row matching the user's requested ChEMBL ID string
    target_row = df_metadata[df_metadata['chembl_id'] == target_id]
    
    if target_row.empty:
        return jsonify({"error": f"ID '{target_id}' not found in the dataset database."}), 404
        
    # 2. Extract its reference structural index position and canonical SMILES string
    target_int_id = int(target_row.index[0])
    target_smiles = str(target_row.iloc[0]['canonical_smiles'])

    try:
        # 3. Query FPSim2 via positional strings to isolate version discrepancies
        
        # Dynamically uses half of the available logical threads to keep the app responsive
        optimal_workers = max(1, os.cpu_count() // 2) 

        # Positional Order: (query_smiles, threshold, n_workers)
        results = fpe.similarity(target_smiles, 0.90, 'tanimoto', optimal_workers)

        # results is a list of tuples: [(row_index, similarity_score), ...]

        matches = []
        for match in results:
            match_int_id = int(match[0])     # Matched positional row index inside .h5 matrix
            score = round(float(match[1]), 6) # Calculated Tanimoto Similarity score
            
            # Skip self-comparison logic using row markers
            if match_int_id == target_int_id:
                continue
                
            # 4. Map the matching integer pointer back to its native metadata values
            real_chembl_id = df_metadata.iloc[match_int_id]['chembl_id']
            smiles = df_metadata.iloc[match_int_id]['canonical_smiles']
            
            matches.append({
                "chembl_id": real_chembl_id,
                "Tanimoto_Similarity": score,
                "canonical_smiles": smiles
            })
            
        # 5. Order results from highest structural affinity to lowest
        matches = sorted(matches, key=lambda x: x['Tanimoto_Similarity'], reverse=True)

        # print(results)

        return jsonify({
            "target_query": target_id,
            "target_smiles": target_smiles,
            "total_matches_found": len(matches),

            "results": matches
        })

    except Exception as e:
        return jsonify({"error": f"Search execution failed: {str(e)}"}), 500
    

# def render_molecule_svg(smiles):
#     """Parses a SMILES string and returns a cleanly structured inline SVG string."""
#     try:
#         mol = Chem.MolFromSmiles(smiles)
#         if not mol:
#             return ""
#         rdDepictor.Compute2DCoords(mol)
#         drawer = rdMolDraw2D.MolDraw2DSVG(350, 350)
#         drawer.DrawMolecule(mol)
#         drawer.FinishDrawing()
#         raw_svg = drawer.GetDrawingText()
#         svg_start = raw_svg.find("<svg")
#         return raw_svg[svg_start:] if svg_start != -1 else raw_svg
#     except Exception:
#         return ""


def render_molecule_svg(smiles):
    """Parses a SMILES string and returns a cleanly sliced inline SVG string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        raw_svg = drawer.GetDrawingText()
        
        # Completely remove any XML definitions or background metadata objects
        svg_start = raw_svg.find("<svg")
        if svg_start != -1:
            clean_svg = raw_svg[svg_start:]
            # Strip out transparent/white outer rects that force duplicate box clipping borders
            clean_svg = clean_svg.replace("<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300' height='300' x='0' y='0'> </rect>", "")
            return clean_svg
        return raw_svg
    except Exception:
        return ""


# def mol_properties()

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000) 

# flask run --host=127.0.0.1 --port=8000









