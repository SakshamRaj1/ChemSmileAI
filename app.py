from flask import Flask, render_template, request, session, jsonify, send_file, redirect, url_for
from flask_caching import Cache
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, rdShapeHelpers , DataStructs, rdFMCS, rdMolAlign
from rdkit.Chem import AllChem, GraphDescriptors, Descriptors3D, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdDepictor
from rdkit import RDLogger 
import py3Dmol

import os 
import csv
import json
import pandas as pd
import numpy as np
from FPSim2 import FPSim2Engine

from chembl_compound_details_local import get_all_details_from_smiles_locally, get_pref_name_from_smiles_locally

from mol_desc_text import describe_molecule
from predict_iupac_by_smiles import predict_iupac_name_by_smiles
from download_section import download_data_csv, download_data_sdf, download_data_mol, download_data_mol2
from download_section import save_as_2d_sdf, save_as_2d_mol, save_as_2d_mol2, download_pdf, download_svg, download_png
from download_section import save_as_3d_sdf, save_as_3d_mol, save_as_3d_mol2, save_as_json, save_as_html

from pathlib import Path
import gdown

# gdown.download(
#     "https://drive.google.com/file/d/1mJkjwsLKaOKHwcRlbmUi0wz2vqjXWJcT/",
#     "chembl_smiles_prefname_synonyms.csv",
#     quiet=False
# )

# gdown.download(
#     "https://drive.google.com/file/d/1oBadN_FES73IIPKpJTYOVfBXvlBDxhGe/",
#     "molecules_library.h5",
#     quiet=False
# )


file_id1 = '1mJkjwsLKaOKHwcRlbmUi0wz2vqjXWJcT'
output_path1 = 'chembl_smiles_prefname_synonyms.csv'
file_id2 = '1oBadN_FES73IIPKpJTYOVfBXvlBDxhGe'
output_path2 = 'molecules_library.h5'

gdown.download(id=file_id1, output=output_path1, quiet=False)
gdown.download(id=file_id2, output=output_path2, quiet=False)

# https://www.ebi.ac.uk/chembl/api/data/status.json 
RDLogger.DisableLog('rdApp.*') 
# This silences warnings and error logs from RDKit

app = Flask(__name__)
app.secret_key = "chemsmile_ai_secure_key"
app.json.sort_keys = False 
# This line preserves the order of keys in the JSON output, which is important for maintaining the order of results in the response.

print("ChemSmileAI is launching!")

# Configure simple in-memory cache
cache = Cache(app, config={'CACHE_TYPE': 'SimpleCache', 'CACHE_DEFAULT_TIMEOUT': 3600})


# --- 1. PRELOAD DATA AT STARTUP (<1 second) ---
print("Initializing Instant Fingerprint Search Engine...")
# Load the pre-compiled binary database file
fpe = FPSim2Engine('molecules_library.h5')

# Load the source CSV metadata file
# Preserving the natural integer row index matching the database matrix

df_metadata = pd.read_csv(r'chembl_smiles_prefname_synonyms.csv', dtype = {"pref_name": "string", "all_synonyms": "string"})

print("Search engine is hot and ready.")

print("ChemSmileAI is working for you!")

# using hardware-accelerated FPSim2 matrix screening

#
# < - - - - - - -- - - - - - - - - - - - - - ChemSmileAI Home Page - - - - - - - - - - - - - - - - - - - - - - - >
#

def structure3d(smiles):
    if not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        result = AllChem.EmbedMolecule(mol, params)

        if result != 0:
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)

        if result != 0:
            return None

        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)

        mol_block = Chem.MolToMolBlock(mol)

        viewer = py3Dmol.view(width="100%", height="80%")
        viewer.addModel(mol_block, "mol")

        viewer.setBackgroundColor("black", 0)

        viewer.setStyle({
            "stick": {"radius": 0.15},
            "sphere": {"scale": 0.25}
        })

        viewer.zoomTo()

        n = mol.GetNumAtoms()

        if n <= 10:
            viewer.zoom(3.5)
        elif n <= 20:
            viewer.zoom(1.8)
        elif n <= 30:
            viewer.zoom(1.2)
        elif n <= 40:
            viewer.zoom(0.9)
        else:
            viewer.zoom(0.6)

        viewer.spin("y", 0.5)

        return viewer._make_html()

    except Exception as e:
        print(e)
        return None


@app.route("/", methods=["GET", "POST"])
def index():
    data = None
    error = None
    
    if request.method == "POST":
        smiles = request.form.get("smiles", "").strip() or request.args.get("smiles", "").strip()
        print("Home page received input:", smiles)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                error = "Invalid connection schema. Could not parse structure."
                print(error)
                error = "Invalid SMILES"
            else:
                predicted_iupac = predict_iupac_name_by_smiles(smiles)
                compound_details = get_all_details_from_smiles_locally(smiles)
                all_synonyms = compound_details['all_synonyms']
                compound_details['all_synonyms'] = [name.strip() for name in all_synonyms.split(",") if name.strip()]
                describe_molecule_text = describe_molecule(smiles)
                atom_count = mol.GetNumAtoms()
                structure3dim = structure3d(smiles) if atom_count <= 30 else None

                data = {
                    "svg": render_molecule_svg(smiles),
                    "compound_name": compound_details['pref_name'],
                    "all_synonyms": compound_details['all_synonyms'],
                    "predicted_chembl_id": compound_details['chemblid'],        
                    "structure3d": structure3dim,
                    "smiles": Chem.MolToSmiles(mol, canonical=True),
                    "formula": rdMolDescriptors.CalcMolFormula(mol),
                    "mw": round(Descriptors.MolWt(mol), 2),
                    "logp": round(Crippen.MolLogP(mol), 2),
                    "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
                    "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                    "lipinski": compute_lipinski(mol), 
                    "iupac_name": predicted_iupac, #later I will use describe_molecule_text to fetch iupac also
                    "describe_molecule_parent": describe_molecule_text['parent'],
                    "describe_molecule_features": describe_molecule_text['principal_feature']
                }
        else:
            error = "Please provide a valid molecular token identifier."

    return render_template("index.html", data=data, error=error)


@app.route("/receive_smiles_for_home", methods=["POST"])
# route to send the smile data to the terminal from ketcher
def receive_smiles_for_home():
    data = request.get_json()
    # print(data)
    print(f"Smiles from ketcher for Home page is {data["smiles"]}")

    return {"status": "ok"}




#
# < - - - - - - - - - - - - - - - - - - - - - Molecular Analysis NavBar Route - - - - - - - - - - - - - - - - - - - - - - - >
#

@app.route("/analysis", methods=["GET", "POST"])
def analysis():
    return render_template("molecule_analysis.html")

@app.route("/molecule_analysis", methods=["GET", "POST"])
def molecule_analysis():
    data = None
    error = None
    
    # if request.method == "POST":
    smiles = request.form.get("smiles1", "").strip() or request.args.get("smiles", "").strip()
    smarts = request.form.get("smarts", "").strip() # Fetch from Substructure search input card

    # Checks if SMILES is present.
    data = cache.get(smiles)

    # Smiles drawn on ketcher or typed in searchbar and received after form submission
    print("Received SMILES after form submission:", smiles)

    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            error = "Invalid connection schema. Could not parse structure."
        else:
            # Substruct search (runs if smarts is provided)
            substruct_result = substruct_search(mol, smarts) if smarts else {"present": False}
            compound_name= get_all_details_from_smiles_locally(smiles)
            predicted_iupac_name = predict_iupac_name_by_smiles(smiles)

            home_name = compound_name['pref_name'] if compound_name['pref_name'] != "No preferred name found" else predicted_iupac_name if predicted_iupac_name != "No preferred IUPAC name found" else "No preferred name found"

            data = {
                # Function to display the Molecule's Name on the Molecular ANalysis.html
                "home_mol_name" : home_name,
                # Calling function for accessing the sub-functions
                # 1. Molecular Property
                "molecular_property" : molecular_property(mol, compound_name, predicted_iupac_name),
                # 2. Lipinski Property of 5
                "lipinski5": lipinski_5(mol),
                #3. Molecular Fingerprint
                "morganfp": morgan_fp(mol),
                #4. Atom Analysis
                "atom_analysis": atom_analysis(mol),
                #5. Bond Analysis
                "bond_analysis": bond_analysis(mol),
                # 6. Substruct Search
                "substruct_search": substruct_result,
                # "delete_fragment": delete_substruct_fragment,
                # "replace_fragment": replace_substruct_fragment, 
                # 7. identify Groups
                "identify_groups": identify_groups(mol),
                # 8. Molecular Complexity
                "molecular_complexity": molecular_complexity(mol),
                # 9. 2D Structure - SVG Image
                "svg": render_molecule_svg(smiles),
                # 10. 3D Structure
                "html_3d": generate_3d_html(smiles),
                #11. Molecular Descriptors
                "molecular_descriptors": molecular_descriptors(mol),
                #11. Molecular Scaffold Analysis
                "scaffold_analysis": scaffold_analysis(mol)
            }

            cache.set(smiles, data) # Caching 'data' so that other routes can use it
            
            # Save JSON cache for exports
            json_filename = f"cache_{hash(smiles)}.json"
            json_filename
            with open(json_filename, "w", encoding="utf-8") as f:
                json.dump(data, f)

            # Store active SMILES & cache file in session
            session["cached_json_file"] = json_filename
            session["active_smiles"] = smiles

            # 2. REDIRECT to GET request (Clears the POST payload!)
            # return redirect(url_for("molecule_analysis"))
            
    else:
        print("Smiles Input Error: Please submit a valid SMILES structure above to query properties.")

    return render_template("molecule_analysis.html", data=data, error=error)


@app.route("/receive_smiles_for_home", methods=["POST"])
# route to send the smile data to the terminal from ketcher
def receive_smiles():
    data = request.get_json()
    # print(data)
    print(f"Smiles from ketcher for Molecule Analysis is {data["smiles"]}")

    return {"status": "ok"}


def molecular_property(mol, compound_name, predicted_iupac_name):
    # data = request.get_json()
    # smiles = data["smiles"]
    smiles = Chem.MolToSmiles(mol)
    
    return{
        "compound_name": compound_name['pref_name'],
        "common_name": compound_name['all_synonyms'],
        # "iupac_name": get_iupac_name_by_smiles(smiles),
        "predicted_iupac_name": predicted_iupac_name,
        "predicted_chembl_id": compound_name['chemblid'],
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

def lipinski_5(mol):
    return{
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "lipinski": compute_lipinski(mol)
    }

def morgan_fp(mol):
    mf = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useChirality=True)
    morganfp = list(mf)
    morganfp_2048 = morganfp[:2049]
    onbits = len(mf.GetOnBits())
    idxonbits = tuple(mf.GetOnBits())
    return{
        "morganfingerprint": morganfp_2048,
        "lenmf": len(mf),
        "onbits": onbits,
        "idxonbits": idxonbits
    }

def atom_analysis(mol):
    atom_data_list = []
    
    for atom in mol.GetAtoms():
        # Extracted metrics formatted as clean key-value dictionary items
        atom_metrics = {
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
            "aromatic": "Yes" if atom.GetIsAromatic() else "No",
            "chiral": str(atom.GetChiralTag()),
            "isotope": atom.GetIsotope()
        }
        atom_data_list.append(atom_metrics)
        
    return atom_data_list

def bond_analysis(mol):
    bond_data_list = []
    
    for bond in mol.GetBonds():
        # Extracted metrics formatted as clean key-value dictionary items
        bond_metrics = {
            "begin_atom_index": bond.GetBeginAtomIdx(),
            "end_atom_index": bond.GetEndAtomIdx(),
            "bond": bond.GetBondTypeAsDouble(),
            "bond_type": str(bond.GetBondType()),
            "is_ring": "Yes" if bond.IsInRing() else "No"
        }
        # Never use hyphens in dictionary keys for Jinja templating.
        bond_data_list.append(bond_metrics)

    return bond_data_list
        
def substruct_search(mol, smarts_str):

    if not smarts_str:
        return {
            "present": False, 
            "indices": [], 
            "svg": None, 
            "count_matches": None,
            "error": "No SMARTS pattern entered."
            }

    patt = Chem.MolFromSmarts(smarts_str)

    if patt is None:
        return {
            "present": False, 
            "indices": [], 
            "svg": None, 
            "count_matches": None,
            "error": "Invalid SMARTS pattern syntax."
            }

    # 1. Gather matches (uniquify=False captures all symmetry-equivalent mappings)
    # matches = mol.GetSubstructMatches(patt, uniquify=False)

    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return {
            "present": False, 
            "indices": [], 
            "svg": None, 
            "count_matches": None,
            "error": "Substructure not found in this molecule."
            }

    # 2. Flatten and sort atom indices for highlighting
    highlight_atoms = sorted(set(atom for match in matches for atom in match))

    # 3. Generate the highlighted SVG image string
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)

    # Optional styling adjustments for high-contrast dark theme compatibility
    options = drawer.drawOptions()
    # options.clearBackground = False  
    # Allows transparent backgrounds

    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=highlight_atoms
    )
    drawer.FinishDrawing()

    # 4. Extract drawing text as standard UTF-8 string text
    svg_text = drawer.GetDrawingText()

    return {
        "pattern": smarts_str,
        "present": True,
        "count_matches": len(matches),
        "highlight_atoms": len(highlight_atoms),
        "indices": list(matches), # Returns tuple of tuples of matching paths
        "svg": svg_text,
        "error": None
    }


def delete_fragment(mol, smarts_str):
    if not smarts_str:
        return {
            "present": False, 
            "svg": None, 
            "error": "No SMARTS pattern entered."
            }
    
    patt = Chem.MolFromSmarts(smarts_str)
    
    if patt is None:
        return {
            "present": False, 
            "svg": None, 
            "error": "Invalid SMARTS pattern to be deleted."
            }

    new_mol = AllChem.DeleteSubstructs(mol,patt)
    new_smiles = Chem.MolToSmiles(new_mol)

    new_deleted_svg = render_molecule_svg(new_smiles)

    deleted = True if new_deleted_svg else False

    return {
        "new_smiles": new_smiles,
        "pattern": smarts_str,
        "present": True,
        "deleted": deleted,
        "svg": new_deleted_svg,
        "error": None
    }


def replace_fragment(mol, smarts_str, new_smarts_str):
    # 1. Input validations
    if not smarts_str:
        return {
            "present": False,
            "replaced": False,
            "error": "No SMARTS pattern entered to replace."
        }
        
    if not new_smarts_str:
        return {
            "present": False,
            "replaced": False,
            "error": "No NEW replacement SMARTS pattern entered."
        }

    # 2. Parse target pattern (patt1) - Try SMARTS first, then fallback to SMILES
    patt1 = Chem.MolFromSmarts(smarts_str)
    if patt1 is None:
        patt1 = Chem.MolFromSmiles(smarts_str)
        
    if patt1 is None:
        return {
            "present": False,
            "replaced": False,
            "error": "Invalid SMARTS/SMILES pattern to be replaced."
        }

    # 3. Parse replacement pattern (patt2) - Try SMARTS first, then fallback to SMILES
    patt2 = Chem.MolFromSmiles(new_smarts_str)
    if patt2 is None:
        patt2 = Chem.MolFromSmiles(new_smarts_str)
        
    if patt2 is None:
        return {
            "present": False,
            "replaced": False,
            "error": "Invalid SMARTS/SMILES Replacement Fragment."
        }

    # 4. Verify target exists inside the molecule
    if not mol.HasSubstructMatch(patt1):
        return {
            "present": False,
            "replaced": False,
            "error": "Target pattern not found in the target molecule."
        }

    try:
        # 5. Perform replacement operation on all matches
        new_mols = AllChem.ReplaceSubstructs(mol, patt1, patt2, replaceAll=True)
        if not new_mols:
            return {
                "present": False,
                "replaced": False,
                "error": "Substructure replacement failed."
            }

        new_mol = new_mols[0]
        
        # Sanitize to fix valence and implicit hydrogens after modification
        Chem.SanitizeMol(new_mol)
        new_smiles = Chem.MolToSmiles(new_mol)
        
        # 6. Generate SVG for the clean new molecule
        new_replaced_svg = render_molecule_svg(new_smiles)

        # 7. Generate Highlighted SVG for newly inserted fragment
        # Get matches of patt2 in the updated molecule
        matches = new_mol.GetSubstructMatches(patt2, uniquify=False)
        highlight_atoms = sorted(set(atom for match in matches for atom in match)) if matches else []

        drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)
        
        # Handle transparent drawing preparation
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            new_mol,
            highlightAtoms=highlight_atoms
        )
        drawer.FinishDrawing()
        highlighted_svg_text = drawer.GetDrawingText()

        return {
            "present": True,
            "replaced": True,
            "pattern": smarts_str,
            "new_pattern": new_smarts_str,
            "new_smiles": new_smiles,
            "highlighted_svg": highlighted_svg_text,
            "svg": new_replaced_svg,
            "error": None
        }

    except Exception as e:
        print(f"Error in replace_fragment: {e}")
        return {
            "present": False,
            "replaced": False,
            "error": f"Failed to modify structure: {str(e)}"
        }

@app.route("/delete_replace", methods=["POST"])
def delete_replace():
    smiles = request.form.get("smiles1", "").strip()
    smarts1 = request.form.get("smarts1", "").strip()
    newsmarts1 = request.form.get("newsmarts1", "").strip()
    # action_target = request.form.get("action_target", "").strip()

    mol = Chem.MolFromSmiles(smiles) if smiles else None

    # delete_result = {"present": False}
    # replace_result = {"present": False}

    master_data = cache.get(smiles)

    error = None

    if not master_data:
        master_data = {
            "smiles": smiles,
            "molecular_property": molecular_property(mol),
            "lipinski5": lipinski_5(mol),
            "substruct_search": {"pattern": smarts1},
            "delete_fragment": {"present": False},
            "replace_fragment": {"present": False}
        }

    # 3. Process requested fragment action
    if "deleteaction" in request.form:
        master_data["delete_fragment"] = delete_fragment(mol, smarts1)
        
    elif "replaceaction" in request.form:
        master_data["replace_fragment"] = replace_fragment(mol, smarts1, newsmarts1)

    # Preserve the active search input pattern
    # Ensure 'substruct_search' dictionary exists without wiping out previous search results
    if "substruct_search" not in master_data or not isinstance(master_data["substruct_search"], dict):
        master_data["substruct_search"] = {"pattern": smarts1, "present": False}
    else:
        # Update ONLY the pattern string, keeping 'present', 'svg', etc. intact
        master_data["substruct_search"]["pattern"] = smarts1

    # 4. Save/Update server cache with updated delete/replace results
    cache.set(smiles, master_data)

    # 5. Render template using cached data (No heavy re-computations!)
    return render_template(
        "molecule_analysis.html", 
        data=master_data, 
        error=error, 
        active_tab="substruct_search"
    )

    # return redirect(
    #     url_for("molecule_analysis", 
    #             data=master_data, 
    #             error=error, 
    #             active_tab="substruct_search"))


def identify_groups(mol):
    FUNCTIONAL_GROUPS = {
    "Alcohol": "[#6][OX2H]",
    "Aldehyde": "[CX3H1](=O)[#6]",
    "Alkene": "C=C",
    "Alkyne": "C#C",
    "Amide": "C(=O)N",
    "Amine (Primary)": "[NX3;H2;!$(NC=O)]",
    "Amine (Tertiary)": "[NX3;H0;!$(NC=O)]",
    "Aromatic Ring": "c1ccccc1",
    "Carbonyl": "C=O",
    "Carboxylic Acid": "C(=O)[OX2H1]",
    "Ether": "[OD2]([#6])[#6]",
    "Ester": "C(=O)O[#6]",
    "Halogen": "[F,Cl,Br,I]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Nitrile": "C#N",
    "Nitro": "[N+](=O)[O-]",
    "Phenol": "c[OX2H]",
    "Thiol": "[SX2H]",
    "Sulfide": "[SX2]"
}
    results = {}

    for name, smarts in FUNCTIONAL_GROUPS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            results[name] = {
                "present": False,
                "error": "Invalid SMARTS"
            }
            continue

        matches = mol.GetSubstructMatches(pattern)

        results[name] = {
            "present": len(matches) > 0,
            "count": len(matches),
            # "indices": matches
        }    
        
        # sno = [i for i in range(len(FUNCTIONAL_GROUPS))]
        # present = [group_data['present'] for group_name, group_data in results.items()]
        # count = [group_data['count'] for group_name, group_data in results.items()]
        # indices = [group_data['indices'] for group_name, group_data in results.items()]
            
    return results
  

def molecular_complexity(mol):
    complexity = GraphDescriptors.BertzCT(mol) 
    j_value = GraphDescriptors.BalabanJ(mol) 
    AvgIpc = GraphDescriptors.AvgIpc(mol, dMat=None, forceDMat=False) 
    tpsa = Descriptors.TPSA(mol)

    return{
        "BertzCT": complexity,
        "BalabanJ": j_value,
        "AvgIpc": AvgIpc,
        "tpsa": tpsa
    }


def render_molecule_svg(smiles):
    """Parses a SMILES string and returns a cleanly structured inline SVG string."""
    # Using this function in Molecualr Ananysis tab (2D structure card)
    # Also using it on the Home page
    # Parses a SMILES string and returns a cleanly sliced inline SVG string.
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        raw_svg = drawer.GetDrawingText()
        svg_start = raw_svg.find("<svg")
        return raw_svg[svg_start:] if svg_start != -1 else raw_svg
    except Exception:
        return ""

def generate_3d_html(smiles):
    """Generates a 3D Ball and Stick py3Dmol HTML string."""
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return None

        mol = Chem.AddHs(mol)
        if mol.GetNumAtoms() == 0:
            return None

        embed_success = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        
        if embed_success != 0:
            print("Using Random Coords to generate #d conformeer")
            embed_success = AllChem.EmbedMolecule(
                mol, AllChem.ETKDGv3(),
                useRandomCoords=True,
                )

        if embed_success != 0:
            return None

        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)

        # Does MMFF have the necessary parameters for every part of this molecule?"

        # Merck Molecular Force Field   |	Universal Force Field
        # Organic/drug-like molecules	|   Broad range of elements
        # Generally better Accuracy for typical organic molecules  |	Usually less accurate
        # More limited coverage |   Much broader coverage
        
        mol_block = Chem.MolToMolBlock(mol)

        viewer = py3Dmol.view(width=500, height=300)
        viewer.addModel(mol_block, 'mol')
        viewer.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        viewer.zoomTo()
        
        return viewer._make_html()
    except Exception as e:
        print(f"Error generating 3D molecule: {e}")
        return None


def molecular_descriptors3d(mol):
    if not mol:
        return None

    try:
        # Adding explicit hydrogens (Crucial for accurate 3D geometry)
        mol = Chem.AddHs(mol)
        # Generated the 3D conformer
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        # Optimizing the geometry
        AllChem.MMFFOptimizeMolecule(mol)

        if mol is None or mol.GetNumAtoms() == 0:
            return None
        
        descriptors_dict = Descriptors3D.CalcMolDescriptors3D(mol)
        #3D Descriptors
        return{
            "PMI1": descriptors_dict["PMI1"], 
            "PMI2": descriptors_dict["PMI2"], 
            "PMI3": descriptors_dict["PMI3"], 
            "NPR1": descriptors_dict["NPR1"], 
            "NPR2": descriptors_dict["NPR2"], 
            "RadiusOfGyration": descriptors_dict["RadiusOfGyration"], 
            "InertialShapeFactor": descriptors_dict["InertialShapeFactor"], 
            "Eccentricity": descriptors_dict["Eccentricity"], 
            "Asphericity": descriptors_dict["Asphericity"],
            "SpherocityIndex": descriptors_dict["SpherocityIndex"] , 
            "PBF": descriptors_dict["PBF"] 
        }

    except Exception as e:
        print(f"Error generating 3D molecule: {e}")
        return None


def molecular_descriptors(mol):
    try:
        complexity = GraphDescriptors.BertzCT(mol) 
        j_value = GraphDescriptors.BalabanJ(mol) 
        b_cut = "NA"
        b_cut = rdMolDescriptors.BCUT2D(mol)
        if rdMolDescriptors.BCUT2D(mol) != "NA":
            b_cut = rdMolDescriptors.BCUT2D(mol)
        else:
            b_cut = "NA"

    except:
        ''

    return{
        # 0D Descriptors
        "MolecularWeight": round(Descriptors.MolWt(mol), 4),
        "LogP": round(Descriptors.MolLogP(mol), 4),
        "HBondDonors": Descriptors.NumHDonors(mol),
        "HBondAcceptors": Descriptors.NumHAcceptors(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
        "ValenceElectrons": Descriptors.NumValenceElectrons(mol),
        "RadialEelectrons": Descriptors.NumRadicalElectrons(mol),
        "AromaticRings": Descriptors.NumAromaticRings(mol),

        # 1D Descriptors
        "MolarRefractivity": round(Crippen.MolMR(mol), 4),
        "AromaticRings": Descriptors.NumAromaticRings(mol),
        # "HBondDonors": Descriptors.NumHDonors(mol),
        # "HBondAcceptors": Descriptors.NumHAcceptors(mol),
        # "RotatableBonds": Descriptors.NumRotatableBonds(mol),
        # "TPSA": Descriptors.TPSA(mol),
        # "Heavy AtomCount": Descriptors.HeavyAtomCount(mol),
        # "ValenceElectrons": Descriptors.NumValenceElectrons(mol),
        
        # 2D Descriptors
        "BalabanJIndex": j_value,
        "BertzCT": complexity,
        "BCUT2D": b_cut,
        "Chi0v": Descriptors.Chi0v(mol),
        "Chi1v": Descriptors.Chi1v(mol),
        "Chi2v": Descriptors.Chi2v(mol),
        "Chi3v": Descriptors.Chi3v(mol),
        "Chi4v": Descriptors.Chi4v(mol),
        "Chi0n": Descriptors.Chi0n(mol),
        "Chi1n": Descriptors.Chi1n(mol),
        "Chi2n": Descriptors.Chi2n(mol),
        "Chi3n": Descriptors.Chi3n(mol),
        "Chi4n": Descriptors.Chi4n(mol),
        "Kappa1": Descriptors.Kappa1(mol),
        "Kappa2": Descriptors.Kappa2(mol),
        "Kappa3": Descriptors.Kappa3(mol),
        
        # For 3d descriptors, I created a separate function for better calculation.
        "molecular_descriptors3d": molecular_descriptors3d(mol)
    }

def scaffold_analysis(mol):
    try:
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)

        if scaffold_mol is None or scaffold_mol.GetNumAtoms() == 0:
            return {
                "scaffold_mol": None,
                "generic_scaffold": None,
                "scaffold_smiles": None,
                "generic_scaffold_smiles": None
            }
        if scaffold_mol.GetNumAtoms() == 0:
            return {
                    "scaffold_mol": "No atoms in mol",
                    "generic_scaffold": "No atoms in mol",
                    "scaffold_smiles": None,
                    "generic_scaffold_smiles": None
                }

        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        generic_scaffold_mol = MurckoScaffold.MakeScaffoldGeneric(scaffold_mol)

        return {
            "scaffold_smiles": Chem.MolToSmiles(scaffold_mol),
            "generic_scaffold_smiles": Chem.MolToSmiles(generic_scaffold_mol)
        }

    except Exception:
        return {
            "scaffold_smiles": None,
            "generic_scaffold_smiles": None
        }

# Download Route on Molecule Analysis Page
@app.route("/download", methods=["POST"])
def handle_download():
    smiles = request.form.get("smiles1", "").strip()
    action = request.form.get("action")
    mol = Chem.MolFromSmiles(smiles)
    try: 
        # Fast JSON file download
        if action == "data_json":
            json_file = session.get("cached_json_file")
            if json_file and os.path.exists(json_file):
                download_name="molecule_analysis.json"
                print(f"Molecule data saved as {download_name}")
                return send_file(
                    json_file, as_attachment=True, download_name = download_name)
            return "Cached JSON data missing. Please re-submit SMILES.", 400
        
        # Fast CSV download using cached JSON data
        if action == "data_csv":
            json_file = session.get("cached_json_file")

            if json_file and os.path.exists(json_file):
                csv_path = download_data_csv(json_file)
                return send_file(
                    csv_path,
                    as_attachment=True,
                    download_name="molecule_analysis.csv",
                )
            return "Cache missing. Please re-submit SMILES.", 400

        if action == "data_sdf":
            json_file = session.get("cached_json_file")

            if json_file and os.path.exists(json_file):
                # Generate the in-memory BytesIO buffer
                sdf_buffer = download_data_sdf(json_file)
                
                # Send directly to user's browser as a file download
                return send_file(
                    sdf_buffer,
                    as_attachment=True,
                    download_name="molecule_analysis.sdf",
                    mimetype="chemical/x-mdl-sdfile"
                )
            return "Cache missing. Please re-submit SMILES.", 400

        if action == "data_mol":
            json_file = session.get("cached_json_file")

            if json_file and os.path.exists(json_file):
                # Generate the in-memory BytesIO buffer
                sdf_buffer = download_data_mol(json_file)
                
                # Send directly to user's browser as a file download
                return send_file(
                    sdf_buffer,
                    as_attachment=True,
                    download_name="molecule_analysis.mol",
                    mimetype="chemical/x-mdl-molfile"
                )
            return "Cache missing. Please re-submit SMILES.", 400
        
    except Exception as e:
        print(e)
        return str(e), 500

    # Non-CSV actions still use RDKit mol object directly
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    action_handlers = {
        'data_sdf': download_data_sdf,
        'data_mol': download_data_mol,
        # 'data_mol2': download_data_mol2,
        'data_csv': download_data_csv,

        "2d_sdf": save_as_2d_sdf,
        "2d_mol": save_as_2d_mol,
        "2d_mol2": save_as_2d_mol2,
        "2d_pdf": download_pdf,
        "2d_svg": download_svg,
        "2d_png": download_png,

        "3d_sdf": save_as_3d_sdf,
        "3d_mol": save_as_3d_mol,
        "3d_mol2": save_as_3d_mol2,
        "3d_html": save_as_html,
        "3d_json": save_as_json,
    }

    handler = action_handlers.get(action)
    if handler and mol:
        return handler(mol)

    if handler == "download_data_sdf" and mol:
        return handler(mol)

    return "Invalid action or structure", 400

def save_smiles_to_sdf(mol, output_file):

    output_file="smiles.sdf"

    smiles = smiles

    mol.SetProp("_Name","mol-2d")
    AllChem.Compute2DCoords(mol)

    mol.SetProp("_Name","mol-3d")
    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    if status != 0:
        raise RuntimeError("3D coordinate generation failed.")

    # Optimize geometry
    AllChem.MMFFOptimizeMolecule(mol)

    # Calculate molecular properties
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mol_weight = Descriptors.MolWt(mol)
    exact_mass = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    num_atoms = mol.GetNumAtoms()
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    # Store properties in SDF
    mol.SetProp("_Name", "mol-2d mol-3d")
    mol.SetProp("SMILES", smiles)
    mol.SetProp("Formula", formula)
    mol.SetProp("MolecularWeight", f"{mol_weight:.4f}")
    mol.SetProp("ExactMass", f"{exact_mass:.4f}")
    mol.SetProp("LogP", f"{logp:.4f}")
    mol.SetProp("TPSA", f"{tpsa:.4f}")
    mol.SetProp("NumAtoms", str(num_atoms))
    mol.SetProp("NumHeavyAtoms", str(num_heavy_atoms))
    mol.SetProp("Charge", str(Chem.GetFormalCharge(mol)))

    # Write SDF
    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()

    print(f"SDF written to {output_file}")


def mol_summary(mol):
    return 'Coming soon..'

# @app.route("/ketcher")
# def ketcher():
#     return render_template("ketcher.html")

@app.route("/contact", methods=["GET", "POST"])
def contact():
    return render_template("contact.html")


@app.route("/about", methods=["GET", "POST"])
def about():
    return render_template("about.html")



# < - - - - - - -- - - - - - - - - - - - - - Similarity NavBar Route - - - - - - - - - - - - - - - - - - - - - - - >
#
# 

@app.route("/similarity", methods=["GET", "POST"])
def similarity():
    return render_template("similarity.html")

def render_similarmolecule_svg(smiles):
    """Parses a SMILES string and returns a cleanly structured inline SVG string."""
    # Using this function in Molecualr Ananysis tab (2D structure card)
    # Also using it on the Home page
    # Parses a SMILES string and returns a cleanly sliced inline SVG string.
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        raw_svg = drawer.GetDrawingText()
        svg_start = raw_svg.find("<svg")
        return raw_svg[svg_start:] if svg_start != -1 else raw_svg
    except Exception:
        return ""

@app.route('/search', methods=['GET', 'POST'])
def similarity_search():
    data = None
    error = None

    # Use request.form for POST requests, fallback to request.args for GET if needed
    raw_smiles_chemblid = (
        request.form.get("chembl_id1")
        if request.method == "POST"
        else request.args.get("chembl_id1", "")
        )

    # raw_smiles_chemblid = request.args.get('chembl_id1', '').strip()
    clean_str = raw_smiles_chemblid.strip()
    print("Input for similarity search engine:", clean_str)

    check_str = clean_str.upper()
    # 1. Check for ChEMBL ID
    if clean_str[6:].isdigit() or check_str.startswith("CHEMBL"): 
        try:
            clean_str = clean_str.upper()
            clean_str = "".join(clean_str.split()) # allowing chembl 1000 to be accepted despite having a space in between
            print(f"Similarity Search Engine get CHEMBL ID:", clean_str)
            clean_str = clean_str.upper()
            # if clean_str.startswith("CHEMBL") and clean_str[6:].isdigit(): 
                # format of valid chembl id
            data= similarity_search_by_chemid(clean_str)
            return render_template("similarity.html", data=data)
        
        except:
            # error = "This CHEMBL ID does not exists or Invalid CHEMBL ID"
            print(f"Similarity Search Engine received invalid input:", clean_str)
            return render_template("similarity.html", data=data, error=error)
    
    # 2. Check for SMILES String
    # Common structural syntax characters and lowercase aromatic elements
    smiles_indicators = {'=', '#', '@', '(', ')', '[', ']', 'c', 'n', 'o', 's'}
    smiles_special = {'=', '#', '@', '(', ')', '[', ']', '.', '/', '\\', '%', '+' , '-'}
    organic_atoms = {'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I', 'Cl', 'Br', 'b', 'c', 'n', 'o', 'p', 's'}
    elements = {
        'H', 'He',
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
        'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
        'Es', 'Fm', 'Md', 'No', 'Lr',
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
        'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
        }
    
    # If it contains SMILES symbols smile related symbols or organic atoms or 
    # any element of periodic table or contains numbers (for rings)

    clean_str = raw_smiles_chemblid.strip()

    has_smiles_chars = any(char in smiles_indicators for char in clean_str)
    has_smiles_special_chars = any(char in smiles_special for char in clean_str)
    has_organic_smiles_chars = any(char in organic_atoms for char in clean_str)
    has_elements = any(char in elements for char in clean_str)
    has_digits = any(char.isdigit() for char in clean_str)

    
    if has_smiles_chars or has_smiles_special_chars or has_organic_smiles_chars or has_elements or has_digits:
        try:
            print(f"Similarity Search Engine get SMILES: {clean_str}")
            data = similarity_search_by_smiles(clean_str)
            return render_template("similarity.html", data=data)
        
        except Exception as e:
            print(f"Search execution failed: {str(e)}")
            print(f"Similarity Search Engine received invalid input:", {clean_str})
            error = "Invalid Input"
            return render_template("similarity.html", data=data, error=error)

    # error = f"{clean_str} is invalid input. Please enter a valid CHEMBL ID or SMILES string."

    return render_template("similarity.html", data=data, error="Unrecognized input format. Please enter a valid CHEMBL ID or SMILES string.")


def similarity_search_by_smiles(clean_str):
    """
    Looks up a target molecule directly by its SMILES string, standardizes it, 
    and runs a high-speed Tanimoto coefficient matrix similarity query.
    """
    error = None
    data = None

    default_threshold = 0.90  # Default similarity threshold if not provided
    new_threshold = request.form.get("myRange") if request.method == "POST" else default_threshold
    threshold = float(new_threshold) if new_threshold else default_threshold
    
    # 1. Fetch raw SMILES input string from the URL parameters

    # raw_smiles_chemblid = request.args.get('chembl_id1', '').strip()
    # clean_str = raw_smiles_chemblid.strip()
    raw_smiles = clean_str

    if not raw_smiles:
        error = "No SMILES string provided"
        return error
    # return jsonify({"error": "No SMILES string provided"}), 400

    # 2. Canonicalize the input SMILES to ensure consistency with the database
    try:
        mol = Chem.MolFromSmiles(raw_smiles)
        if mol is None:
            error = "Invalid SMILES string structure"
            # Failed to parse invalid SMILES
            print(error)
            data = {
                "target_query": None,
                "target_smiles": raw_smiles,
                "total_matches_found": None,
                "error": error
            }
            return data
            # return jsonify({"error": "Invalid SMILES string structure format."}), 400
        target_smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
      print(f"Failed to parse SMILES: {str(e)}")

    # 3. Find the molecule's integer index position inside df_metadata if it exists
    # (Used strictly to handle the self-comparison filter step below)
    target_row = df_metadata[df_metadata['canonical_smiles'] == target_smiles]
    # storing the index of the target molecule if it exists in the metadata dataframe, otherwise None
    target_int_id = int(target_row.index[0]) if not target_row.empty else None
    
    # If it's a completely novel molecule not in the dataframe metadata, 
    # Assign a placeholder query name
    target_id = str(target_row.iloc[0]['chembl_id']) if not target_row.empty else "NOVEL_QUERY"

    try:
        # 4. Query FPSim2 dynamically using half of the available logical threads
        optimal_workers = max(1, os.cpu_count() // 2) 

        # Positional Order: (query_smiles, threshold, n_workers)
        # results = fpe.similarity(target_smiles, 0.90, 'tanimoto', optimal_workers)
        results = fpe.similarity(target_smiles, threshold, 'tanimoto', optimal_workers)

        matches = []
        for match in results:
            match_int_id = int(match[0])     # Matched row index inside .h5 matrix
            score = round(float(match[1]), 6) # Calculated Tanimoto Similarity score
            
            # 5. Skip self-comparison logic using row markers
            if target_int_id is not None and match_int_id == target_int_id:
                continue
                
            # 6. Map the matching integer pointer back to its native metadata values
            real_chembl_id = df_metadata.iloc[match_int_id]['chembl_id']
            smiles = df_metadata.iloc[match_int_id]['canonical_smiles']
            # chembl_name_of_similar_mol = get_pref_name_from_smiles_locally(smiles) if new_threshold >= '0.80' else ''

            matches.append({
                "chembl_id": real_chembl_id,
                "Tanimoto_Similarity": score,
                "canonical_smiles": smiles,
                # "chembl_name": get_pref_name_from_smiles_locally(smiles),
                # "chembl_name": chembl_name_of_similar_mol,
                "svg_of_similar": render_similarmolecule_svg(smiles)
            })
            
        # 7. Order results from highest structural affinity to lowest
        matches = sorted(matches, key=lambda x: x['Tanimoto_Similarity'], reverse=True)

        data = {
            "target_query": target_id,
            "target_smiles": target_smiles,
            "target_chembl_name": get_pref_name_from_smiles_locally(target_smiles),
            "total_matches_found": len(matches),
            "results": matches,
            "svg": render_molecule_svg(raw_smiles), # svg of target similar smiles
            "error": error,
            "threshold": threshold
        }
        return data

    except Exception as e:
        print(f"Search execution failed: {str(e)}")
        return None


def similarity_search_by_chemid(clean_str):
    """
    Looks up a target molecule by ChEMBL ID, maps it to its internal SMILES string, and runs a high-speed Tanimoto coefficient matrix similarity query.
    """
    data = None
    error = None

    default_threshold = 0.90  # Default similarity threshold if not provided
    new_threshold = request.form.get("myRange") if request.method == "POST" else default_threshold
    threshold = float(new_threshold) if new_threshold else default_threshold
    
    target_id = clean_str

    if not target_id:
        print(f"No ChEMBL ID provided")

    # 1. Fetch the metadata row matching the user's requested ChEMBL ID string
    target_row = df_metadata[df_metadata['chembl_id'] == target_id]
    
    if target_row.empty:
        if clean_str[6:].isdigit():
            print(f"ID '{target_id}' not found in the dataset database")
            error_tuple = (target_id, "does not exist")
            # error_tuple = ('CHEMBL1000005', 'does not exist or is invalid.')
            # Join the items together with a space
            error = f"{error_tuple[0]} {error_tuple[1]}"
            data = {
                    "target_query": target_id,
                    "target_smiles": None,
                    "total_matches_found": None,
                    "error": error
                    }
            return data
        print(f"ID '{target_id}' is invalid CHEMBL ID")
        error_tuple = (target_id, "is invalid CHEMBL ID")
        error = f"{error_tuple[0]} {error_tuple[1]}"
        data = {
                "target_query": target_id,
                "target_smiles": None,
                "total_matches_found": None,
                "error": error
                }
        return data
        
    # 2. Extract its reference structural index position and canonical SMILES string
    target_int_id = int(target_row.index[0])
    target_smiles = str(target_row.iloc[0]['canonical_smiles'])

    # 2. Canonicalize the input SMILES to ensure consistency with the database
    try:
        mol = Chem.MolFromSmiles(target_smiles)
        if mol is None:
            error = "Invalid SMILES string structure"
            # Failed to parse invalid SMILES
            print(error)
            data = {
                "target_query": None,
                "target_smiles": target_smiles,
                "total_matches_found": None,
                "error": error
            }
            return data
            # return jsonify({"error": "Invalid SMILES string structure format."}), 400
        target_smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        print(f"Failed to parse SMILES: {str(e)}")

    try:
        # 3. Query FPSim2 via positional strings to isolate version discrepancies
        
        # Dynamically uses half of the available logical threads to keep the app responsive
        optimal_workers = max(1, os.cpu_count() // 2) 

        # Positional Order: (query_smiles, threshold, n_workers)
        # results = fpe.similarity(target_smiles, 0.90, 'tanimoto', optimal_workers)
        results = fpe.similarity(target_smiles, threshold, 'tanimoto', optimal_workers)

        # results is a list of tuples: [(row_index, similarity_score), ...]

        matches = []
        for match in results:
            match_int_id = int(match[0])     # Matched positional row index inside .h5 matrix
            score = round(float(match[1]), 6) # Calculated Tanimoto Similarity score
            
            # Skip self-comparison logic using row markers
            if match_int_id == target_int_id:
                continue
                
            # 4. Map the matching integer pointer back to its native metadata values
            real_chembl_id = "No name found"
            real_chembl_id = df_metadata.iloc[match_int_id]['chembl_id'] 
            smiles = df_metadata.iloc[match_int_id]['canonical_smiles']
            # chembl_name_of_similar_mol = get_pref_name_from_smiles_locally(smiles) if new_threshold >= '0.80' else ''
            
            matches.append({
                "chembl_id": real_chembl_id if real_chembl_id != "No name found" else "Does not exists in CHEMBL database",
                "Tanimoto_Similarity": score,
                "canonical_smiles": smiles,
                # "chembl_name": get_pref_name_from_smiles_locally(smiles),
                # "chembl_name": chembl_name_of_similar_mol,
                "svg_of_similar": render_similarmolecule_svg(smiles)
            })
            
        # 5. Order results from highest structural affinity to lowest
        matches = sorted(matches, key=lambda x: x['Tanimoto_Similarity'], reverse=True)

        # print(results)
        data = {
            "target_query": target_id,
            "target_smiles": target_smiles,
            "target_chembl_name": get_pref_name_from_smiles_locally(target_smiles),
            "total_matches_found": len(matches),
            "results": matches,
            "svg": render_molecule_svg(target_smiles), # svg of target smiles
            "error": error,
            "threshold": threshold
        }

        return data

    except Exception as e:
       print(f"Search execution failed: {str(e)}")
       return None
   


# < - - - - - - -- - - - - - - - - - - - - - Compare NavBar Route - - - - - - - - - - - - - - - - - - - - - - - >
#
# 

@app.route("/compare", methods=["GET", "POST"])
def compare():
    return render_template("compare.html")

@app.route('/comparison', methods=["GET", "POST"])
def comparison():
    data = None
    error = None

    smiles1 = (
            request.form.get("chembl_id1")
            if request.method == "POST"
            else request.args.get("chembl_id1", "")
            ).strip()

    smiles2 = (
            request.form.get("chembl_id2")
            if request.method == "POST"
            else request.args.get("chembl_id2", "")
            ).strip()

    print("Input for comparison engine:", smiles1, "and", smiles2)

    if not smiles1 or not smiles2:
        error = "Both SMILES strings are required for comparison."
        data = {
            "smiles1": smiles1,
            "smiles2": smiles2,
            "data1": None,
            "data2": None,
            "error": error
        }
        return render_template("compare.html", data=data, error=error)

    if smiles1 == smiles2:
        error = "Both SMILES strings are identical. \nPlease provide two different SMILES strings for comparison."
        data = {
            "smiles1": smiles1,
            "smiles2": smiles2,
            "data1": None,
            "data2": None,
            "error": error
        }
        return render_template("compare.html", data=data, error=error)

    if smiles1 and smiles2:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 is None:
            error = f"Invalid SMILES string for molecule 1. Please check the input."
            data = {
                "smiles1": smiles1,
                "smiles2": smiles2,
                "data1": None,
                "data2": None,
                "error": error
            }
            return render_template("compare.html", data=data, error=error)
        if mol2 is None:
            error = f"Invalid SMILES string for molecule 2. Please check the input."
            data = {
                "smiles1": smiles1,
                "smiles2": smiles2,
                "data1": None,
                "data2": None,
                "error": error
            }
            return render_template("compare.html", data=data, error=error)

        if mol2 and mol2 is None:
            error = f"Invalid SMILES string for molecule 1 and 2."
            data = {
                "smiles1": smiles1,
                "smiles2": smiles2,
                "data1": None,
                "data2": None,
                "error": error
            }
            return render_template("compare.html", data=data, error=error)

    if mol1 and mol2:
        try:
            data = compare_molecules(mol1, mol2)
            return render_template("compare.html", data=data, error=error)
        except Exception as e: 
            print(f"Comparison execution failed: {str(e)}")
            error = "An error occurred during the comparison. \nPlease check the input SMILES strings."
            data = {
                    "smiles1": smiles1,
                    "smiles2": smiles2,
                    "data1": None,
                    "data2": None,
                    "error": error
                }
            return render_template("compare.html", data=data, error=error)

    return render_template("compare.html", data=data, error=error)


def molecular_property_for_compare(smiles):
    compound_details = get_all_details_from_smiles_locally(smiles)
    mol = Chem.MolFromSmiles(smiles)
    ring_info_count = mol.GetRingInfo().NumRings()
    all_synonyms = compound_details['all_synonyms']
    compound_details['all_synonyms'] = [name.strip() for name in all_synonyms.split(",") if name.strip()]
    return {
        "compound_name": compound_details['pref_name'],
        "common_names": compound_details['all_synonyms'],
        "predicted_iupac_name": predict_iupac_name_by_smiles(smiles),
        "predicted_chembl_id": compound_details['chemblid'],
        "svg": render_molecule_svg(smiles),
        "elemental_formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "HBondDonors": Descriptors.NumHDonors(mol),
        "HBondAcceptors": Descriptors.NumHAcceptors(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "ring_count": ring_info_count,
        "lipinski": compute_lipinski(mol),
    }

def render_molecule_svg_for_compare(smiles):
    """Parses a SMILES string and returns a cleanly structured inline SVG string."""
    # Using this function in Molecualr Ananysis tab (2D structure card)
    # Also using it on the Home page
    # Parses a SMILES string and returns a cleanly sliced inline SVG string.
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        raw_svg = drawer.GetDrawingText()
        svg_start = raw_svg.find("<svg")
        return raw_svg[svg_start:] if svg_start != -1 else raw_svg
    except Exception:
        return ""
    
def morgan_fp_for_compare(mol):
    mf = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048, useChirality=True)
    onbits = len(mf.GetOnBits())
    idxonbits = tuple(mf.GetOnBits())
    return{
        "lenmf": len(mf),
        "onbits": onbits,
        "idxonbits": idxonbits
    }

def substruct_compare(molecule1, molecule2):
    matches = molecule1.GetSubstructMatches(molecule2, useChirality=True)

    if not matches:
        return {
            "present": False, 
            "indices": None, 
            "svg": None, 
            "count_matches": None,
            "error": "Substructure not found in this molecule."
            }
    
    if matches:
        return {
                "pattern": molecule2,
                "present": True,
                "count_matches": len(matches),
                "indices": str(matches), # Returns string of tuples of matching paths
                "error": None
            }


def compare_molecules(mol1, mol2):
    smiles1 = Chem.MolToSmiles(mol1, canonical=True)
    smiles2 = Chem.MolToSmiles(mol2, canonical=True)
    # rmsd_value = rmsd_same_molecule(smiles1, smiles2) or rmsd_mcs(smiles1, smiles2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048, useChirality=True)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048, useChirality=True)
    tanimoto_similarity = round(DataStructs.TanimotoSimilarity(fp1, fp2), 4)
    dice_similarity = round(DataStructs.DiceSimilarity(fp1, fp2), 4)
    cosine_similarity = round(DataStructs.CosineSimilarity(fp1, fp2), 4)
    similarity = [tanimoto_similarity, dice_similarity, cosine_similarity]
    avg_sim = sum(similarity)/3
    avg_similarity = round(avg_sim, 4)
    
    data1 = {
        "smiles1": smiles1,
        "svg1": render_molecule_svg_for_compare(smiles1),
        "molecular_property1": molecular_property_for_compare(smiles1),
        "aromatic_ring1": Descriptors.NumAromaticRings(mol1),
        "morgan_fp1": morgan_fp_for_compare(mol1),
        "molecular_complexity1": molecular_complexity(mol1),
        "tanimoto_similarity": tanimoto_similarity,
        "dice_similarity": dice_similarity,
        "cosine_similarity": cosine_similarity,
        "average_similarity": avg_similarity,
        "substructure1": substruct_compare(mol1, mol2),
        # "rmsd": rmsd_value
    }

    data2 = {
        "smiles2": smiles2,
        "svg2": render_molecule_svg_for_compare(smiles2),
        "molecular_property2": molecular_property_for_compare(smiles2),
        "morgan_fp2": morgan_fp_for_compare(mol2),
        "aromatic_ring2": Descriptors.NumAromaticRings(mol2),
        "molecular_complexity2": molecular_complexity(mol2),
        "substructure2": substruct_compare(mol2, mol1),
    }

    return{
        "data1": data1,
        "data2": data2
    }





if __name__ == "__main__":
    app.run(host="0.0.0.0", port=9000) 

# flask run --host=127.0.0.1 --port=8000









