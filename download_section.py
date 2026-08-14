#  Download Section #


# Download Data
#     type="submit" name="action" value="data_sdf">SDF
#     type="submit" name="action" value="data_mol">MOL
#     type="submit" name="action" value="data_json">JSON>
#     type="submit" name="action" value="data_csv">CSV<

#     2D Structure
#     type="submit" name="action" value="2d_sdf">SDF
#     type="submit" name="action" value="2d_mol">MOL
#     type="submit" name="action" value="2d_mol2">MOL2
#     type="submit" name="action" value="2d_pdf">PDF
#     type="submit" name="action" value="2d_svg">SVG
#     type="submit" name="action" value="2d_png">PNG

#     3D Structure
#     type="submit" name="action" value="3d_sdf">SDF
#     type="submit" name="action" value="3d_mol">MOL
#     type="submit" name="action" value="3d_mol2">MOL
#     type="submit" name="action" value="3d_html">HTML
#     type="submit" name="action" value="3d_json">JSON
# </div>

from flask import send_file
from openbabel import openbabel
import cairosvg
import csv
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors3D, Draw
import py3Dmol
import os 
import io
import json
from rdkit.Chem.Draw import rdMolDraw2D
from flask import request


#################################### Download Data ####################################

# Download Data (SDF)
def download_data_sdf(json_filepath, filename="molecule_analysis.sdf"):
    try:
        with open(json_filepath, "r+t", encoding="utf-8") as f:
            data = json.load(f)

        smiles = data["molecular_property"]["smiles"]

        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)

        # Add Hydrogens & Generate 3D Coordinates
        mol = Chem.AddHs(mol)
        status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if status == 0:
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            # Fallback to 2D coordinates if 3D embedding fails
            AllChem.Compute2DCoords(mol)

        # Name
        mol.SetProp("_Name", data["home_mol_name"])

        # Set Molecular Properties
        for key, value in data.get("molecular_property", {}).items():
            if key != "svg":
                mol.SetProp(key, str(value))

        # Set Lipinski Properties
        for key, value in data.get("lipinski5", {}).items():
            mol.SetProp(f"Lipinski_{key}", str(value))

        # Set Molecular Descriptors
        for key, value in data.get("molecular_descriptors", {}).items():
            mol.SetProp(str(key), str(value))

        # Write SDF to an in-memory Text Stream (StringIO)
        sio = io.StringIO()
        writer = Chem.SDWriter(sio)
        writer.write(mol)
        writer.close()

        # Get raw SDF string and convert to BytesIO for Flask send_file
        sdf_text = sio.getvalue()
        bio = io.BytesIO(sdf_text.encode("utf-8"))

        print(f"Molecule data saved as {filename}")

        return bio
        
    except Exception as e:
        print(e)
        return str(e), 500

    
# Download Data (MOL)
def download_data_mol(json_filepath, filename="molecule_analysis.mol"):
    try:
        """
        Reads JSON data and generates a .MOL file content in memory.
        """
        with open(json_filepath, "r", encoding="utf-8") as f:
            data = json.load(f)

        smiles = data["molecular_property"]["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES pattern found in JSON.")

        # Add Hydrogens & Generate 3D Coordinates
        mol = Chem.AddHs(mol)
        status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if status == 0:
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            # Fallback to 2D coordinates if 3D embedding fails
            AllChem.Compute2DCoords(mol)

        # Set Molecule Name (header of MOL block)
        mol.SetProp("_Name", data.get("home_mol_name", "Molecule"))

        # Convert Mol object directly to a MOL block string
        mol_block = Chem.MolToMolBlock(mol)

        # Convert string to BytesIO buffer for Flask send_file
        bio = io.BytesIO(mol_block.encode("utf-8"))

        print(f"Molecule data saved as {filename}")

        return bio

    except Exception as e:
        print(e)
        return str(e), 500

    
# Download Data (MOL2)
def download_data_mol2(mol, filename="molecule_analysis.mol2"):
    if mol is None:
        return False
    try:
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("sdf", "mol2")
        mol = openbabel.OBMol()
        sdf_file = save_as_2d_sdf(mol)
        sdf_file
        conv.ReadFile(mol, "molecule_analysis.mol")
        conv.WriteFile(mol, "molecule_analysis.mol2")
        print(f"Molecule data saved as {filename}")
        return send_file( filename, as_attachment=True,download_name=os.path.basename(filename)
            )
        
    except Exception as e:
        print(e)
        return str(e), 500


# Download Data (JSON)
# def download_data_json(json_filepath, export_json="molecule_analysis_export.json"
# ):
#     """Loads cached calculation data from JSON for download."""
#     with open(json_filepath, "r", encoding="utf-8") as f:
#         # data = json.load(f)

#         print(f"Molecule data saved as {export_json}")
#     return export_json


# Download Data (CSV)
def download_data_csv(json_filepath, export_csv="molecule_analysis_export.csv"
):
    """Loads cached calculation data from JSON and writes to CSV without recomputing."""
    with open(json_filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    with open(export_csv, mode="w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        for title, section_data in data.items():
            section_title = title.replace("_", " ").upper()
            writer.writerow([f"--- {section_title} ---"])

            if isinstance(section_data, dict):
                writer.writerow(["Property", "Value"])  # Header for key-value
                for key, val in section_data.items():
                    writer.writerow([key, val])
            elif isinstance(section_data, list):
                if len(section_data) > 0 and isinstance(section_data[0], dict):
                    # List of dicts: write table headers first, then rows
                    headers = list(section_data[0].keys())
                    writer.writerow(headers)
                    for item in section_data:
                        writer.writerow([item.get(h, "") for h in headers])
                else:
                    # Simple list of values
                    for item in section_data:
                        writer.writerow([item])

            # Handle Primitive Data (Strings, Numbers, Booleans)
            else:
                writer.writerow(["Value", section_data])

            writer.writerow([])  # Spacing row

        # return filename
        print(f"Molecule data saved as {export_csv}")
        return export_csv


#################################### 2D Structure 1(SDF) ####################################

# 2D Structure 1(SDF) 
def save_as_2d_sdf(mol, filename="mol-2d.sdf"):
    if mol is None:
        return False

    try:
        mol = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_Name", "mol-2d")

        writer = Chem.SDWriter(filename)
        writer.write(mol)
        writer.close()

        print(f"Molecule data saved as {filename}")

        return send_file(
        filename, as_attachment=True, download_name=os.path.basename(filename)
        )
    
    except Exception as e:
        print(e)
        return str(e), 500

# 2D Structure 2(MOL)
def save_as_2d_mol(mol, filename="mol-2d.mol"):
    if mol is None:
        return False

    try:
        mol = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_Name", "mol-2d")
        Chem.MolToMolFile(mol, filename)
        print(f"Molecule data saved as {filename}")
        return send_file(
        filename, as_attachment=True, download_name=os.path.basename(filename)
        )

    except Exception as e:
        print(e)
        return str(e), 500

# 2D Structure 3(MOL2)
def save_as_2d_mol2(mol, filename="mol-2d.mol2"):
  if mol is None:
        return False
  try:
      conv = openbabel.OBConversion()
      conv.SetInAndOutFormats("sdf", "mol2")
      mol = openbabel.OBMol()
      sdf_file = save_as_2d_sdf(mol)
      sdf_file
      conv.ReadFile(mol, "mol-2d.sdf")
      conv.WriteFile(mol, "mol-2d.mol2")
      print(f"Molecule data saved as {filename}")
      return send_file(
        filename, as_attachment=True,download_name=os.path.basename(filename)
        )
    
  except Exception as e:
    print(e)
    return str(e), 500


# 2D Structure 4(PDF)
def download_pdf(mol, basename="molecule"):
  if mol is None:
        return False

  try:
      mol = Chem.Mol(mol)
      AllChem.Compute2DCoords(mol)
      drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
      drawer.DrawMolecule(mol)
      drawer.FinishDrawing()
      svg = drawer.GetDrawingText()
      with open("molecule.svg", "w") as f:
          f.write(svg)

      cairosvg.svg2pdf(url="molecule.svg", write_to="molecule.pdf")
      print(f"Saved {basename}.pdf")
      return send_file(
        "molecule.pdf", as_attachment=True, download_name ="molecule.pdf"
        )

  except Exception as e:
    print(e)
    return str(e), 500


# 2D Structure 5(SVG)
def download_svg(mol, basename="molecule"):
  if mol is None:
        return False

  try:
      mol = Chem.Mol(mol)
      AllChem.Compute2DCoords(mol)
      drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
      drawer.DrawMolecule(mol)
      drawer.FinishDrawing()

      svg = drawer.GetDrawingText()

      with open("molecule.svg", "w") as f:
          f.write(svg)
          print(f"Saved {basename}.svg")
      return send_file(
        "molecule.svg", as_attachment=True, download_name =" molecule.svg"
        )

  except Exception as e:
    print(e)
    return str(e), 500


# 2D Structure 6(PNG)
def download_png(mol, basename="molecule", size=(500, 500)):
  if mol is None:
        return False

  try:
      mol = Chem.Mol(mol)
      AllChem.Compute2DCoords(mol)
      img = Draw.MolToImage(mol, size=size)
      img.save(f"{basename}.png")
      print(f"Saved {basename}.png")
      return send_file(
        "molecule.png", as_attachment=True, download_name =" molecule.png"
        )

  except Exception as e:
    print(e)
    return str(e), 500

#################################### 3D Structure ####################################

# 3D Structure 1(SDF)
def save_as_3d_sdf(mol, filename="mol-3d.sdf"):
    if mol is None:
        return False

    try:
        mol = Chem.Mol(mol)
        mol = Chem.AddHs(mol)

        mol.SetProp("_Name", "mol-3d")

        # Generate 3D coordinates
        status = AllChem.EmbedMolecule(mol,AllChem.ETKDGv3())

        if status != 0:
            raise RuntimeError("3D coordinate generation failed.")
        
        AllChem.MMFFOptimizeMolecule(mol)

        writer = Chem.SDWriter(filename)
        writer.write(mol)
        writer.close()

        print(f"SDF saved as {filename}")

        return send_file(
        filename, as_attachment=True, download_name=os.path.basename(filename)
        )
    
    except Exception as e:
        print(e)
        return str(e), 500

# 3D Structure 2(MOL)
def save_as_3d_mol(mol, filename="mol-3d.mol"):
    if mol is None:
        return False

    try:
        mol = Chem.Mol(mol)
        mol = Chem.AddHs(mol)

        mol.SetProp("_Name", "mol-3d")

        # Generate 3D coordinates
        status = AllChem.EmbedMolecule(mol,AllChem.ETKDGv3())

        if status != 0:
            raise RuntimeError("3D coordinate generation failed.")
        
        AllChem.MMFFOptimizeMolecule(mol)

        Chem.MolToMolFile(mol, filename)
        print(f"Molecule data saved as {filename}")

        return send_file(
        filename, as_attachment=True, download_name=os.path.basename(filename)
        )

    except Exception as e:
        print(e)
        return str(e), 500
    
# 3D Structure 3(MOL2)
def save_as_3d_mol2(mol, filename="mol-3d.mol2"):
  if mol is None:
        return False
  try:
      conv = openbabel.OBConversion()
      conv.SetInAndOutFormats("sdf", "mol2")
      mol = openbabel.OBMol()
      sdf_file = save_as_3d_sdf(mol)
      conv.ReadFile(mol, sdf_file)
      conv.WriteFile(mol, "mol-3d.mol2")
      print(f"Molecule data saved as {filename}")
      
      return send_file(
        filename, as_attachment=True, download_name=os.path.basename(filename)
        )
    
  except Exception as e:
    print(e)
    return str(e), 500

# 3D Structure 4(HTML)
def save_as_html(mol, filename="3D-molecule.html"):
    if mol is None:
        return False
    try:
        mol = Chem.AddHs(mol)

        # Crucial fallback: Add RandomCoords if standard embedding fails (common for complex rings)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) < 0:
            AllChem.EmbedMolecule(mol, useRandomCoords=True)

        AllChem.MMFFOptimizeMolecule(mol)
        mol_block = Chem.MolToMolBlock(mol)

        # 2. Initialize viewer (Fixed height typo from 4500 to 450)
        viewer = py3Dmol.view(width=800, height=500)
        viewer.addModel(mol_block, 'mol')

        style = {'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}}

        # 4. Apply style and zoom
        viewer.setStyle(style)
        viewer.zoomTo()

        html_content = viewer._make_html()

        with open(filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        print(f"Successfully saved 3D visualization to {filename}")

        return send_file(
            filename, as_attachment=True, download_name=os.path.basename(filename)
            )
    except Exception as e:
        print(e)
        return str(e), 500
    
    
# 3D Structure 5(JSON)
def save_as_json(mol, filename="mol.json"):
    if mol is None:
        return False

    try:
        mol = Chem.Mol(mol)
        mol = Chem.AddHs(mol)

        status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if status != 0:
            raise RuntimeError("3D coordinate generation failed.")

        AllChem.MMFFOptimizeMolecule(mol)

        # RDKit JSON
        json_str = Chem.MolsToJSON([mol])

        with open(filename, "w") as f:
            f.write(json_str)

        print(f"JSON saved as {filename}")

        return send_file(
        filename, as_attachment=True, download_name=os.path.basename(filename)
        )

    except Exception as e:
        print(e)
        return str(e), 500