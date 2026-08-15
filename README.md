# ChemSmileAI: a web app for Molecular Analysis and Similarity Search Engine

ChemSmileAI is a web-based computational cheminformatics platform designed to provide accessible, code-free chemical analysis, molecular property computation, and structural modification workflows. Powered by Python, Flask, and RDKit, the platform translates molecular inputs (SMILES, Molecule Names, ChEMBL IDs, or interactive chemical sketches) into standard molecular descriptors, fingerprint representations, and structural visualizations.

---

## System Overview & Architecture

The application is structured as a full-stack cheminformatics pipeline that integrates client-side interactive drawing engines with server-side chemical algorithms and property calculators.


## Core Capabilities

### 1.Molecule Analysis and Information
* Analyse the molecular structures and inspect how atoms, bonds, and functional groups are arranged, how the fingerprint, complexities, Scaffold Analysis convey, with detailed informatics calculated/predicted using RDKit and see the visual representation. Molecules information such as Name, COmmon Name, IUPAC Name and CHEMBL ID are also listed.

### 2. Molecule Similarity Search Engine
* Find similar compounds from millions of compounds with related structural patterns as of query molecule using POPCOUNT and investigate how closely molecules resemble one another. Quick similarity search from millions of moleculas an get SMILES, CHEMBL ID, Tanimoto Similarity Score, 2D SVG as well.

### 3. Molecule Comparison
* Compare two molecules side by side to evaluate the molecular properties and identify differences or similarities in their structures, descriptors, fingerprints, substructure and other properties. Molecules' information such as Name, Common Name, IUPAC Name and CHEMBL ID are also listed.

### 4. Molecule Sketching
* Interactive 2D structure (Ketcher editor) with a variety of built in features to draw custom molecular structures directly in the browser and output SMILES, MOL or other format data. Use SMILES string directly to view and manipulate the structure.

### 5. Delete/Replace Fragment (Substructure Analysis)
* Investigate shared structural fragments based on input SMILES. Easily delete unwanted fragments or replace chemical core groups based on input SMARTS/SMILES. Highlighted SVG, new SMILES are generated as well.

### 6. Molecule Fingerprint
* Represent molecular structures as fingerprints into 2048 bits. Using Morgan Fingerprints (ECFP4: radius of 2): Length of onbits and Index of onbits are computed.

### 7. Molecular Descriptors
* Calculate over 35+ descriptors (1D, 1D, 2D and 3D descriptors) including molecular weight, hydrogen bond donors/acceptors, surface area, refractivity, and many more. Examine molecular descriptors.

### 8. 2D & 3D Conformer Generation
* Generate high-resolution vector 2D SVG publication-ready diagrams and 3D spatial conformers optimized with ETKDG and MMFF force field calculations. SVG/PNG image export for 2D structures and Interactive 3D conformer viewer.

### 9. Download Feature
* Access to 9 file formats for ease of information of a molecule. SDF, MOL, MOL2, PDF, SVG, PNG, HTML/JS, CSV, and JSON

---

## Technical Stack

* **Backend & Web Framework:** Python, Flask, Jinja2 Template Engine
* **Cheminformatics Core:** RDKit, FPSim2 Engine, Openclatura
* **2D and 3D Structures:** py3Dmol, openbabel, cairosvg
* **Chemical Editor:** Ketcher Interactive Drawing Canvas
* **Frontend Architecture:** HTML5, CSS3 (Glassmorphism UI, Responsive CSS Grid & Flexbox), Bootstrap 5, JavaScript (ES6)
* **Deployment Target:** Render Cloud Application Platform
* **APIs:** No external API or website dependency



---
## Installation & Local Setup

### Prerequisites
* Python 3.13.9+
* `pip` package manager
* Recommended: Conda / Miniconda (for optimal C++ library dependency resolution with RDKit)

### Step 1: Clone the Repository
```
https://github.com/SakshamRaj1/ChemSmileAI.git
cd ChemSmileAI
```

### Step 2: Set Up a Virtual Environment
## Using Conda
```
conda create -n chemsmile python=3.10
conda activate chemsmile

# Or using standard venv
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
```


### Step 3: Install Dependencies
```
pip install -r requirements.txt
```

### Step 4: Run the Application
```
flask run --host=127.0.0.1 --port=8000
```
### Step 5: Access ChemSmileAI locally on your system

Open http://127.0.0.1:8000 (depending on your local setup)

