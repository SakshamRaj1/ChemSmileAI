# ChemSmileAI: a web app for Molecular Analysis and Similarity Search Engine

<div align="center">

[![License](https://img.shields.io/badge/License-MIT-007EC6.svg)](LICENSE)
[![Maintained](https://img.shields.io/badge/Maintained%3F-yes-007EC6.svg)](#)
[![Deployment](https://img.shields.io/badge/Deployment-Render-46E3B7?logo=render&logoColor=white)](https://chemsmileai.onrender.com/)
[![APIs](https://img.shields.io/badge/APIs-None%20(100%25%20Offline%2FLocal)-4c1.svg)](#)

[![Python](https://img.shields.io/badge/Python-3.13.9%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Flask](https://img.shields.io/badge/Flask-3.1.3-000000?logo=flask&logoColor=white)](https://flask.palletsprojects.com/)
[![RDKit](https://img.shields.io/badge/RDKit-2026.3.3-3843d0.svg)](https://www.rdkit.org/)
[![FPSim2](https://img.shields.io/badge/FPSim2-0.7.4-FF6F00.svg)](https://github.com/chembl/FPSim2)
[![Open Babel](https://img.shields.io/badge/Open_Babel-3.2.1-1B365D.svg)](https://openbabel.org/)

[![py3Dmol](https://img.shields.io/badge/py3Dmol-2.5.5-009688.svg)](https://3dmol.csb.pitt.edu/)
[![Openclatura](https://img.shields.io/badge/Openclatura-0.2.0-8E44AD.svg)](https://github.com/lamalab-org/openclatura/)
[![CairoSVG](https://img.shields.io/badge/CairoSVG-2.9.0-E05D44.svg)](https://cairosvg.org/)
[![Ketcher](https://img.shields.io/badge/Ketcher-3.17.0-007EC6.svg)](https://github.com/epam/ketcher)

[![Bootstrap 5](https://img.shields.io/badge/Bootstrap-5-7952B3?logo=bootstrap&logoColor=white)](https://getbootstrap.com/)
[![Pandas](https://img.shields.io/badge/Pandas-3.0.3-150458?logo=pandas&logoColor=white)](https://pandas.pydata.org/)
[![NumPy](https://img.shields.io/badge/NumPy-2.5.0-013243?logo=numpy&logoColor=white)](https://numpy.org/)
[![Gunicorn](https://img.shields.io/badge/Gunicorn-21.2.0-499848?logo=gunicorn&logoColor=white)](https://gunicorn.org/)

</div>

ChemSmileAI is open-sourced web-based computational cheminformatics platform designed to provide accessible, code-free chemical analysis, molecular property computation, and structural modification workflows. ChemSmileAI is a research project for researchers, chemists and enthusiasts for providing a streamlined, web-based platform aimed at simplifying access to chemical data and analysis. Powered by Python, Flask, and RDKit the platform translates molecular inputs (SMILES, CHEMBL IDs, or interactive chemical sketches) into standard molecular descriptors, fingerprint representations, structural visualizations, and similarity search from 2.9 million CHEMBL database compounds.

🌐 **Live Demo:** [Click here to view fully functional web-app](https://tinyurl.com/chemsmileai) [Click here to view web-app without Backend](https://chemsmileai.onrender.com)

<img width="1878" height="896" alt="image" src="https://github.com/user-attachments/assets/e1fb31b9-bace-4945-80da-1e8e692b6bce" />


## System Overview & Architecture

The application is structured as a full-stack cheminformatics pipeline that integrates client-side interactive drawing engines with server-side chemical algorithms and property calculators.


## Core Capabilities

### 1. Molecule Analysis and Information
* Analyse the molecular structures and inspect how atoms, bonds, and functional groups are arranged, how the fingerprint, complexities, Scaffold Analysis convey, with detailed informatics calculated/predicted using RDKit and see the visual representation.
* Molecules information such as Name, Common Name, IUPAC Name and CHEMBL ID are also listed.
* It consists of 12 main functions that contributes to the analysis capabilities.
  
  <img width="1883" height="949" alt="image" src="https://github.com/user-attachments/assets/55f8b11d-faf0-429c-946a-8a976e69d400" />


### 2. Molecule Similarity Search Engine
* Find similar compounds from millions of compounds with related structural patterns as of query molecule using POPCOUNT and investigate how closely molecules resemble one another.
* Quick similarity search from millions of moleculas an get SMILES, CHEMBL ID, Tanimoto Similarity Score, 2D SVG as well.
  
  <img width="1904" height="1080" alt="image" src="https://github.com/user-attachments/assets/94ca57df-dec3-4809-b43c-a3aeb522efb6" />


### 3. Molecule Comparison
* Compare two molecules side by side using 31 functions/properties to evaluate the molecular properties and identify differences or similarities in their structures, descriptors, fingerprints, substructure and other properties. Molecules' information such as Name, Common Name, IUPAC Name and CHEMBL ID are also listed.
  
  <img width="1904" height="1080" alt="image" src="https://github.com/user-attachments/assets/3fbd73ee-9e9b-4458-9ece-8504a93de2af" />


### 4. Molecule Sketching
* Interactive 2D structure (Ketcher editor) with a variety of built in features to draw custom molecular structures directly in the browser and output SMILES, MOL or other format data.
* Use SMILES string directly to view and manipulate the structure.
  
  <img width="1896" height="944" alt="image" src="https://github.com/user-attachments/assets/baab6b96-7b22-4c0d-90ab-f487d27753b8" />


### 5. Delete/Replace Fragment (Substructure Analysis)
* Investigate shared structural fragments based on input SMILES.
* Easily delete unwanted fragments or replace chemical core groups based on input SMARTS/SMILES.
* Highlighted SVG, new SMILES are generated as well.


### 6. Molecule Fingerprint
* Represent molecular structures as fingerprints into 2048 bits. Using Morgan Fingerprints (ECFP4: radius of 2): Length of onbits and Index of onbits are computed.


### 7. Molecular Descriptors
* Calculate over 35+ descriptors (1D, 1D, 2D and 3D descriptors) including molecular weight, hydrogen bond donors/acceptors, surface area, refractivity, and many more. Examine molecular descriptors.


### 8. 2D & 3D Conformer Generation
* Generate high-resolution vector 2D SVG publication-ready diagrams and 3D spatial conformers optimized with ETKDG and MMFF force field calculations.
* SVG/PNG image export for 2D structures and Interactive 3D conformer viewer.


### 9. Download Feature
* Access to 9 file formats for ease of information of a molecule. SDF, MOL, MOL2, PDF, SVG, PNG, HTML/JS, CSV, and JSON

---

## Technical Stack

* **Backend & Web Framework:** Python, Flask, Jinja2 Template Engine
* **Cheminformatics Core:** RDKit, FPSim2 Engine, Openclatura
* **2D and 3D Structures:** py3Dmol, openbabel, cairosvg
* **Chemical Editor:** Ketcher Interactive Drawing Canvas
* **Frontend Architecture:** HTML5, CSS3 (Glassmorphism UI, Responsive CSS Grid & Flexbox), Bootstrap 5, JS (ES6)
* **Current Deployment:** Render Cloud Application Platform (with limitations)
* **APIs:** No external API or website dependency



---
## Installation & Local Setup

### Prerequisites
* Python 3.13.9+
* `pip` package manager
* Recommended: Conda / Miniconda (for optimal C++ library dependency resolution with RDKit)

### Step 1: Clone the Repository
```
# The Backend Data Files and Codes are available on Github but commented in the source code due to Memory limitations of the Deployment Services by Render.
# Once cloned, check for comments at the top of .py files "Uncomment the below lines if cloning!" and uncomment them bofore proceeding to Step 2
```
```
https://github.com/SakshamRaj1/ChemSmileAI.git
cd ChemSmileAI
```

### Step 2: Set Up a Virtual Environment
## Using Conda
```
conda create -n chemsmile python=3.13.9
conda activate chemsmile

# Or using standard venv
python -m venv venv
source venv/bin/activate

# On Windows use:
venv\Scripts\activate
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
```
# Open the localhost (depending on your local setup)
http://127.0.0.1:8000
```
