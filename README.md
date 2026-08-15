# ChemSmileAI: a web app for Molecular Analysis and Similarity Search Engine

ChemSmileAI is a web-based computational cheminformatics platform designed to provide accessible, code-free chemical analysis, molecular property computation, and structural modification workflows. Powered by Python, Flask, and RDKit, the platform translates molecular inputs (SMILES, Molecule Names, ChEMBL IDs, or interactive chemical sketches) into standard molecular descriptors, fingerprint representations, and structural visualizations.

---

## System Overview & Architecture

The application is structured as a full-stack cheminformatics pipeline that integrates client-side interactive drawing engines with server-side chemical algorithms and property calculators.



## Core Capabilities

### 1. Molecular Property Computation & Descriptor Pipeline
* **Physicochemical Descriptors:** Evaluates fundamental properties including Molecular Weight (MW), Wildman-Crippen calculated partition coefficient (LogP), Topological Polar Surface Area (TPSA), and Rotatable Free Bonds.
* **Drug-Likeness Evaluation:** Automates rule-based screening, specifically validating the Lipinski Rule of 5 criteria (Hydrogen Bond Donors/Acceptors, MW <= 500 Da, LogP <= 5) for potential drug candidates.
* **Topological Complexity:** Measures atom/bond distribution matrices, ring structures, and topological complexity indices.

### 2. Fingerprinting & Molecular Similarity Search
* **Bit Vector Generation:** Generates structural fingerprints including Morgan Fingerprints (Circular Fingerprints / ECFP representations) and MACCS structural key pattern sets.
* **Similarity Algorithms:** Executes structural analog searches and distance evaluations across molecules using Tanimoto and Dice similarity coefficients.

### 3. Structural Comparison & Alignment
* **Multi-Target Comparison:** Provides side-by-side comparative analysis of multiple compounds.
* **Property Delta Mapping:** Evaluates variances in LogP, topological polar surface area, and molecular weight across series of analogs.
* **Alignment & RMSD Computation:** Structural coordinate alignment and Root-Mean-Square Deviation (RMSD) evaluation.

### 4. Interactive Chemical Sketching & Structure Editing
* **Ketcher Canvas Integration:** Features an in-browser molecular editor allowing users to draw and modify chemical structures, outputting SMILES strings directly into the computational pipeline.
* **Substructure Search & Fragment Replacement:** Identifies specific sub-graphs or functional groups within target molecules to perform fragment deletion or core replacement (lead optimization/R-group exploration).

### 5. 2D & 3D Conformer Generation
* **2D Vector Graphics:** Generates publication-ready, clean SVG 2D structural depictions with highlighted substructure mappings.
* **3D Coordinate Generation:** Generates 3D spatial conformers via distance geometry (ETKDG) and energy minimization routines using force fields (MMFF94).

---

## Technical Stack

* **Backend & Web Framework:** Python, Flask, Jinja2 Template Engine
* **Cheminformatics Core:** RDKit (Python C++ Bindings)
* **Chemical Editor:** Ketcher Interactive Drawing Canvas
* **Frontend Architecture:** HTML5, CSS3 (Glassmorphism UI, Responsive CSS Grid & Flexbox), Bootstrap 5, JavaScript (ES6)
* **Deployment Target:** Render Cloud Application Platform

---

[ Client Interface / Ketcher Sketcher / Search Input ]
│
▼ (HTTP GET / POST)
[ Flask Application Route Handlers ]
│
┌────────────────┴────────────────┐
▼                                 ▼
[ RDKit Cheminformatics Engine ]  [ Local Metadata Store / Cache ]

Descriptors (MW, LogP, TPSA)    - Compound Name & Synonyms

Fingerprints (Morgan, MACCS)    - IUPAC & ChEMBL Identifiers

2D Depiction (SVG)              - Molecule Parent Entities

Substructure Search & Replace
│
▼
[ Server-Side Jinja2 HTML / SVG Rendering ]

---

## Repository Structure

ChemSmileAI/
├── static/
│   ├── css/               # Glassmorphism themes & responsive stylesheets
│   ├── js/                # Client-side toggles & drawing canvas scripts
│   └── images/            # Platform architecture diagrams and static assets
├── templates/
│   ├── base.html          # Base structural layout, navigation, & offcanvas sidebar
│   ├── index.html         # Main dashboard, search bar, & descriptor results
│   ├── analysis.html      # In-depth property & atom/bond analysis views
│   ├── similarity.html    # Molecular similarity search engine interface
│   ├── compare.html       # Side-by-side molecular comparison workflow
│   ├── about.html         # Project overview & documentation
│   └── contact.html       # Communication and feedback interface
├── app.py                 # Flask server, routing definitions, & RDKit processing pipelines
├── requirements.txt       # Python package dependencies
└── README.md              # Research & technical documentation
