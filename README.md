ChemSmileAI: Computational Cheminformatics and Molecular Analysis PlatformChemSmileAI is a web-based computational cheminformatics platform designed to provide accessible, code-free chemical analysis, molecular property computation, and structural modification workflows. Powered by Python, Flask, and RDKit, the platform translates molecular inputs (SMILES, Molecule Names, ChEMBL IDs, or interactive chemical sketches) into standard molecular descriptors, fingerprint representations, and structural visualizations.System Overview & ArchitectureThe application is structured as a full-stack cheminformatics pipeline that integrates client-side interactive drawing engines with server-side chemical algorithms and property calculators.[ Client Interface / Ketcher Sketcher / Search Input ]
                        │
                        ▼ (HTTP GET / POST)
        [ Flask Application Route Handlers ]
                        │
       ┌────────────────┴────────────────┐
       ▼                                 ▼
[ RDKit Cheminformatics Engine ]  [ Local Metadata Store / Cache ]
 - Descriptors (MW, LogP, TPSA)    - Compound Name & Synonyms
 - Fingerprints (Morgan, MACCS)    - IUPAC & ChEMBL Identifiers
 - 2D Depiction (SVG)              - Molecule Parent Entities
 - Substructure Search & Replace
                        │
                        ▼
    [ Server-Side Jinja2 HTML / SVG Rendering ]
Core Capabilities1. Molecular Property Computation & Descriptor PipelinePhysicochemical Descriptors: Evaluates fundamental properties including Molecular Weight ($\text{MW}$), Wildman-Crippen calculated partition coefficient ($\text{LogP}$), Topological Polar Surface Area ($\text{TPSA}$), and Rotatable Free Bonds.Drug-Likeness Evaluation: Automates rule-based screening, specifically validating the Lipinski Rule of 5 criteria (Hydrogen Bond Donors/Acceptors, $\text{MW} \le 500\text{ Da}$, $\text{LogP} \le 5$) for potential drug candidates.Topological Complexity: Measures atom/bond distribution matrices, ring structures, and topological complexity indices.2. Fingerprinting & Molecular Similarity SearchBit Vector Generation: Generates structural fingerprints including Morgan Fingerprints (Circular Fingerprints / ECFP representations) and MACCS structural key pattern sets.Similarity Algorithms: Executes structural analog searches and distance evaluations across molecules using Tanimoto and Dice similarity coefficients.3. Structural Comparison & AlignmentMulti-Target Comparison: Provides side-by-side comparative analysis of multiple compounds.Property Delta Mapping: Evaluates variances in $\text{LogP}$, topological polar surface area, and molecular weight across series of analogs.Alignment & RMSD Computation: Structural coordinate alignment and Root-Mean-Square Deviation (RMSD) evaluation.4. Interactive Chemical Sketching & Structure EditingKetcher Canvas Integration: Features an in-browser molecular editor allowing users to draw and modify chemical structures, outputting SMILES strings directly into the computational pipeline.Substructure Search & Fragment Replacement: Identifies specific sub-graphs or functional groups within target molecules to perform fragment deletion or core replacement (lead optimization/R-group exploration).5. 2D & 3D Conformer Generation2D Vector Graphics: Generates publication-ready, clean SVG 2D structural depictions with highlighted substructure mappings.3D Coordinate Generation: Generates 3D spatial conformers via distance geometry (ETKDG) and energy minimization routines using force fields (MMFF94).Technical StackBackend & Web Framework: Python, Flask, Jinja2 Template EngineCheminformatics Core: RDKit (Python C++ Bindings)Chemical Editor: Ketcher Interactive Drawing CanvasFrontend Architecture: HTML5, CSS3 (Glassmorphism UI, Responsive CSS Grid & Flexbox), Bootstrap 5, JavaScript (ES6)Deployment Target: Render Cloud Application PlatformRepository StructureChemSmileAI/
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
Installation & Local SetupPrerequisitesPython 3.9+pip package managerRecommended: Conda / Miniconda (for optimal C++ library dependency resolution with RDKit)Step 1: Clone the RepositoryBashgit clone https://github.com/SakshamRaj1/ChemSmileAI.git
cd ChemSmileAI
Step 2: Set Up a Virtual EnvironmentBash# Using Conda
conda create -n chemsmile python=3.10
conda activate chemsmile

# Or using standard venv
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
Step 3: Install DependenciesBashpip install -r requirements.txt
Step 4: Run the ApplicationBashpython app.py
Open [http://127.0.0.1:5000](http://127.0.0.1:5000) (or http://localhost:9000 depending on your local port configuration) in your web browser.Application Route ReferenceRouteHTTP MethodDescription/GET, POSTPrimary search interface and computed descriptor dashboard/analysisGET, POSTDetailed atom/bond metrics, Lipinski compliance, and descriptor breakdown/similarityGET, POSTTanimoto/Dice molecular similarity search engine/compareGET, POSTComparative analysis between two target molecules/aboutGETPlatform architecture and documentation/contactGET, POSTUser inquiries and contact form
