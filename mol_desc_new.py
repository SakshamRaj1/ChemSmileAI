def molecular_descriptors(mol):
    try:
        # 0D Descriptors
        MolecularWeight= Descriptors.MolWt(mol)
        LogP= Descriptors.MolLogP(mol)
        HBondDonors= Descriptors.NumHDonors(mol)
        HBondAcceptors= Descriptors.NumHAcceptors(mol)
        RotatableBonds= Descriptors.NumRotatableBonds(mol)
        tPSA= Descriptors.TPSA(mol)
        HeavyAtomCount= Descriptors.HeavyAtomCount(mol)
        ValenceElectrons= Descriptors.NumValenceElectrons(mol)
        RadialElectrons= Descriptors.NumRadicalElectrons(mol)
        AromaticRings= Descriptors.NumAromaticRings(mol)

        # 1D Descriptors
        MolarRefractivity= Crippen.MolMR(mol)
        AromaticRings= Descriptors.NumAromaticRings(mol)
        HBondDonors= Descriptors.NumHDonors(mol)
        HBondAcceptors= Descriptors.NumHAcceptors(mol)
        RotatableBonds= Descriptors.NumRotatableBonds(mol)
        tPSA= Descriptors.TPSA(mol)
        HeavyAtomCount= Descriptors.HeavyAtomCount(mol)
        ValenceElectrons= Descriptors.NumValenceElectrons(mol)
        
        # 2D Descriptors
        
        BalabanIndex= GraphDescriptors.BalabanJ(mol)
        BertzCT= GraphDescriptors.BertzCT(mol)
        bCUT2D= rdMolDescriptors.BCUT2D(mol)
        Chi0v= Descriptors.Chi0v(mol)
        Chi1v= Descriptors.Chi1v(mol)
        Chi2v= Descriptors.Chi2v(mol)
        Chi3v= Descriptors.Chi3v(mol)
        Chi4v= Descriptors.Chi4v(mol)
        Chi0n= Descriptors.Chi0n(mol)
        Chi1n= Descriptors.Chi1n(mol)
        Chi2n= Descriptors.Chi2n(mol)
        Chi3n= Descriptors.Chi3n(mol)
        Chi4n= Descriptors.Chi4n(mol)
        Kappa1= Descriptors.Kappa1(mol)
        Kappa2= Descriptors.Kappa2(mol)
        Kappa3= Descriptors.Kappa3(mol)
        
        
    except Exception as e:
        return f"Error calculating in: {str(e)}"

    return{
        # 0D Descriptors
        "MolecularWeight": MolecularWeight,
        "LogP": LogP,
        "HBondDonors": HBondDonors,
        "HBondAcceptors": HBondAcceptors,
        "RotatableBonds": RotatableBonds,
        "TPSA": tPSA,
        "HeavyAtomCount": HeavyAtomCount,
        "ValenceElectrons": ValenceElectrons,
        "RadialElectrons": RadialElectrons,
        "AromaticRings": AromaticRings,

        # 1D Descriptors
        "MolarRefractivity": MolarRefractivity,
        "AromaticRings": AromaticRings,
        "HBondDonors": HBondDonors,
        "HBondAcceptors": HBondAcceptors,
        "RotatableBonds": RotatableBonds,
        "TPSA": tPSA,
        "HeavyAtomCount": HeavyAtomCount,
        "ValenceElectrons": ValenceElectrons,
        
        # 2D Descriptors
        
        "BalabanJIndex": BalabanIndex,
        "BertzCT": BertzCT,
        "BCUT2D": bCUT2D,
        "Chi0v": Chi0v,
        "Chi1v": Chi1v,
        "Chi2v": Chi2v,
        "Chi3v": Chi3v,
        "Chi4v": Chi4v,
        "Chi0n": Chi0n,
        "Chi1n": Chi1n,
        "Chi2n": Chi2n,
        "Chi3n": Chi3n,
        "Chi4n": Chi4n,
        "Kappa1": Kappa1,
        "Kappa2": Kappa2,
        "Kappa3": Kappa3,
        
        # For 3d descriptors, I created a separate function for better calculation.
        "molecular_descriptors3d": molecular_descriptors3d(mol)
    }
