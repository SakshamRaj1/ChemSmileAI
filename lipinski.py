def lipinski_table(mol):
    mw = Descriptors.MolWt(mol)
    violation1, violation2, violation3, violation4 = 0, 0, 0, 0

    for i in range(1):  # Loop to allow early exit on first violation

        if mol is not None:
            mw = Descriptors.MolWt(mol)
            if mw > 500: 
                violation1 = 1
                print(f"Fail: MW = ({mw} > 500 Daltons)")
            else:
                print(f"Pass: MW = ({mw} <= 500 Daltons)")
            
            logp = Crippen.MolLogP(mol)
            if logp > 5: 
                violation2 = 1
                print(f"Fail: logp = ({logp} > 5)")
            else:
                print(f"Pass: logp = ({logp} <= 5)")

            hbd = rdMolDescriptors.CalcNumHBD(mol)
            if hbd > 5: 
                violation3 = 1
                print(f"Fail: HBD = ({hbd} > 5)")
            else:
                print(f"Pass: HBD = ({hbd} <= 5)")

            hba = rdMolDescriptors.CalcNumHBA(mol)
            if hba > 10:
                violation4 = 1
                print(f"Fail: HBA = ({hba} > 10)")
            else:
                print(f"Pass: HBA = ({hba} <= 10)")
            
        else:
            print("Invalid molecule. Please check the SMILES string and try again.")

    return f"Pass: All Lipinski's Rule of 5 criteria met" if sum([violation1, violation2, violation3, violation4]) <= 1 else f"Fail: Lipinski's Rule of 5 criteria NOT met {sum([violation1, violation2, violation3, violation4])} Violations"