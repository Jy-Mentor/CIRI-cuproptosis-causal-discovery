import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

def generate_3d_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Failed to parse SMILES")
        return None
    
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    success = AllChem.EmbedMolecule(mol, params)
    if success != 0:
        print("Failed to generate 3D coordinates")
        return None
    
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500, nonBondedThresh=100.0)
    return mol

def process_bcp_ligand(output_dir, cid=5281515, ph=7.4, minimize=True):
    # Ensure output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory created: {output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        return
    
    # Use absolute paths
    output_dir = os.path.abspath(output_dir)
    pdbqt_file = os.path.join(output_dir, "BCP.pdbqt")
    info_file = os.path.join(output_dir, "ligand_info.txt")
    
    print(f"Output files will be saved to: {output_dir}")
    
    # Generate molecule from SMILES
    print("Generating molecule from SMILES...")
    mol = generate_3d_from_smiles("C=C1CCC2C(C1)C2(C)CCC=C(C)C")
    
    if mol is None:
        print("Failed to obtain molecule")
        return
    
    # Additional energy minimization if requested
    if minimize:
        print("Performing energy minimization...")
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500, nonBondedThresh=100.0)
        except Exception as e:
            print(f"Energy minimization failed: {e}")
    
    # Save PDBQT using Open Babel
    try:
        print("Generating PDBQT file with rotatable bonds...")
        import subprocess
        
        # Save temporary PDB file
        temp_pdb = os.path.join(output_dir, "BCP_temp.pdb")
        writer = Chem.PDBWriter(temp_pdb)
        writer.write(mol)
        writer.close()
        
        # Convert to PDBQT using Open Babel
        subprocess.run(["obabel", temp_pdb, "-opdbqt", "-O", pdbqt_file], check=True)
        print(f"PDBQT file saved (via Open Babel): {pdbqt_file}")
        
        # Clean up
        os.remove(temp_pdb)
        
        # Verify TORSDOF value
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('TORSDOF'):
                    torsdof = int(line.strip().split()[-1])
                    print(f"TORSDOF value: {torsdof}")
                    if torsdof > 0:
                        print("✓ Rotatable bonds detected")
                    else:
                        print("⚠ No rotatable bonds detected")
                    break
    except Exception as e:
        print(f"Failed to generate PDBQT: {e}")
        # Try using Meeko's correct API
        try:
            print("Trying Meeko approach...")
            from meeko import MoleculePreparation
            
            # Use Meeko's MoleculePreparation
            preparer = MoleculePreparation()
            preparer.prepare(mol)
            
            # Get PDBQT string
            pdbqt_string = preparer.write_pdbqt_string()
            
            # Write to file
            with open(pdbqt_file, 'w') as f:
                f.write(pdbqt_string)
            
            print(f"PDBQT file saved (via Meeko): {pdbqt_file}")
            
            # Verify TORSDOF value
            with open(pdbqt_file, 'r') as f:
                for line in f:
                    if line.startswith('TORSDOF'):
                        torsdof = int(line.strip().split()[-1])
                        print(f"TORSDOF value: {torsdof}")
                        if torsdof > 0:
                            print("✓ Rotatable bonds detected")
                        else:
                            print("⚠ No rotatable bonds detected")
                        break
        except Exception as e2:
            print(f"Meeko approach failed: {e2}")
            return
    
    # Generate ligand info
    print("Generating ligand information...")
    try:
        with open(info_file, 'w') as f:
            f.write("BCP Ligand Information\n")
            f.write("====================\n\n")
            f.write(f"SMILES: {Chem.MolToSmiles(mol)}\n")
            f.write(f"InChIKey: {Chem.MolToInchiKey(mol)}\n")
            f.write(f"Molecular Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}\n")
            f.write(f"Molecular Weight: {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}\n")
            f.write(f"LogP: {Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]:.2f}\n")
            f.write(f"Rotatable Bonds: {Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)}\n")
            f.write(f"H-Bond Donors: {Chem.rdMolDescriptors.CalcNumHBD(mol)}\n")
            f.write(f"H-Bond Acceptors: {Chem.rdMolDescriptors.CalcNumHBA(mol)}\n")
        print(f"Ligand info saved: {info_file}")
    except Exception as e:
        print(f"Failed to save ligand info: {e}")
    
    print("\nProcess completed!")
    print(f"Generated files:")
    print(f"- PDBQT: {pdbqt_file}")
    print(f"- Ligand info: {info_file}")

def main():
    parser = argparse.ArgumentParser(description='Prepare BCP ligand for molecular docking')
    parser.add_argument('--output_dir', default='C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\output\\ligands\\', help='Output directory')
    parser.add_argument('--cid', type=int, default=5281515, help='PubChem CID')
    parser.add_argument('--ph', type=float, default=7.4, help='pH for protonation')
    parser.add_argument('--minimize', type=bool, default=True, help='Perform energy minimization')
    args = parser.parse_args()
    
    process_bcp_ligand(args.output_dir, args.cid, args.ph, args.minimize)

if __name__ == "__main__":
    main()