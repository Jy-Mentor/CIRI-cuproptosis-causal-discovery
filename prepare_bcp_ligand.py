import os
import argparse
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
import pandas as pd

def download_pubchem_sdf(cid, output_file):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        with open(output_file, 'wb') as f:
            f.write(response.content)
        return True
    except Exception as e:
        print(f"Failed to download from PubChem: {e}")
        return False

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
    raw_sdf = os.path.join(output_dir, "BCP_raw.sdf")
    optimized_sdf = os.path.join(output_dir, "BCP_3D.sdf")
    pdb_file = os.path.join(output_dir, "BCP.pdb")
    mol2_file = os.path.join(output_dir, "BCP.mol2")
    pdbqt_file = os.path.join(output_dir, "BCP.pdbqt")
    info_file = os.path.join(output_dir, "ligand_info.txt")
    
    print(f"Output files will be saved to: {output_dir}")
    
    # Step 1: Generate from SMILES (bypassing SDF issues)
    print("Generating molecule from SMILES...")
    mol = generate_3d_from_smiles("C=C1CCC2C(C1)C2(C)CCC=C(C)C")
    
    # Also try to download from PubChem for reference
    success = download_pubchem_sdf(cid, raw_sdf)
    if success:
        print("PubChem SDF downloaded successfully (for reference)")
    else:
        print("Failed to download PubChem SDF, using SMILES-generated structure")
    
    if mol is None:
        print("Failed to obtain molecule")
        return
    
    # Step 2: Optimization and format conversion
    mol = Chem.AddHs(mol)
    
    if minimize:
        print("Performing energy minimization...")
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500, nonBondedThresh=100.0)
        except Exception as e:
            print(f"Energy minimization failed: {e}")
    
    # Save SDF
    writer = Chem.SDWriter(optimized_sdf)
    writer.write(mol)
    writer.close()
    
    # Save PDB
    writer = Chem.PDBWriter(pdb_file)
    writer.write(mol)
    writer.close()
    
    # Save MOL2
    try:
        writer = Chem.MolWriter.Mol2Writer(mol2_file)
        writer.write(mol)
        writer.close()
    except Exception as e:
        print(f"Failed to save MOL2: {e}")
    
    # Save PDBQT using Meeko
    try:
        # Configure Meeko to preserve rotatable bonds
        preparer = MoleculePreparation(
            merge_non_rot_atom_types=False,
            remove_salts=False,
            hydrate=False
        )
        preparer.prepare(mol)
        preparer.write_pdbqt(pdbqt_file)
        
        # Verify TORSDOF value
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('TORSDOF'):
                    torsdof = int(line.strip().split()[-1])
                    print(f"TORSDOF value: {torsdof}")
                    if torsdof == 0:
                        print("Warning: No rotatable bonds detected. Checking molecule structure...")
                    break
    except Exception as e:
        print(f"Failed to save PDBQT: {e}")
    
    # Step 3: BCP specific analysis
    rotatable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    print(f"Number of rotatable bonds: {rotatable_bonds}")
    
    # Step 4: Generate ligand info
    with open(info_file, 'w') as f:
        f.write("BCP Ligand Information\n")
        f.write("====================\n\n")
        f.write(f"SMILES: {Chem.MolToSmiles(mol)}\n")
        f.write(f"InChIKey: {Chem.MolToInchiKey(mol)}\n")
        f.write(f"Molecular Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}\n")
        f.write(f"Molecular Weight: {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}\n")
        f.write(f"LogP: {Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]:.2f}\n")
        f.write(f"Rotatable Bonds: {rotatable_bonds}\n")
        f.write(f"H-Bond Donors: {Chem.rdMolDescriptors.CalcNumHBD(mol)}\n")
        f.write(f"H-Bond Acceptors: {Chem.rdMolDescriptors.CalcNumHBA(mol)}\n")
    
    # Validate PDB file
    validate_pdb(pdb_file)
    
    print("\nProcess completed!")
    print(f"Generated files:")
    print(f"- Raw SDF: {raw_sdf}")
    print(f"- Optimized SDF: {optimized_sdf}")
    print(f"- PDB: {pdb_file}")
    print(f"- MOL2: {mol2_file}")
    print(f"- PDBQT: {pdbqt_file}")
    print(f"- Ligand info: {info_file}")

def validate_pdb(pdb_file):
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith('ATOM') or line.startswith('HETATM')]
        if len(atom_lines) == 0:
            print("Warning: No ATOM records found in PDB file")
        else:
            print(f"PDB file contains {len(atom_lines)} atom records")
            
            # Check coordinate range
            coords = []
            for line in atom_lines:
                if len(line) >= 54:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append((x, y, z))
                    except:
                        pass
            
            if coords:
                min_x = min([c[0] for c in coords])
                max_x = max([c[0] for c in coords])
                min_y = min([c[1] for c in coords])
                max_y = max([c[1] for c in coords])
                min_z = min([c[2] for c in coords])
                max_z = max([c[2] for c in coords])
                
                print(f"Coordinate range: X [{min_x:.2f}, {max_x:.2f}], Y [{min_y:.2f}, {max_y:.2f}], Z [{min_z:.2f}, {max_z:.2f}]")
    except Exception as e:
        print(f"Failed to validate PDB file: {e}")

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