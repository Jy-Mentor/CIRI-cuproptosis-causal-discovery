import os
import time
import argparse
import requests
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select

class CleanSelect(Select):
    def accept_residue(self, residue):
        resname = residue.get_resname()
        if resname in ['HOH', 'H2O', 'WAT']:
            return False
        if residue.id[0] != ' ':  # Skip heteroatoms
            return False
        return True

    def accept_atom(self, atom):
        if atom.get_altloc() not in [' ', 'A']:
            return False
        return True

def download_pdb(pdb_id, output_dir):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    filename = os.path.join(output_dir, f"{pdb_id.lower()}.pdb")
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        with open(filename, 'w') as f:
            f.write(response.text)
        
        if os.path.getsize(filename) == 0:
            os.remove(filename)
            return False, f"Empty file downloaded for {pdb_id}"
        
        if 'ATOM' not in response.text:
            os.remove(filename)
            return False, f"No ATOM records in {pdb_id}"
        
        return True, filename
    except Exception as e:
        return False, str(e)

def parse_pdb_metadata(pdb_file):
    metadata = {
        'PDB_ID': os.path.basename(pdb_file).split('.')[0].upper(),
        'Title': '',
        'Resolution': '',
        'Organism': '',
        'Expression_System': '',
        'Ligands_Original': '',
        'Chain_Count': 0,
        'Residue_Count': 0
    }
    
    chains = set()
    residues = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HEADER'):
                pass
            elif line.startswith('TITLE'):
                metadata['Title'] += line[10:].strip() + ' '
            elif line.startswith('COMPND'):
                if 'ORGANISM' in line:
                    metadata['Organism'] += line.split('ORGANISM:')[1].strip() + ' '
                elif 'EXPRESSION_SYSTEM' in line:
                    metadata['Expression_System'] += line.split('EXPRESSION_SYSTEM:')[1].strip() + ' '
            elif line.startswith('SOURCE'):
                pass
            elif line.startswith('REMARK 2'):
                if 'RESOLUTION' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'ANGSTROMS':
                            metadata['Resolution'] = parts[i-1]
                            break
            elif line.startswith('HET'):
                lig = line[17:20].strip()
                if lig not in ['HOH', 'H2O', 'WAT'] and lig:
                    metadata['Ligands_Original'] += lig + ', '
            elif line.startswith('ATOM'):
                chain = line[21:22].strip()
                chains.add(chain)
                res_id = (chain, line[22:26].strip(), line[26:27].strip())
                if res_id not in residues:
                    residues.append(res_id)
    
    metadata['Chain_Count'] = len(chains)
    metadata['Residue_Count'] = len(residues)
    metadata['Title'] = metadata['Title'].strip()
    metadata['Organism'] = metadata['Organism'].strip()
    metadata['Expression_System'] = metadata['Expression_System'].strip()
    metadata['Ligands_Original'] = metadata['Ligands_Original'].strip().rstrip(',')
    
    return metadata

def clean_pdb(input_file, output_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_file)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, CleanSelect())
    
    return len(list(structure.get_atoms()))

def main():
    parser = argparse.ArgumentParser(description='Download and clean PDB structures')
    parser.add_argument('--output_dir', default='C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\output\\receptors\\', help='Output directory')
    parser.add_argument('--pdb_ids', default='3O3U,5N2S,6QCB,1ZMC,2FLU,3RZF,8VMC,3P1M,9SOQ,6DVH', help='Comma-separated PDB IDs')
    args = parser.parse_args()
    
    output_dir = args.output_dir
    pdb_ids = [pdb.strip() for pdb in args.pdb_ids.split(',')]
    
    os.makedirs(output_dir, exist_ok=True)
    
    metadata_list = []
    
    log_file = os.path.join(output_dir, 'download.log')
    with open(log_file, 'w') as log:
        log.write(f"Download started at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        for pdb_id in pdb_ids:
            log.write(f"\nProcessing {pdb_id}...\n")
            
            success, result = download_pdb(pdb_id, output_dir)
            if not success:
                log.write(f"Failed to download {pdb_id}: {result}\n")
                metadata = {
                    'PDB_ID': pdb_id,
                    'Title': '',
                    'Resolution': '',
                    'Organism': '',
                    'Expression_System': '',
                    'Ligands_Original': '',
                    'Chain_Count': 0,
                    'Residue_Count': 0,
                    'Download_Status': f"Failed: {result}"
                }
                metadata_list.append(metadata)
                time.sleep(0.5)
                continue
            
            input_file = result
            output_file = input_file
            
            try:
                initial_atoms = clean_pdb(input_file, output_file)
                final_atoms = len(list(PDBParser(QUIET=True).get_structure('structure', output_file).get_atoms()))
                
                file_size = os.path.getsize(output_file) / 1024  # KB
                log.write(f"Downloaded: {pdb_id} ({file_size:.2f} KB)\n")
                log.write(f"Atoms before cleaning: {initial_atoms}\n")
                log.write(f"Atoms after cleaning: {final_atoms}\n")
                
                metadata = parse_pdb_metadata(output_file)
                metadata['Download_Status'] = 'Success'
                metadata_list.append(metadata)
                
            except Exception as e:
                log.write(f"Error processing {pdb_id}: {str(e)}\n")
                metadata = {
                    'PDB_ID': pdb_id,
                    'Title': '',
                    'Resolution': '',
                    'Organism': '',
                    'Expression_System': '',
                    'Ligands_Original': '',
                    'Chain_Count': 0,
                    'Residue_Count': 0,
                    'Download_Status': f"Error: {str(e)}"
                }
                metadata_list.append(metadata)
            
            time.sleep(0.5)
        
        log.write(f"\nDownload finished at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    if metadata_list:
        df = pd.DataFrame(metadata_list)
        csv_file = os.path.join(output_dir, 'pdb_metadata.csv')
        df.to_csv(csv_file, index=False)
        print(f"Metadata saved to: {csv_file}")
    
    print(f"Log file saved to: {log_file}")
    print("Process completed!")

if __name__ == "__main__":
    main()