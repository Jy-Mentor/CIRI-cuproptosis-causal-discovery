#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import requests
import os
import time

# Create directory for saving PDB files
pdb_dir = "output/new_receptors"
os.makedirs(pdb_dir, exist_ok=True)

# Protein to PDB ID mapping based on user's table (including both options when available)
protein_pdb_map = {
    "RAGE": ["3O3U", "5CJY"],
    "ADORA1": ["5N2S"],
    "CTSD": ["6QCB", "1LYW"],
    "LIAS": ["1ZMC"],
    "NFE2L2": ["2FLU"],
    "IKKβ": ["3RZF", "4KIK"],
    "FASN": ["6C7P", "8VMC"],
    "FDX1": ["3P1M"],
    "TIMP1": ["9SOQ"],
    "TSPO": ["5DVH"]
}

# PDB database API URL
pdb_api_url = "https://files.rcsb.org/download/"

print("Starting download of protein PDB files...")
print("PDB files will be saved to:", os.path.abspath(pdb_dir))
print()

success_count = 0
fail_count = 0
failed_proteins = []

for protein, pdb_ids in protein_pdb_map.items():
    for pdb_id in pdb_ids:
        print("Downloading", protein, "(PDB ID:", pdb_id, ")...")
        
        try:
            # Download PDB file
            pdb_file_url = pdb_api_url + pdb_id + ".pdb"
            pdb_response = requests.get(pdb_file_url, timeout=10)
            pdb_response.raise_for_status()
            
            # Save PDB file (use protein name + pdb id for multiple options)
            if len(pdb_ids) > 1:
                pdb_file_path = os.path.join(pdb_dir, protein.lower() + "_" + pdb_id.lower() + ".pdb")
            else:
                pdb_file_path = os.path.join(pdb_dir, protein.lower() + ".pdb")
            
            with open(pdb_file_path, "wb") as f:
                f.write(pdb_response.content)
            
            print("Successfully downloaded", protein, "PDB:", pdb_file_path)
            success_count += 1
            
        except Exception as e:
            print("Error downloading", protein, "(", pdb_id, "):", str(e))
            fail_count += 1
            failed_proteins.append(f"{protein} ({pdb_id})")
        
        # Avoid too frequent requests
        time.sleep(1)
        print()

print("\nDownload complete!")
print("Successfully downloaded:", success_count, "proteins")
print("Failed to download:", fail_count, "proteins")

if failed_proteins:
    print("\nProteins that failed to download:")
    for protein in failed_proteins:
        print("-", protein)

print("\nPDB files saved to:", os.path.abspath(pdb_dir))
