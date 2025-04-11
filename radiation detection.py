from Bio.PDB import PDBParser, PPBuilder
from Bio import Entrez
import ssl
import pandas as pd
import re
from time import sleep

# Parse the protein file
parser = PDBParser()
structure = parser.get_structure("1a2b", r"C:\Users\Nikassh\Downloads\1a2b.pdb")

for model in structure:
    for chain in model:
        print(f"Chain {chain.id} has {len(chain)} residues")

# Extract sequence and count secondary structures
ppb = PPBuilder()
sequence = ""
helix_count = 0
sheet_count = 0

for pp in ppb.build_peptides(structure[0]):
    sequence += str(pp.get_sequence())
    helix_count += sequence.count("H")  # Placeholder
    sheet_count += sequence.count("E")  # Placeholder

print(f"Protein sequence length: {len(sequence)}")
print(f"Estimated helices: {helix_count}, Estimated sheets: {sheet_count}")

# Load the downloaded Excel file from GeneLab
data = pd.read_excel(r"C:\Users\Nikassh\AppData\Local\Programs\Python\Python313\osdr-environmental-data-rr-radiation-data-5-5-23.xlsx")

# Extract total absorbed doses (convert mGy to Gy)
doses_mgy = data["Total Absorbed Dose (Experiment)"].dropna().tolist()
doses_gy = [dose / 1000 for dose in doses_mgy]  # Convert mGy to Gy
print(f"Extracted doses (Gy): {doses_gy}")

# Estimate damages based on low-dose effects
num_samples = max(len(doses_gy), 5)
damages = [5, 10, 15, 20, 25][:num_samples]  # Default range for now
print(f"Estimated damages (%): {damages}")

# Create updated radiation dataset
num_samples = max(len(doses_gy), len(damages), 5)
data = {
    "Dose (Gy)": doses_gy[:num_samples] + [doses_gy[-1] + 0.001] * (num_samples - len(doses_gy)) if doses_gy else [0.005, 0.010, 0.015, 0.020, 0.025],
    "Damage (%)": damages[:num_samples] + [damages[-1] + 5] * (num_samples - len(damages)) if damages else [5, 10, 15, 20, 25],
    "Residue_Count": [218] * num_samples,
    "Helix_Count": [helix_count] * num_samples,
    "Sheet_Count": [sheet_count] * num_samples
}

df = pd.DataFrame(data)
df.to_csv("radiation_dataset.csv", index=False)
print("\nUpdated dataset saved to radiation_dataset.csv:")
print(df)
