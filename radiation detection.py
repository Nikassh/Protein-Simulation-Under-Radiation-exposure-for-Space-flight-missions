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
print(f"Raw doses (mGy) as strings: {doses_mgy}")
doses_gy = [float(re.search(r'\d+\.?\d*', str(dose)).group()) / 1000 for dose in doses_mgy if re.search(r'\d+\.?\d*', str(dose))]
print(f"Converted doses (Gy): {doses_gy}")

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


#AI training

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import re
    

# Load the Excel file
file_path = r"C:\Users\Nikassh\AppData\Local\Programs\Python\Python313\osdr-environmental-data-rr-radiation-data-5-5-23.xlsx"
df = pd.read_excel(file_path)
print("Excel data preview:\n", df.head())

# Extract and convert doses (mGy to Gy)
doses_mgy = df["Total Absorbed Dose (Experiment)"].dropna().tolist()
doses_gy = [float(re.search(r'\d+\.?\d*', str(dose)).group()) / 1000 for dose in doses_mgy if re.search(r'\d+\.?\d*', str(dose))]
print(f"Converted doses (Gy): {doses_gy}")

# Create placeholder damages (cyclic 5–25% for now)
num_samples = len(doses_gy)
damages = [5 + min((d * 10 / max(doses_gy)), 10) for d in doses_gy]
print(f"Estimated damages (%): {damages[:num_samples]}")

# Features (Dose) and Target (Damage)
X = np.array(doses_gy).reshape(-1, 1)  # 2D array for scikit-learn
y = np.array(damages[:num_samples])

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize and train the model
model = LinearRegression()
model.fit(X_train, y_train)

new_dose = np.array([[0.03]])
predicted_damage = model.predict(new_dose)
print(f"Predicted damage for 0.03 Gy: {predicted_damage[0]}%")


# Predict on test set
y_pred = model.predict(X_test)

# Evaluate
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)
print(f"Mean Squared Error: {mse}")
print(f"R-squared Score: {r2}")
print(f"Model Coefficients: {model.coef_}, Intercept: {model.intercept_}")

# Visualize
plt.scatter(X_train, y_train, color='blue', label='Training data')
plt.plot(X_test, y_pred, color='red', label='Regression line')
plt.xlabel("Dose (Gy)")
plt.ylabel("Damage (%)")
plt.title("Radiation Dose vs. Protein Damage (Excel Data)")
plt.legend()
plt.show()

# Test with a new dose
new_dose = np.array([[0.03]])
predicted_damage = model.predict(new_dose)
print(f"Predicted damage for 0.03 Gy: {predicted_damage[0]}%")













