from Bio.PDB import PDBParser, PPBuilder, Selection
from Bio.SeqUtils import ProtParam
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import re

# Parse the protein file
parser = PDBParser()
structure = parser.get_structure("AF-P04637-F1-model_v4", r"C:\Users\Nikassh\Downloads\AF-P04637-F1-model_v4.pdb")

# Function to estimate secondary structure (moved up to avoid NameError)
def estimate_secondary_structure(seq):
    helix_count = 0
    sheet_count = 0
    for i in range(len(seq) - 4):
        if seq[i] in "AILVF" and seq[i+3] in "AILVF":
            helix_count += 1
        if seq[i] in "AILVF" and seq[i+1] in "DEKRH":
            sheet_count += 1
    return helix_count // 5, sheet_count // 3

# Filter chains with too many missing residues
ppb = PPBuilder()
valid_chains = []
for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        print(f"Chain {chain.id} has {len(residues)} residues")
        if len(residues) > 30:  # Adjusted threshold for smaller chains
            valid_chains.append(chain)

# Rebuild sequence from valid chains and calculate chain-specific features
sequence = ""
chain_features = []
for chain in valid_chains:
    chain_sequence = ""
    peptides = ppb.build_peptides(chain)
    if not peptides:
        print(f"Skipping Chain {chain.id} due to empty sequence (likely discontinuous).")
        continue
    for pp in peptides:
        chain_sequence += str(pp.get_sequence())
    if not chain_sequence:
        print(f"Skipping Chain {chain.id} due to empty sequence after peptide building.")
        continue
    h, s = estimate_secondary_structure(chain_sequence)
    chain_analyzer = ProtParam.ProteinAnalysis(chain_sequence)
    chain_aa = chain_analyzer.amino_acids_percent
    chain_sensitive_aa = chain_aa.get("C", 0) + chain_aa.get("Y", 0)
    chain_features.append({
        "chain_id": chain.id,
        "sensitive_aa": chain_sensitive_aa,
        "helix_count": h,
        "sheet_count": s,
        "sequence": chain_sequence
    })
    sequence += chain_sequence

print(f"Filtered protein sequence length: {len(sequence)}")
print("Chain-specific features:")
for cf in chain_features:
    print(f"Chain {cf['chain_id']}: Sensitive AA: {cf['sensitive_aa']:.2f}%, Helices: {cf['helix_count']}, Sheets: {cf['sheet_count']}")

# Check if we have valid chains
if not chain_features:
    raise ValueError("No valid chains found with usable sequences. Please check the PDB file for issues.")

# Load the Excel file
data = pd.read_excel(r"C:\Users\Nikassh\AppData\Local\Programs\Python\Python313\osdr-environmental-data-rr-radiation-data-5-5-23.xlsx")

# Extract total absorbed doses (convert mGy to Gy)
doses_mgy = data["Total Absorbed Dose (Experiment)"].dropna().tolist()
doses_gy = [float(re.search(r'\d+\.?\d*', str(dose)).group()) / 1000 for dose in doses_mgy if re.search(r'\d+\.?\d*', str(dose))]
print(f"Converted doses (Gy): {doses_gy}")

# Radiation type with more variability
radiation_type = np.random.choice([0, 0.5, 1], len(doses_gy), p=[0.4, 0.2, 0.4])
print(f"Radiation types encoded: {radiation_type.tolist()}")

# Exposure time correlated with dose
exposure_times = np.array(doses_gy) * 1000 + np.random.normal(0, 1, len(doses_gy))
print(f"Exposure times (days): {exposure_times.tolist()}")

# Realistic damage values with chain-specific factors
def dose_response(dose, radiation_type, exposure_time, sensitive_aa, max_damage=15, k=500, x0=0.008):
    base_damage = max_damage / (1 + np.exp(-k * (dose - x0)))
    rad_factor = 1.2 if radiation_type == 0 else 1.0 if radiation_type == 1 else 0.9
    time_factor = 1 + np.log1p(exposure_time) * 0.05
    aa_factor = 1 + sensitive_aa * 0.02
    adjusted_damage = base_damage * rad_factor * time_factor * aa_factor
    noise = np.random.normal(0, 1.5)
    return max(5, min(max_damage, adjusted_damage + noise))

# Assign chain-specific features to each data point
damages = []
sensitive_aas = []
helix_counts = []
sheet_counts = []
for i in range(len(doses_gy)):
    cf = np.random.choice(chain_features)
    sensitive_aas.append(cf["sensitive_aa"])
    helix_counts.append(cf["helix_count"])
    sheet_counts.append(cf["sheet_count"])
    damage = dose_response(doses_gy[i], radiation_type[i], exposure_times[i], cf["sensitive_aa"])
    damages.append(damage)

# Features
X_extended = np.column_stack((
    np.array(doses_gy).reshape(-1, 1),
    np.array(radiation_type).reshape(-1, 1),
    np.array(sensitive_aas),
    exposure_times,
    np.array(helix_counts),
    np.array(sheet_counts)
))

# Target
y = np.array(damages)

# Augment with synthetic data
n_synthetic = 100
synthetic_doses = np.random.uniform(min(doses_gy), max(doses_gy), n_synthetic)
synthetic_rad_types = np.random.choice([0, 0.5, 1], n_synthetic, p=[0.4, 0.2, 0.4])
synthetic_exposure = np.array(synthetic_doses) * 1000 + np.random.normal(0, 1, n_synthetic)

synthetic_damages = []
synthetic_sensitive_aas = []
synthetic_helix_counts = []
synthetic_sheet_counts = []
for i in range(n_synthetic):
    cf = np.random.choice(chain_features)
    synthetic_sensitive_aas.append(cf["sensitive_aa"])
    synthetic_helix_counts.append(cf["helix_count"])
    synthetic_sheet_counts.append(cf["sheet_count"])
    damage = dose_response(synthetic_doses[i], synthetic_rad_types[i], synthetic_exposure[i], cf["sensitive_aa"])
    synthetic_damages.append(damage)

X_synthetic = np.column_stack((
    synthetic_doses,
    synthetic_rad_types,
    np.array(synthetic_sensitive_aas),
    synthetic_exposure,
    np.array(synthetic_helix_counts),
    np.array(synthetic_sheet_counts)
))
X_extended = np.vstack((X_extended, X_synthetic))
y = np.concatenate((y, synthetic_damages))

# Split data
X_train_ext, X_test_ext, y_train, y_test = train_test_split(X_extended, y, test_size=0.2, random_state=42)

# Hyperparameter tuning for Random Forest with stricter regularization
param_grid = {
    'n_estimators': [100, 200],
    'max_depth': [3, 5],
    'min_samples_split': [10, 20],
    'min_samples_leaf': [8, 12]
}
grid_search = GridSearchCV(RandomForestRegressor(random_state=42), param_grid, cv=5, scoring='r2', n_jobs=-1)
grid_search.fit(X_extended, y)
print(f"Best parameters: {grid_search.best_params_}")
print(f"Best R-squared: {grid_search.best_score_:.2f}")
rf_model = grid_search.best_estimator_

# Evaluate
y_pred_rf = rf_model.predict(X_test_ext)
mse_rf = mean_squared_error(y_test, y_pred_rf)
r2_rf = r2_score(y_test, y_pred_rf)
print(f"Updated Random Forest MSE: {mse_rf:.2f}, R-squared: {r2_rf:.2f}")

# Feature importance
feature_names = ["Dose", "Radiation Type", "Sensitive AA", "Exposure Time", "Helix Count", "Sheet Count"]
importances = rf_model.feature_importances_
for name, importance in zip(feature_names, importances):
    print(f"{name}: {importance:.3f}")

# Visualize
plt.scatter(X_train_ext[:, 0], y_train, color='blue', label='Training data')
plt.scatter(X_test_ext[:, 0], y_pred_rf, color='green', label='Random Forest Prediction', alpha=0.6)
plt.xlabel("Dose (Gy)")
plt.ylabel("Protein Damage (%)")
plt.title("Radiation Dose vs. Protein Damage")
plt.legend()
plt.show()

# Chain-specific damage prediction
chain_damages = []
for cf in chain_features:
    input_data = np.array([[0.01, 0.5, cf["sensitive_aa"], 5.0, cf["helix_count"], cf["sheet_count"]]])
    damage = rf_model.predict(input_data)[0]
    chain_damages.append((cf["chain_id"], damage))

# Plot chain-specific damage
chains, damages = zip(*chain_damages)
plt.bar(chains, damages, color='purple')
plt.xlabel("Chain")
plt.ylabel("Predicted Damage (%) at 0.01 Gy")
plt.title("Protein Damage Across Chains")
plt.show()

# Dose-response curve for Chain A
chain_a = next(cf for cf in chain_features if cf["chain_id"] == 'A')
dose_range = np.linspace(0.001, 0.015, 50)
chain_a_damages = []
for dose in dose_range:
    input_data = np.array([[dose, 0.5, chain_a["sensitive_aa"], 5.0, chain_a["helix_count"], chain_a["sheet_count"]]])
    damage = rf_model.predict(input_data)[0]
    chain_a_damages.append(damage)

plt.plot(dose_range, chain_a_damages, color='red', marker='o')
plt.xlabel("Dose (Gy)")
plt.ylabel("Predicted Damage (%)")
plt.title("Dose-Response Curve for Chain A")
plt.show()

# Save dataset
data_dict = {
    "Dose (Gy)": X_extended[:, 0],
    "Radiation_Type": X_extended[:, 1],
    "Sensitive_AA": X_extended[:, 2],
    "Exposure_Time": X_extended[:, 3],
    "Helix_Count": X_extended[:, 4],
    "Sheet_Count": X_extended[:, 5],
    "Damage (%)": y
}
df = pd.DataFrame(data_dict)
df.to_csv("radiation_dataset_updated.csv", index=False)
print("\nUpdated dataset saved to radiation_dataset_updated.csv:")
print(df)
