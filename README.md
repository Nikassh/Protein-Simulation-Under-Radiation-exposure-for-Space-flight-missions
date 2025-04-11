# 🧬 Protein Radiation Damage Estimation

This project simulates and estimates **radiation damage on protein structures** using Biopython and real radiation dose data from **NASA GeneLab**.

Perfectly brewed for those looking into **bioinformatics + space research + radiation biology** combos. Inspired by future work in space biosciences at agencies like **NASA** and **ISRO**.

---

## 🚀 Project Overview

- Parses a **PDB protein structure** (example: `1a2b.pdb`)
- Calculates structural features like:
  - Total sequence length
  - Helix count
  - Sheet count
- Extracts **absorbed radiation doses** from real **GeneLab datasets**
- Converts doses from mSv → Gy
- Simulates **radiation-induced protein damage** using basic modeling

---

## 🧠 Why This Project?

Astronauts and biological specimens face constant exposure to cosmic radiation. This prototype attempts to:
- Understand secondary structure changes due to space radiation
- Lay the foundation for **predictive modeling** of protein damage in space missions
- Connect **space bioinformatics** with **molecular-level insights**

---

## 🧰 Technologies Used

| Tool        | Purpose                                  |
|-------------|-------------------------------------------|
| Python      | Main programming language                 |
| Biopython   | Parsing PDB files, analyzing structures   |
| Pandas      | Reading GeneLab Excel data               |
| NASA GeneLab | Real-world radiation exposure dataset    |

---

## 📂 File Structure

```bash
.
├── radiation_detection.py    # Main analysis script
├── 1a2b.pdb                  # Example protein structure
├── osdr-radiation-data.xlsx  # Environmental dose dataset
└── README.md                 # This beautiful file
