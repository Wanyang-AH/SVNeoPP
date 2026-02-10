#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np

# === Load Data ===
print("🔍 Loading cleaned peptide_summary.csv ...")
df = pd.read_csv("peptide_summary.csv")

# Drop unnecessary columns (if any remnants)
drop_cols = [c for c in df.columns if c.lower().startswith("final")]
df.drop(columns=drop_cols, errors="ignore", inplace=True)

# === Step 1: Basic Descriptive Statistics ===
desc = df[["N_score", "C_score", "max_internal_score", "internal_cleavage_count"]].describe()

# === Step 2: Distribution Skewness and Kurtosis ===
skew = df[["N_score", "C_score", "max_internal_score", "internal_cleavage_count"]].skew()
kurt = df[["N_score", "C_score", "max_internal_score", "internal_cleavage_count"]].kurt()

# === Step 3: Correlation Matrix ===
corr = df[["N_score", "C_score", "max_internal_score", "internal_cleavage_count"]].corr(method="pearson")

# === Step 4: Relationship Check (N + C vs Internal features) ===
df["N_plus_C"] = df["N_score"] + df["C_score"]
corr_NC_internal = df[["N_plus_C", "max_internal_score", "internal_cleavage_count"]].corr()

# === Step 5: Simulate Internal_Risk under varying α, β ===
alpha_values = [0.5, 1.0, 1.5]
beta_values = [0.1, 0.3, 0.5]
sim_results = []
for a in alpha_values:
    for b in beta_values:
        risk = a * df["max_internal_score"] + b * np.log1p(df["internal_cleavage_count"])
        sim_results.append({
            "alpha": a,
            "beta": b,
            "mean_risk": np.mean(risk),
            "std_risk": np.std(risk),
            "p25": np.percentile(risk, 25),
            "p50": np.percentile(risk, 50),
            "p75": np.percentile(risk, 75)
        })
risk_table = pd.DataFrame(sim_results)

# === Step 6: Write everything to a text file ===
with open("internal_risk_diagnostic.md", "w") as f:
    f.write("### Internal Risk Coefficient Diagnostics\n\n")
    f.write("Data Shape: {}\n\n".format(df.shape))
    f.write("== Basic Descriptive Statistics ==\n")
    f.write(desc.to_string())
    f.write("\n\n== Skewness ==\n")
    f.write(skew.to_string())
    f.write("\n\n== Kurtosis ==\n")
    f.write(kurt.to_string())
    f.write("\n\n== Pearson Correlation Matrix ==\n")
    f.write(corr.to_string())
    f.write("\n\n== Correlation with (N + C) ==\n")
    f.write(corr_NC_internal.to_string())
    f.write("\n\n== Simulated Internal_Risk Summary (varying α, β) ==\n")
    f.write(risk_table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    f.write("\n\n# Notes:\n")
    f.write("- Use median and interquartile range to evaluate appropriate penalty scale.\n")
    f.write("- Choose α, β where risk distribution aligns with observed N/C range (0–1).\n")
    f.write("- High correlation between max_internal_score and ICC may suggest smaller β.\n")
    f.write("- Weak correlation implies β should increase for balanced penalization.\n")

print("✅ Diagnostic report saved as 'internal_risk_diagnostics.md'.")

