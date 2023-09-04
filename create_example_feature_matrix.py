import pandas as pd
import numpy as np

# Define the metabolites and their names
metabolites = ["Metabolite1", "Metabolite2", "Metabolite3", "Metabolite4"]

# Create a list of sample names with concentrations from 1 to 1000 mM in 50 mM steps
concentrations = list(range(0, 1100, 100))

# Create an empty DataFrame to store the data
data = pd.DataFrame(columns=["metabolite"] + [f"{sample_name} mM" for sample_name in concentrations])

# Generate sigmoidal intensity values for each metabolite in each sample
for i, metabolite in enumerate(metabolites):
    metabolite_data = {"metabolite": metabolite}
    for concentration in concentrations:
        # Generate sigmoidal intensity values (you can adjust parameters as needed)
        intensity = int((1 / (1 + np.exp(-(concentration - 400+(50*(i+1))) / 100))) * 100000 * (i+1))  # Sigmoidal function example
        metabolite_data[f"{concentration} mM"] = intensity
    data.loc[metabolite] = metabolite_data

data["control"] = [5000, 400, 70000, 8000]
data["treatment"] = [60000, 7500, 250000, 7500]
# Save the DataFrame to a TSV file
data.to_csv("metabolite_intensity_data.tsv", sep="\t", index=False)