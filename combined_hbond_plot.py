import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# File paths to the hbond.xvg files
simulations = [
    "hbond1.xvg",  # First simulation
    "hbond2.xvg",  # Second simulation
    "hbond3.xvg",  # Third simulation
]

# Initialize an empty list to store the hydrogen bond data
data_combined = []

# Loop through each simulation's hbond file
for i, sim_file in enumerate(simulations, 1):
    data = np.loadtxt(sim_file, comments=["#", "@"])  # Ignore headers
    time = data[:, 0]
    hbond_count = data[:, 1]

    for t, h in zip(time, hbond_count):
        data_combined.append([t, h, f"Simulation {i}"])

# Create DataFrame
df_combined = pd.DataFrame(data_combined, columns=["Time (ns)", "Hbond_Count", "Simulation"])

plt.figure(figsize=(10, 6))

# Color map for clarity
colors = {
    "Simulation 1": "red",
    "Simulation 2": "blue",
    "Simulation 3": "green"
}

# Plot line for each simulation
for sim in df_combined["Simulation"].unique():
    sim_data = df_combined[df_combined["Simulation"] == sim]
    time = sim_data["Time (ns)"]
    hbond = sim_data["Hbond_Count"]

    # Optional: smooth with rolling average
    smoothed = hbond.rolling(window=10, center=True).mean()

    plt.plot(time, smoothed, label=sim, color=colors[sim], linewidth=1.5)

plt.xlabel("Time (ns)")
plt.ylabel("Hydrogen Bonds")
plt.title("Hydrogen Bonds vs Time")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("hbond_lineplot_clean.png", dpi=300)
plt.show()
