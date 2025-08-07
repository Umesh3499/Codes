import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ========= Simulation Info ========= #
simulations = [
    ("md1.tpr", "final_pbc1.xtc"),
    ("md2.tpr", "final_pbc2.xtc"),
    ("md3.tpr", "final_pbc3.xtc"),
]

output_dir = "distance_outputs"
os.makedirs(output_dir, exist_ok=True)

colors = ['r', 'g', 'b']
labels = ['Simulation 1', 'Simulation 2', 'Simulation 3']

plt.figure(figsize=(12, 6))

# ==== SHADED REGIONS FOR DISTANCE RANGES ====
plt.axhspan(0, 2.0, facecolor='skyblue', alpha=0.3, label="< 2.0 Å (Very Close Contact)")
plt.axhspan(2.0, 3.5, facecolor='lightgreen', alpha=0.3, label="2.0–3.5 Å (Equilibrium Contact)")
plt.axhspan(3.5, 10, facecolor='lightcoral', alpha=0.2, label="> 3.5 Å (Weak/No Contact)")

for i, (tpr, xtc) in enumerate(simulations, 1):
    print(f"▶ Processing LEU146 in Simulation {i}...")

    u = mda.Universe(tpr, xtc)

    ligand = u.select_atoms("resname UNL")
    leu146 = u.select_atoms("resname LEU and resid 146")
    res_atoms = leu146.atoms

    frame_nums = []
    distances_min = []

    for ts in u.trajectory:
        frame = ts.frame
        frame_nums.append(frame)

        dist_matrix = distances.distance_array(res_atoms.positions, ligand.positions)
        min_dist_idx = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
        distance = dist_matrix[min_dist_idx]

        distances_min.append(distance)

    # Save to CSV
    df = pd.DataFrame({"Frame": frame_nums, "Distance": distances_min})
    df.to_csv(f"{output_dir}/leu146_distance_sim{i}.csv", index=False)

    # Plot
    plt.plot(frame_nums, distances_min, label=labels[i-1], color=colors[i-1], linewidth=1.5)

# ==== Finalize Plot ====
plt.xlabel("Frame Number", fontsize=12)
plt.ylabel("Minimum Distance (Å)", fontsize=12)
plt.title("Distance Between LEU146 and Ligand (UNL) Over Time", fontsize=14)
plt.ylim(0, max(distances_min) + 1)
plt.legend(loc="upper right", fontsize=10)
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("leu146_distance_zones_plot.png", dpi=300)
plt.show()

print("✅ Shaded plot for LEU146 saved as 'leu146_distance_zones_plot.png'.")

