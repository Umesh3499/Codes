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

plt.figure(figsize=(10, 6))

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
        res_atom = res_atoms[min_dist_idx[0]]
        lig_atom = ligand[min_dist_idx[1]]
        distance = dist_matrix[min_dist_idx]

        distances_min.append(distance)

    # Save to CSV
    df = pd.DataFrame({"Frame": frame_nums, "Distance": distances_min})
    df.to_csv(f"{output_dir}/leu146_distance_sim{i}.csv", index=False)

    # Plot
    plt.plot(frame_nums, distances_min, label=labels[i-1], color=colors[i-1])

# Finalize plot
plt.xlabel("Frame Number")
plt.ylabel("Minimum Distance (Å)")
plt.title("Distance Between LEU146 and Ligand (UNL) Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("leu146_distance_vs_frame_3sims.png", dpi=300)
plt.show()

print("✅ Distance vs. Frame Number plot saved as 'leu146_distance_vs_frame_3sims.png'.")

