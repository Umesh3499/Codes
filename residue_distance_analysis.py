import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
import os

# ========= Simulation Info ========= #
simulations = [
    ("md1.tpr", "final_pbc1.xtc"),
    ("md2.tpr", "final_pbc2.xtc"),
    ("md3.tpr", "final_pbc3.xtc"),
]

residue_list = [
    ("VAL", 95),
    ("LYS", 34),
    ("GLU", 93),
    ("ASP", 98),
    ("ASP", 157),
]

colors = ['r', 'g', 'b']
labels = ['Simulation 1', 'Simulation 2', 'Simulation 3']
output_dir = "distance_outputs"
os.makedirs(output_dir, exist_ok=True)

for resname, resid in residue_list:
    plt.figure(figsize=(12, 6))

    # ==== SHADED CONTACT ZONES ====
    plt.axhspan(0, 2.0, facecolor='skyblue', alpha=0.3, label="< 2.0 Å (Very Close Contact)")
    plt.axhspan(2.0, 3.5, facecolor='lightgreen', alpha=0.3, label="2.0–3.5 Å (Equilibrium Contact)")
    plt.axhspan(3.5, 10, facecolor='lightcoral', alpha=0.2, label="> 3.5 Å (Weak/No Contact)")

    for i, (tpr, xtc) in enumerate(simulations, 1):
        print(f"▶ Processing {resname}{resid} in Simulation {i}...")

        u = mda.Universe(tpr, xtc)

        ligand = u.select_atoms("resname UNL")
        residue = u.select_atoms(f"resname {resname} and resid {resid}")
        res_atoms = residue.atoms

        time_ns = []
        distances_min = []

        for ts in u.trajectory:
            time_ns.append(ts.time / 1000)  # Convert ps to ns

            dist_matrix = distances.distance_array(res_atoms.positions, ligand.positions)
            min_dist_idx = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
            distance = dist_matrix[min_dist_idx]

            distances_min.append(distance)

        # Optional: save distances if needed
        # np.savetxt(f"{output_dir}/{resname}{resid}_distance_sim{i}.txt",
        #            np.column_stack([time_ns, distances_min]),
        #            header="Time_ns Distance")

        # Plot for this sim
        plt.plot(time_ns, distances_min, label=labels[i-1], color=colors[i-1], linewidth=1.5)

    # ==== Finalize Plot for this residue ====
    plt.xlabel("Time (ns)", fontsize=12)
    plt.ylabel("Minimum Distance (Å)", fontsize=12)
    plt.title(f"Minimum Distance: {resname}{resid} - Ligand (UNL)", fontsize=14)
    plt.ylim(0, 10)
    plt.legend(loc="upper right", fontsize=10)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()

    output_file = f"{output_dir}/{resname.lower()}{resid}_distance_plot.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"✅ Plot saved: {output_file}")

