import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# List of (topology, trajectory) pairs
simulations = [
    ("npt1.gro", "final_pbc1.xtc"),
    ("npt2.gro", "final_pbc2.xtc"),
    ("npt3.gro", "final_pbc3.xtc")
]

labels = ["Set 1", "Set 2", "Set 3"]
colors = ["blue", "green", "red"]

plt.figure(figsize=(10, 6))

# Loop through each topology/trajectory pair
for (top, traj), label, color in zip(simulations, labels, colors):
    u = mda.Universe(top, traj)
    complex_structure = u.select_atoms("protein or resname UNL")

    # Time in ns
    dt = u.trajectory.dt / 1000  # ps to ns
    time_ns = np.arange(len(u.trajectory)) * dt

    # Radius of Gyration in Å
    Rg_values = np.array([complex_structure.radius_of_gyration() * 10 for ts in u.trajectory])

    # Plot
    plt.plot(time_ns, Rg_values, label=label, color=color, linewidth=1.5)

# Plot settings
plt.xlabel("Time (ns)", fontsize=12)
plt.ylabel("Radius of Gyration (Å)", fontsize=12)
plt.title("Radius of Gyration Over Time (Protein + Ligand)", fontsize=14)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("combined_radius_of_gyration.png", dpi=300)
plt.show()

