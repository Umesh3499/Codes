import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
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

    # Reference from same topology
    ref = mda.Universe(top)
    ref_complex = ref.select_atoms("protein or resname UNL")

    # RMSD analysis
    rmsd_analysis = RMSD(complex_structure, reference=ref_complex, select="all")
    rmsd_analysis.run()

    time_ns = rmsd_analysis.rmsd[:, 1] / 1000  # Convert ps to ns
    rmsd_values = rmsd_analysis.rmsd[:, 2]

    plt.plot(time_ns, rmsd_values, label=label, color=color, linewidth=1.5)

# Plot settings
plt.xlabel("Time (ns)", fontsize=12)
plt.ylabel("RMSD (Å)", fontsize=12)
plt.title("RMSD of Protein-Ligand Complex (3 Simulation Sets)", fontsize=14)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("combined_complex_rmsd.png", dpi=300)
plt.show()

