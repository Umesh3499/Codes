import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rms import RMSF

# ====== List of (TPR, XTC) pairs ====== #
simulations = [
    ("md1.tpr", "final_pbc1.xtc"),
    ("md2.tpr", "final_pbc2.xtc"),
    ("md3.tpr", "final_pbc3.xtc")
]

labels = ["Set 1", "Set 2", "Set 3"]
colors = ["blue", "green", "red"]

plt.figure(figsize=(10, 6))

# ====== Loop through each run ====== #
for (top, traj), label, color in zip(simulations, labels, colors):
    u = mda.Universe(top, traj)

    ca_atoms = u.select_atoms("name CA")
    residue_ids = ca_atoms.resids

    # Compute RMSF
    rmsf = RMSF(ca_atoms).run()

    # Plot RMSF
    plt.plot(residue_ids, rmsf.rmsf, label=label, color=color, linewidth=1.5)

# ====== Plot Styling ====== #
plt.xlabel("Residue Number", fontsize=12)
plt.ylabel("RMSF (Å)", fontsize=12)
plt.title("Residue-wise RMSF (Cα Atoms) from 3 Simulations", fontsize=14)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("combined_rmsf_plot.png", dpi=300)
plt.show()

