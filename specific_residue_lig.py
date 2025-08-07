import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd

# ======= Load Trajectory ======= #
u = mda.Universe("npt.tpr", "final_pbc.xtc")  # Update with your filenames

# Select ligand and LYS34
ligand = u.select_atoms("resname 6ZV")  # Update 'LIG' with your ligand residue name
lys34 = u.select_atoms("resname LYS and resid 34")

# ======= Initialize Data Storage ======= #
data = []

# ======= Loop Over Frames ======= #
for ts in u.trajectory:
    frame = ts.frame  # Frame number

    # Get all atoms of LYS34
    res_atoms = lys34.atoms

    # Calculate distances between LYS34 atoms and ligand atoms
    dist_matrix = mda.analysis.distances.distance_array(res_atoms.positions, ligand.positions)

    # Find the minimum distance and corresponding atom indices
    min_dist_idx = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
    res_atom = res_atoms[min_dist_idx[0]]  # Closest LYS34 atom
    lig_atom = ligand.atoms[min_dist_idx[1]]  # Closest ligand atom
    distance = dist_matrix[min_dist_idx]  # Minimum distance

    # Check if the closest atom is a hydrogen
    if res_atom.name.startswith("H"):
        # Get the parent atom if hydrogen
        parent_atom = res_atom.bonded_atoms[0] if res_atom.bonded_atoms else None
        parent_atom_name = parent_atom.name if parent_atom else "Unknown"
    else:
        parent_atom_name = res_atom.name  # Not a hydrogen

    # Store the closest interaction for this frame
    data.append([frame, "LYS", 34, res_atom.name, parent_atom_name, lig_atom.name, distance])

# ======= Convert to DataFrame ======= #
df = pd.DataFrame(data, columns=["Frame", "Residue", "ResidueID", "ResAtom", "ParentAtom", "LigAtom", "Distance"])
df.to_csv("lys34_interaction.csv", index=False)

print("✅ Closest interaction per frame recorded! Results saved in 'lys34_interaction.csv'.")

#====== graph_plot========#
import pandas as pd
import matplotlib.pyplot as plt

# ======= Load Interaction Data ======= #
df = pd.read_csv("lys34_interaction.csv")

# ======= Filter Data: Use Parent Atom If Hydrogen ======= #
# If the interacting atom is hydrogen, use the parent atom for plotting
df["PlotAtom"] = df.apply(lambda row: row["ParentAtom"] if row["ResAtom"].startswith("H") else row["ResAtom"], axis=1)

# ======= Plot Parent Atom Interaction Distance Over Time ======= #
plt.figure(figsize=(8, 5))
plt.plot(df["Frame"], df["Distance"], marker="o", linestyle="-", color="blue", label="Parent-Ligand Distance (Å)")

# ======= Customize Plot ======= #
plt.xlabel("Frame")
plt.ylabel("Distance (Å)")
plt.title("Interaction Distance Between LYS34 Parent Atom and Ligand Over Time")
plt.grid(True)
plt.legend()
plt.savefig("lys34_interaction_plot.png", dpi=300)
plt.show()

print("✅ Graph plotted successfully! Saved as 'lys34_interaction_plot.png'.")

