import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Config for Multiple Simulations ---
simulations = [
    {"topology": "npt1.gro", "trajectory": "final_pbc1.xtc", "label": "Sim1"},
    {"topology": "npt2.gro", "trajectory": "final_pbc2.xtc", "label": "Sim2"},
    {"topology": "npt3.gro", "trajectory": "final_pbc3.xtc", "label": "Sim3"},
]

ligand_resname = "UNL"
cutoff_distance = 3.5  # Angstroms

# --- Master Residue Set ---
all_residues = set()
residue_counts_by_sim = {}

# --- Loop Over Simulations ---
for sim in simulations:
    universe = mda.Universe(sim["topology"], sim["trajectory"])
    timestep = universe.trajectory.dt / 1000  # ps to ns
    ligand_atoms = universe.select_atoms(f"resname {ligand_resname}")
    protein_atoms = universe.select_atoms(f"protein and not resname {ligand_resname}")
    
    residue_count = {}

    for ts in universe.trajectory:
        distances = mda.lib.distances.distance_array(ligand_atoms.positions, protein_atoms.positions)
        close_atoms_indices = np.where(distances <= cutoff_distance)

        interacting_residues = set()
        for atom_index in close_atoms_indices[1]:
            residue = protein_atoms[atom_index].residue
            res_key = f"{residue.resname}{residue.resid}"
            interacting_residues.add(res_key)

        for residue in interacting_residues:
            all_residues.add(residue)
            residue_count[residue] = residue_count.get(residue, 0) + 1

    residue_counts_by_sim[sim["label"]] = residue_count

# --- Create DataFrame for Plotting ---
all_residues = sorted(all_residues)
df_plot = pd.DataFrame(index=all_residues)

for sim in simulations:
    sim_label = sim["label"]
    counts = residue_counts_by_sim[sim_label]
    df_plot[sim_label] = [counts.get(res, 0) for res in all_residues]

# --- Plotting ---
df_plot.plot(kind='bar', figsize=(16, 6))
plt.title("Residue-Ligand Interaction Frequency Across Simulations")
plt.xlabel("Residue")
plt.ylabel("Number of Frames Interacting")
plt.xticks(rotation=90)
plt.tight_layout()
plt.legend(title="Simulation")
plt.savefig("multi_sim_residue_interactions.png", dpi=300)
plt.show()

