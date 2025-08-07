import sys
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ======= Get Residue Info from CLI ======= #
res_name = sys.argv[1].upper()
res_id = int(sys.argv[2])
ligand_resname = "UNL"  # Update if needed

# ======= Simulation Files ======= #
simulations = [
    {"name": "Sim 1", "tpr": "npt1.tpr", "xtc": "final_pbc1.xtc"},
    {"name": "Sim 2", "tpr": "npt2.tpr", "xtc": "final_pbc2.xtc"},
    {"name": "Sim 3", "tpr": "npt3.tpr", "xtc": "final_pbc3.xtc"},
]

# ======= Collect Distance Data ======= #
all_data = []

for sim in simulations:
    print(f"🔄 Processing {sim['name']} for {res_name}{res_id}...")

    u = mda.Universe(sim["tpr"], sim["xtc"])
    ligand = u.select_atoms(f"resname {ligand_resname}")
    residue = u.select_atoms(f"resname {res_name} and resid {res_id}")

    if residue.n_atoms == 0:
        print(f"⚠️ No atoms found for {res_name}{res_id} in {sim['name']}. Skipping.")
        continue

    distances_list = []

    for ts in u.trajectory:
        dist_matrix = distances.distance_array(residue.positions, ligand.positions)
        min_distance = np.min(dist_matrix)
        distances_list.append(min_distance)

    all_data.extend([
        {"Distance": d, "Simulation": sim["name"]} for d in distances_list
    ])

# ======= Plotting ======= #
if not all_data:
    print(f"❌ No distance data for {res_name}{res_id}. No plot generated.")
    sys.exit()

plot_df = pd.DataFrame(all_data)

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.histplot(
    data=plot_df,
    x="Distance",
    hue="Simulation",
    kde=True,
    stat="probability",
    common_norm=False,
    bins=50,
    palette="Set1",
    alpha=0.6
)

plt.xlabel("Distance (Å)")
plt.ylabel("Probability")
plt.title(f"Distance Distribution: {res_name}{res_id} - Ligand")
plt.tight_layout()
output_file = f"{res_name.lower()}{res_id}_distance_distribution.png"
plt.savefig(output_file, dpi=300)
plt.close()

print(f"✅ Saved: {output_file}")

