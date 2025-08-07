import MDAnalysis as mda
import numpy as np
import pandas as pd

# --- Configuration ---
trajectory_file = "final_pbc.xtc"  # Replace with your trajectory file
topology_file = "npt.gro"          # Replace with your topology file
ligand_resname = "UNL"             # Replace with your ligand residue name
cutoff_distance = 3.5              # Distance cutoff in Angstroms
output_file = "interacting_residues.txt"  # Output file to store results

# --- Load Files ---
universe = mda.Universe(topology_file, trajectory_file)
timestep = universe.trajectory.dt / 1000  # Convert ps to ns

# --- Select Ligand Atoms ---
ligand_atoms = universe.select_atoms(f"resname {ligand_resname}")

# --- Initialize Data ---
residue_list = []  # Master list of all interacting residues
frame_data = []    # Store data for each frame
residue_count = {}  # Track residue occurrence count across frames

# --- Iterate Through Frames ---
for ts in universe.trajectory:
    frame_number = universe.trajectory.frame + 1  # Frame number starts from 1
    time_ns = frame_number * timestep  # Time in nanoseconds

    # --- Select Protein Atoms (Excluding Ligand) ---
    protein_atoms = universe.select_atoms(f"protein and not resname {ligand_resname}")

    # --- Calculate Distances ---
    distances = mda.lib.distances.distance_array(ligand_atoms.positions, protein_atoms.positions)
    
    # --- Identify Interacting Residues ---
    interacting_residues = set()
    close_atoms_indices = np.where(distances <= cutoff_distance)

    for atom_index in close_atoms_indices[1]:
        residue = protein_atoms[atom_index].residue
        res_key = f"{residue.resname}{residue.resid}"
        interacting_residues.add(res_key)

    # --- Update Master List of Residues ---
    for residue in interacting_residues:
        if residue not in residue_list:
            residue_list.append(residue)
            residue_count[residue] = 0

    # --- Count Occurrence of Interacting Residues ---
    for residue in interacting_residues:
        residue_count[residue] += 1

    # --- Prepare Frame Data ---
    frame_residues = []
    for residue in residue_list:
        if residue in interacting_residues:
            frame_residues.append(residue)
        else:
            frame_residues.append("None")

    # --- Add Count of Interacting Residues ---
    num_interacting_residues = len(interacting_residues)
    frame_residues.append(num_interacting_residues)

    frame_data.append([f"Frame {frame_number}", round(time_ns, 2)] + frame_residues)

# --- Prepare DataFrame ---
columns = ["Frame", "Time (ns)"] + residue_list + ["No_of_Interacting_Residues"]
df = pd.DataFrame(frame_data, columns=columns)

# --- Add Total Count Row ---
total_count_row = ["Total Count", ""] + [residue_count.get(residue, 0) for residue in residue_list] + [""]
df.loc[len(df.index)] = total_count_row

# --- Save to File ---
df.to_csv(output_file, sep="\t", index=False)
print(f"✅ Analysis complete. Results written to {output_file}")

