# Load the structure and trajectory
load MD.gro
load_traj MD.xtc

# Select the protein
select my_protein, polymer

# Select cations within 3.4 Å of the protein
select cations_near_protein, (resname CHL and (br. all within 3.4 of my_protein))

# Select anions within 3.4 Å of the protein
select anions_near_protein, (resname HSO4 and (br. all within 3.4 of my_protein))

# Select water within 3.4 Å of the protein
select water_near_protein, (resname SOL and (br. all within 3.4 of my_protein))

# Remove everything except the protein, cations, and anions
remove not (my_protein or cations_near_protein or anions_near_protein or water_near_protein)

# Print the number of cations and anions near the protein
print ("Anions near protein:", cmd.count_atoms("anions_near_protein"))
print ("Cations near protein:", cmd.count_atoms("cations_near_protein"))
print ("water near protein:", cmd.count_atoms("water_near_protein"))

# Visualization settings
show cartoon, my_protein
color gray90, my_protein

show sticks, cations_near_protein
color red, cations_near_protein

show sticks, anions_near_protein
color green, anions_near_protein

ls show sticks, water_near_protein
color blue, water_near_protein

# Optional: Save the selection for further analysis
save filtered_system.pdb, (my_protein or cations_near_protein or anions_near_protein)

# Optional: Create a high-quality image
set ray_trace_mode, 1
set ray_shadows, 1
ray
png protein_with_ions.png

hide spheres, resname GER
hide sticks, resname CHL
set ray_opaque_background, off

