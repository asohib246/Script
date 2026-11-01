from pymatgen.io.vasp import Poscar
from pymatgen.core import Species, PeriodicSite
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core.structure import Structure
import os

# Step 1: Load the Structure from POSCAR (LATP structure)
poscar = Poscar.from_file("POSCAR")  # Replace with your POSCAR file for LATP
structure = poscar.structure

# Step 2: Assign Oxidation States (for LATP)
oxidation_states = {"Li": +1, "Ti": +4, "P": +5, "O": -2}
structure.add_oxidation_state_by_element(oxidation_states)

# Step 3: Identify Ti Sites
Ti_indices = [i for i, site in enumerate(structure.sites) if site.specie.symbol == "Ti"]

# Step 4: Read the updated 36f Wyckoff positions for Li atoms from the 36f-position.txt file
file_path_36f_positions = '36f-position.txt'

# Read the 36f positions from the file
additional_Li_positions = []
with open(file_path_36f_positions, 'r') as file:
    lines = file.readlines()
    for line in lines:
        # The line format is: number x y z (so we split and take the x, y, z)
        parts = line.split()
        x, y, z = map(float, parts[1:4])  # Extract x, y, z coordinates
        additional_Li_positions.append([x, y, z])

# Step 5: Create Output Directory for POSCAR Files
output_dir = "Al_doped_Ti_LATP_with_Li_all_combinations"
os.makedirs(output_dir, exist_ok=True)

# Step 6: Loop Over Ti Sites and Replace One with Al, Add Each Li Atom Position, Calculate Ewald Energy
ewald_results = []  # To store results
structure_counter = 1  # This will keep track of the structure number

for Ti_index in Ti_indices:
    # Make a copy of the original structure to modify
    structure_copy = structure.copy()
    
    # Replace the Ti atom with an Al atom at the Ti_index
    structure_copy.replace(Ti_index, "Al")  # Replace Ti with Al
    
    # Loop through each updated 36f Wyckoff Li position
    for j, Li_position in enumerate(additional_Li_positions):
        # Step 7: Remove the last Li atom (if any) from the structure to ensure only one Li is added
        if len(structure_copy.sites) > len(structure.sites):  # If a Li atom exists
            structure_copy.remove_sites([len(structure_copy.sites) - 1])  # Remove the last Li atom added
        
        # Add the single Li atom at the current 36f Wyckoff position
        structure_copy.append("Li", Li_position)
        
        # Update oxidation state for Al (assuming Al is +3) and Li (+1)
        structure_copy.add_oxidation_state_by_element({"Al": +3, "Li": +1, "Ti": +4, "P": +5, "O": -2})
        
        # Step 8: Perform Ewald Summation for the modified structure
        ewald = EwaldSummation(structure_copy)
        ewald_energy = ewald.total_energy
        
        # Step 9: Save the modified structure to a POSCAR file with the correct naming scheme (using hyphens instead of spaces)
        modified_poscar_filename = os.path.join(output_dir, f"Structure-{structure_counter}_Ti-{Ti_index+1}_replaced_by_Al_with_additional_Li_at_the_{j+1}_position_of_36f")
        Poscar(structure_copy).write_file(modified_poscar_filename)
        
        # Store results (index, energy, filename)
        ewald_results.append((structure_counter, ewald_energy, modified_poscar_filename))
        
        # Print the details of the modification
        print(f"Structure-{structure_counter} Ti-{Ti_index+1} replaced by Al with additional Li at the {j+1} position of 36f: Ewald Energy = {ewald_energy:.4f} eV")
        
        # Increment the structure counter
        structure_counter += 1

# Step 10: Sort the results by Ewald Energy (ascending order)
ewald_results.sort(key=lambda x: x[1])  # Sort by Ewald energy

# Step 11: Save Results to File
results_file = "Al_doped_Ti_LATP_with_Li_all_combinations_ewald_energies.txt"
with open(results_file, "w") as f:
    f.write("Structure_Index\tEwald_Energy\tModified_POSCAR_Filename\n")
    for idx, energy, filename in ewald_results:
        f.write(f"{idx}\t{energy:.6f}\t{filename}\n")

# Print the number of generated structures
print(f"Generated {len(ewald_results)} Al-doped Ti structures with additional Li atom for all combinations.")
print(f"Ewald energies saved to {results_file}.")
