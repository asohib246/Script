from pymatgen.io.vasp import Poscar
from pymatgen.analysis.ewald import EwaldSummation

# Step 1: Load the Structure from POSCAR
poscar = Poscar.from_file("POSCAR")  # Replace with your POSCAR file path
structure = poscar.structure

# Step 2: Assign Oxidation States (for LATP)
oxidation_states = {"Li": +1, "Ti": +4, "P": +5, "O": -2}
structure.add_oxidation_state_by_element(oxidation_states)

# Step 3: Check if all sites have oxidation states
for site in structure.sites:
    if site.specie.symbol not in oxidation_states:
        print(f"Warning: No oxidation state defined for {site.specie.symbol} at {site.coords}.")

# Step 4: Calculate Ewald Energy for the structure
ewald = EwaldSummation(structure)

# Step 5: Calculate the total Ewald energy
ewald_energy = ewald.total_energy

# Output the Ewald energy
print(f"Ewald energy for the structure: {ewald_energy:.4f} eV")
