# Cryptolepine_Intercalation_Paper
Output files from the Cryptolepine Paper alongside ASE (Python), Qchem input scripts and GROMACS .mdp files/scripts

# MD
Cotains ten .mdp files: four used in the RMSD calculations performed in the paper: Min_Energy, Equil_Temperature and Equil_Pressure were used to prepare the system before MD.mdp was used to run and produce the results. Four, labled with the ending '_FEP' were used to run the FEP calculations. Two are used in the Umbrella Sampling Calculation, 'US' and 'MD_Steered', used to run the Umbrella sampling and Steered MD respectfully. These are scripts for their specific structure and different structures contained minor changes in the .mdp files. 

Finally it contains two forcefield files: "CRYC.itp" and "charmm36.itp". These are used to provide the MD parameters for Cryptolepine in our MD simulations.

# MD Structures
Contains two folders; FEP and RMSD+US. These are the input structures for the respective calculations. FEP contains a salt concentration and was based on the RMSD+US structures. 

# Optimisation - Ring Model
Contains the final output geometries of all NWChem optimsation calculations ran on a ring model of DNA. 

# Python
Contains Thermo_Calc which is how we predicted the occupation of states via the Boltzmann Probability distribution as well as an input script example used to run NWChem via ASE for our optimisations. The NWChem script is set up to run the CC system.

# Qchem
Contains an example input script we used to run the Energy Decomposition Analysis (EDA) calculation via Qchem. The example in the script is for CC.
