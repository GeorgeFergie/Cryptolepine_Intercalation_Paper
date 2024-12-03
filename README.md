# Cryptolepine_Intercalation_Paper
Output files from the Cryptolepine Paper alongside ASE (Python), CP2K and GROMACS .mdp files/scripts

# CP2K
Contains an example CP2K optimisation script used for the optimisations performed in the paper.

# MD
Cotains ten .mdp files: four used in the RMSD calculations performed in the paper: Min_Energy, Equil_Temperature and Equil_Pressure were used to prepare the system before MD.mdp was used to run and produce the results. Four, labled with the ending '_FEP' were used to run the FEP calculations. Two are used in the Umbrella Sampling Calculation, 'US' and 'MD_Steered', used to run the Umbrella sampling and Steered MD respectfully. These are scripts for their specific structure and different structures contained minor changes in the .mdp files. 

Finally it contains two forcefield files: "CRYC.itp" and "charmm36.itp". These are used to provide the MD parameters for Cryptolepine in our MD simulations.

# MD Structures
Contains two folders; FEP and RMSD+US. These are the input structures for the respective calculations. FEP contains a salt concentration and was based on the RMSD+US structures. 

# NEB
Contains the final output geometries of all the NEB calculations performed in the paper.

# Optimisation
Contains 3 directories which describe the size of the system: Small, Small_BB and Large_BB. Small contains the optimised geometry of all small structures (Cryptolepine and two base pairs above and below with no connecting backbone). Small_BB is the same as Small but now including the connecting backbone between the base pairs. Large_BB contains the largest geometry optimisations performed with 4 base pairs, two above and below a Cryptolepine molecule with all connecting backbone.

# Python
Contains the NEB script ran to calculate the NEB in ASE via CP2K. Additionally, it contants Thermo_Calc which is how we predicted the occupation of states via the Boltzmann Probability distribution
