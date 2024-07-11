# from ase.calculators.cp2k import CP2K
from ase.calculators.cp2k_noshell import CP2K
from ase.build import molecule
from ase.optimize import BFGS
import ase.io
from ase.neb import NEB, NEBTools
from ase.visualize import view
import time
import copy
import matplotlib.pyplot as plt

def n_plot(xlab, ylab, xs=14, ys=14):
    plt.minorticks_on()
    plt.tick_params(axis='both', which='major', labelsize=ys - 2, direction='in', length=6, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=ys - 2, direction='in', length=4, width=2)
    plt.tick_params(axis='both', which='both', top=True, right=True)
    if xlab is not None:
        plt.xlabel(xlab, fontsize=xs)
    if ylab is not None:
        plt.ylabel(ylab, fontsize=ys)
    plt.tight_layout()
    return None

def cp2k_calc(xc_f="BLYP",
              e_cut=450,
              r_cut=50,
              ngrids=3,
              eps_scf=1.0e-7,
              charge=0,
              multiplicity=1,
              uks=False,
              basis_set="TZV2P-MOLOPT-GTH",
              basis_set_file="BASIS_MOLOPT",
              poisson_solver="AUTO",
              potential_file="POTENTIAL",
              f_disp=None,
              f_solv=None,
              f_choice=None,
              max_scf=100, ):
    """
    Simple CP2K calculator

    Basis set choices
    TZVP-MOLOPT-GTH
    DZVP-MOLOPT-SR-GTH
    TZV2P-MOLOPT-GTH

    Dispersion choices
    https://manual.cp2k.org/trunk/CP2K_INPUT/ATOM/METHOD/XC/VDW_POTENTIAL.html#list_POTENTIAL_TYPE
    https://manual.cp2k.org/trunk/CP2K_INPUT/ATOM/METHOD/XC/VDW_POTENTIAL/NON_LOCAL.html
    https://manual.cp2k.org/trunk/CP2K_INPUT/ATOM/METHOD/XC/VDW_POTENTIAL/PAIR_POTENTIAL.html


    """
    # Determine the pseudo potential
    if xc_f is None:
        potential = None
    elif xc_f.upper() == "B3LYP":
        potential = "GTH-{}".format("BLYP")
    else:
        potential = "GTH-{}".format("BLYP")

    # Check if dispersion correction is requested
    if f_disp is not None and f_disp is not False:
        inpt_disp = '''
                    &VDW_POTENTIAL
                        POTENTIAL_TYPE PAIR_POTENTIAL
                        &PAIR_POTENTIAL
                            REFERENCE_FUNCTIONAL {}
                            PARAMETER_FILE_NAME dftd3.dat
                        &END PAIR_POTENTIAL
                    &END VDW_POTENTIAL'''.format(xc_f.upper())
    else:
        inpt_disp = ""

    # Check if implicit solvent is requested
    if f_solv is not None and f_solv is not False:
        # Check if presetting
        if f_solv.upper() == "WATER":
            dielec_const = 78.5
        elif f_solv.upper() == "PROTEIN":
            dielec_const = 8.0
        else:
            dielec_const = f_solv

        # Make the input
        inpt_solv = '''
            &SCCS  T
                RELATIVE_PERMITTIVITY {}
            &END SCCS'''.format(str(dielec_const))
    else:
        inpt_solv = ""

    if f_choice.upper() == "NEB":
        inpt = '''
        &FORCE_EVAL
            &DFT
                CHARGE +1
                MULTIPLICITY 1
                &MGRID
                    NGRIDS 4
                    CUTOFF 1300
                    REL_CUTOFF 90
                    COMMENSURATE
                &END MGRID
                &SCF
                    SCF_GUESS RESTART
                    EPS_SCF 5.0E-8
                    &OT  T
                        MINIMIZER  CG
                        STEPSIZE   0.15
                        PRECONDITIONER FULL_ALL
                    &END OT
                    &OUTER_SCF  T
                        EPS_SCF 5.0E-8
                    &END OUTER_SCF
                &END SCF
                &XC
                    DENSITY_CUTOFF     1.0E-8
                    GRADIENT_CUTOFF    1.0E-8
                    TAU_CUTOFF         1.0E-8
                    &XC_FUNCTIONAL BLYP
                    &END XC_FUNCTIONAL
                    &VDW_POTENTIAL
                        POTENTIAL_TYPE PAIR_POTENTIAL
                        &PAIR_POTENTIAL
                            REFERENCE_FUNCTIONAL BLYP
                            PARAMETER_FILE_NAME dftd3.dat
                        &END PAIR_POTENTIAL
                    &END VDW_POTENTIAL
                &END XC
                  &POISSON
                    &IMPLICIT
                      BOUNDARY_CONDITIONS PERIODIC
                      MAX_ITER 30
                      &DIELECTRIC
                          DIELECTRIC_CONSTANT 80
                      &END DIELECTRIC
                    &END IMPLICIT
                    PERIODIC XYZ
                  &END POISSON
            &END DFT
        &END FORCE_EVAL
        '''.format(charge, multiplicity, e_cut, r_cut)
    elif f_choice.upper() == "DFTB":
        # This uses DFTB input https://www.cp2k.org/howto:gfn1xtb
        inpt = '''
        &FORCE_EVAL
            &DFT
                CHARGE {}
                MULTIPLICITY {}
                &QS
                    METHOD XTB
                    &XTB
                        CHECK_ATOMIC_CHARGES F
                        DO_EWALD  T
                    &END XTB
                &END QS
                &SCF
                    SCF_GUESS RESTART
                    EPS_SCF {}
                    &OT ON
                        PRECONDITIONER FULL_SINGLE_INVERSE
                        MINIMIZER DIIS
                    &END
                    &OUTER_SCF
                        EPS_SCF {}
                        MAX_SCF {}
                    &END OUTER_SCF
                &END SCF{}
            &END DFT
            &SUBSYS
                &CELL
                     A 2.9881000000000e+01 0.0000000000000e+00 0.0000000000000e+00
                     B 0.0000000000000e+00 2.9246000000000e+01 0.0000000000000e+00
                     C 0.0000000000000e+00 0.0000000000000e+00 3.1873000000000e+01
                &END CELL
            &SUBSYS
        &END FORCE_EVAL
        '''.format(charge, multiplicity, eps_scf, eps_scf, max_scf, inpt_solv)
        # Set the xc functional, and basis set to None
        basis_set = None
        max_scf = None
        basis_set_file = None
        potential_file = None
        potential = None
    else:
        # Reverts to the default input
        inpt = '''
        &FORCE_EVAL
            &DFT
            CHARGE {}
            MULTIPLICITY {}
                &MGRID
                    NGRIDS {}
                    CUTOFF {}
                    REL_CUTOFF {}
                &END MGRID
                &SCF
                    SCF_GUESS RESTART
                    EPS_SCF {}
                    &OUTER_SCF
                        EPS_SCF {}
                        MAX_SCF {}
                    &END
                &END SCF
                &XC
                    &XC_FUNCTIONAL {}
                    &END XC_FUNCTIONAL{}
                &END XC{}
            &END DFT
        &END FORCE_EVAL
        '''.format(charge, multiplicity, ngrids, e_cut, r_cut, eps_scf, eps_scf, max_scf, xc_f, inpt_disp, inpt_solv)

    return CP2K(inp=inpt,
                max_scf=max_scf,
                xc=None,
                cutoff=None,
                print_level='LOW',
                basis_set=basis_set,
                basis_set_file=basis_set_file,
                poisson_solver=poisson_solver,
                charge=0,
                pseudo_potential=potential,
                potential_file=potential_file,
                auto_write=False,
                debug=False,
                uks=uks,
                )


f_max = 0.1
steps = 500
e_cut = 1300
N = 5
i=0
r_cut = 90
xc_f = "BLYP"
f_neb = True
f_neb_plot = True
calc = cp2k_calc(f_choice="NEB") #f_choice="DFTB"

initial = ase.io.read('CC_Start.xyz', index='-1')
final = ase.io.read('CC_Middle.traj', index='-1')

initial.center(vacuum=12.0) # Center structures
final.center(vacuum=12.0)

UpdateName_Init = 'CC_Start_1.traj'
UpdateName_Fin  = 'CC_Middle_1.traj'
f_name_neb = 'CC_Start_to_Middle_NEB.traj'

initial.calc = copy.copy(calc) # Get NEB calculator
final.calc = copy.copy(calc)
final.set_cell(initial.get_cell())

BFGS(initial, trajectory=UpdateName_Init).run(steps=steps, fmax=f_max) # Run optimisations on start and end structure to ensure convergence to force maximum
BFGS(final, trajectory=UpdateName_Fin).run(steps=steps, fmax=f_max)

initialNew = ase.io.read(UpdateName_Init, index=-1) # Read final image of optimisation
finalNew = ase.io.read(UpdateName_Fin, index=-1)

if f_neb:
    # Make the initial band path
    images = [initialNew]
    images += [initialNew.copy() for ii in range(N - 2)]
    images += [finalNew]

    print("Number of images: " + str(len(images)))
    # Set the calculator for all images
    for image in images:
        i+=1
        if i == N:
            image.calc = copy.copy(calc)
            #image.set_cell(initial.get_cell())
            image.get_potential_energy()
            image.get_forces()

        else:
            image.calc = copy.copy(calc)
            image.get_potential_energy()
            image.get_forces()
        
    neb = NEB(images, climb=False)
    neb.interpolate(method='idpp')

    BFGS(neb, trajectory=f_name_neb).run(fmax=f_max)

if f_neb_plot:
    N = NEBTools(ase.io.read(f_name_neb, index=':'))._guess_nimages()

    neb_path = ase.io.read(f_name_neb, index=str(-N) + ':')
    view(neb_path)


    nebtools = NEBTools(neb_path)

    # Get the calculated barrier and the energy change of the reaction
    Ef, dE = nebtools.get_barrier()
    print("Barrier: %.2f eV, reaction energy: %.2f eV" % (Ef, dE))
    # Get the maximum force
    max_force = nebtools.get_fmax()
    print("Maximum force: %.2f eV/A" % max_force)

    # Create a figure like that coming from ASE-GUI
    fig = nebtools.plot_band()
    n_plot(None, None, xs=14, ys=14)
    fig.savefig('NEB_CC_Start_to_End.pdf')
    plt.show()

exit()
