# /// script
# dependencies = [
#    "pymatgen",
#    "ase ==3.26.0",
#    "numpy",
#    "json",
#    "glob",
# ]
# ///

import ase
from ase.io import write
from ase.io.vasp import read_vasp
import json
import glob
import os
import numpy as np


def read_emodulus_data(filename):
    #########
    # INPUT
    # filename: str, location of file with OUTCAR stress data
    #
    # RETURN
    # stress: np.array(3x3) with symmetric stress tensor#
    ###########

    fil = open(filename)
    lnum = 0
    line = []
    while lnum < 3:
        line = fil.readline()
        lnum = lnum + 1
    e_mod = np.zeros((6, 6))
    for i in range(6):
        line = fil.readline().strip().split()
        for j in range(6):
            e_mod[i, j] = float(line[j + 1])
    e_mod = e_mod * 0.1  # convert to GPa from kBar
    print(e_mod)
    return e_mod


def read_stress_data(filename):
    #########
    # INPUT
    # filename: str, location of file with OUTCAR stress data
    #
    # RETURN
    # stress: np.array(3x3) with symmetric stress tensor#
    ###########

    fil = open(filename)
    lnum = 0
    line = []
    while lnum < 14:
        line = fil.readline()
        lnum = lnum + 1
    line = fil.readline().strip().split()
    stress = np.zeros((3, 3))
    stress[0, 0] = float(line[2])
    stress[1, 1] = float(line[3])
    stress[2, 2] = float(line[4])
    stress[0, 1] = float(line[5])
    stress[1, 0] = stress[0, 1]
    stress[1, 2] = float(line[6])
    stress[2, 1] = stress[1, 2]
    stress[0, 2] = float(line[7])
    stress[2, 0] = stress[0, 2]
    stress = stress * 0.1  # convert to GPa from kBar

    return stress


def get_vasp_stress(data_dir="./"):
    # INPUTs
    # defect_dir is string to the path
    #
    #   RETURNs
    # pymatgen outcar object
    #
    stress_data_file = "{}/residual_stress.outcar".format(data_dir)
    if not (os.path.exists(stress_data_file)):
        os.system(
            "grep 'FORCE on cell' {}/OUTCAR -A 30 >{}".format(
                data_dir, stress_data_file
            )
        )
    residual_stress = read_stress_data(stress_data_file)
    return residual_stress


def write_cij_json(cij, filename="cij.json"):
    cij_dict = {"cij": cij.tolist()}
    with open(filename, "w") as fil:
        json.dump(cij_dict, fil, indent=4)


def read_cij_json(filename):
    cij_dict = json.load(open(filename, "r"))
    cij = np.array(cij_dict["cij"])
    return cij


def read_cij(elast_dir):
    if os.path.exists("cij.json"):
        cij = read_cij_json
    else:
        emod_filename = "{}/emodulus_data.outcar".format(elast_dir)
        if not (os.path.exists(emod_filename)):
            outcar_path = "{}/OUTCAR.elastic".format(elast_dir)
            os.system(
                "grep 'TOTAL ELASTIC MODULI' {} -A 8 >{}".format(
                    outcar_path, emod_filename
                )
            )
        cij = read_emodulus_data(emod_filename)
        write_cij_json(cij, "cij.json")
    return cij


def read_outcar_energy(path):
    line = os.system("grep 'ENERGIE' {} -A 5 >temp.energy".format(path))
    fl = open("temp.energy", "r")
    for i in range(3):
        line = fl.readline()
    line = line.strip().split()
    print(line)
    energy = float(line[-2])
    return energy


def write_aneto_input(
    residual_stress, elast_tensor, cell, alat, filename="input_elast"
):
    #############
    #   generates input scrpt in NAETO format given the
    #   necessary input data
    #   by default names file input elast, but this will
    #   overwrite other files with that name
    #############
    fil = open(filename, "w")
    fil.write(" &input\n")
    # residual stress which is the stres due to the dipole in the unit cell
    # not dipole itself. in GPa (not kB in VASP)
    for i in range(3):
        for j in range(3):
            fil.write(
                "   CVoigt({},{})={:.2f}\n".format(i + 1, j + 1, elast_tensor[i, j])
            )
    for i in [4, 5, 6]:
        fil.write("   CVoigt({},{})={:.2f}\n".format(i, i, elast_tensor[i - 1, i - 1]))
    fil.write("\n")
    # residual stress which is the stres due to the dipole in the unit cell
    # not dipole itself. in GPa (not kB in VASP)
    for i in range(3):
        for j in range(3):
            fil.write(
                "   sigma_res({},{})={}\n".format(i + 1, j + 1, residual_stress[i, j])
            )
    fil.write("\n")
    # lattice vectors in fractional coordinates( where 1= alat)
    # not the same as normal fractional coordinates since may not be cubic
    for i in range(3):
        for j in range(3):
            fil.write("   A{}_ref({})={:.6f}\n".format(i + 1, j + 1, cell[i, j]))
    # a lattice vector of structure
    fil.write("   alat={:.3f}\n".format(alat))
    fil.write(" &end\n")
    fil.close()


def save_reference_data(scaled_cell, alat, total_energy, json_name="ref_data.json"):
    ##########
    # Saves reference data to a json
    #
    # INPUTS:
    # scaled_cell: 3x3 numpy array or nested 3x3 list
    # alat: float
    # total_energy: float
    # json_name: str, default 'ref_data.json'
    #
    ##########

    ref_data = {
        "scaled_cell": scaled_cell.tolist(),
        "alat": alat,
        "total_energy": total_energy,
    }
    with open(json_name, "w") as fl:
        json.dump(ref_data, fl, indent=4)


def get_prim_alat(path_to_poscar):

    atoms = read_vasp(path_to_poscar)
    cell = np.array(atoms.cell)
    alat = np.linalg.norm(cell[0])
    return alat


def get_reference_cell_data(path_to_poscar):
    #############
    # INPUT
    # path_to_poscar: str, where POSCAR of reference
    #    supercell can be found in vasp POSCAR format
    #
    # RETURNS
    # cell: 3x3 numpy array or list of list
    #
    #############

    atoms = read_vasp(path_to_poscar)
    cell = np.array(atoms.cell)
    return cell


def scale_cell(cell, alat):
    ###########
    # INPUT
    # cell: 3x3 np array or list of lists
    # alat: float
    #
    # RETURN:
    # scaled_cell: 3x3 np array or list of lists
    ###########

    scaled_cell = cell / alat
    return scaled_cell


def load_ref_data(ref_data_string):
    if ref_data_string[-4:] != "json":
        try:
            os.path.exists("{}/OUTCAR.super".format(ref_data_string))
        except:
            print("No OUTCAR.super found. Please provide path to json\
                \nor dir with OUTCAR.super, POSCAR.prim, and POSCAR.super")
        cell_super = get_reference_cell_data("{}/POSCAR.super".format(ref_data_string))
        alat = get_prim_alat("{}/POSCAR.prim".format(ref_data_string))
        scaled_cell = scale_cell(cell_super, alat)

        ref_energy = read_outcar_energy("{}/OUTCAR.super".format(ref_data_string))
        save_reference_data(scaled_cell, alat, ref_energy)

    with open('ref_data.json', "r") as fl:
        ref_data = json.load(fl)
    return ref_data


def setup_aneto(
    defect_dir,
    ref_dict,
    e_mod=[],
    save_dir="aneto_input_files",
    label="test",
    elastic_dir="./",
):

    #    ref_data = load_ref_data(ref_data_string)
    defect_stress = get_vasp_stress(defect_dir)
    if len(e_mod) == 0:
        e_mod = read_cij(elastic_dir)
    if not (os.path.exists(save_dir)):
        os.mkdir(save_dir)
    write_aneto_input(
        defect_stress,
        e_mod,
        np.array(ref_dict["scaled_cell"]),
        ref_dict["alat"],
        "{}/input_elast_{}".format(save_dir, label),
    )
