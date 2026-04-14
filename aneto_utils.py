# /// script
# dependencies = [
#    "pymatgen",
#    "ase = 3.26.0",
#    "json",
#    "glob",
# ]
# ///

from ase.io import write
from pymatgen.io.vasp.outputs import OUTCAR
import json
import glob

def read_input_args():

def vasp_outcar(defect_dir):
    #   INPUTs
    # defect_dir is string to the path
    #
    #   RETURNs
    # pymatgen outcar object
    #

    outcar = Outcar("{}OUTCAR".format(defect_dir))
    return outcar

def get_vasp_stress():
    #   INPUTs
    # defect_dir is string to the path
    #
    #   RETURNs
    # pymatgen outcar object
    #
    os.system("grep 'FORCE on cell' OUTCAR -A 30 >force.outcar")
    outcar = Outcar("{}OUTCAR".format(defect_dir))
    stress_tensors = outcar.stress

def read_vasp_cij(elast_dir):
    if exits('cij.json')
       cij = read_cij_json
    else
       if not(exists()):
           outcar_path = "{}/OUTCAR".format(elast_dir)
           os.system("grep 'TOTAL ELASTIC MODULI' {} -A 8 >Emod.outcar".format(outcar))
       cij = read_Emod_file
       write_cij_json(cij, 'cij.json')
    return cij

def read_outcar_energy(path):
    line = os.system( "grep 'TOTEN' OUTCAR").strip().split()
    energy = float(line[-2])
    return energy

def write_aneto_input(defect_stress, elast_tensor, cell):


def read_config_file():


