# /// script
# dependencies =[
# "ase==3.26.0",
# "pytest",
# "numpy",
# ]
# ///

import aneto_utils as au
import numpy as np
import os

#@pytest.fixture
#def
#def test_read_emodulus_data():
#    test_file_name = 'Emod.outcar'
#    test_modulus = au.read_emodulus_data(test_file_name)
#    ref_modulus_kB =np.array([[1767.2697, 649.8595, 876.9918, 1.7076,  -0.0844, -0.2260],
#                              [649.8595, 1767.0907, 876.7797, 1.4175,-0.4230, 0.2622],
#                              [876.9918, 876.7797, 1488.2037, -2.6521, -0.3766, -0.0616],
#                              [1.7076, 1.4175, -2.6521, 368.5446, -0.4371, -0.4127],
#                              [-0.0844, -0.4230, -0.3766, -0.4371, 651.1969, -0.1194],
#                              [-0.2260, 0.2622, -0.0616, -0.4127, -0.1194, 651.0030]
#                              ]
#            )
#    ref_modulus_GPa = 0.1*ref_modulus_kB
#    assert test_modulus.all == ref_modulus_GPa.all
#
def test_read_emodulus_file():
#    test_file_name = 'Emod.outcar'
    test_modulus = au.read_vasp_cij('./test_files')
    ref_modulus_kB =np.array([[1767.2697, 649.8595, 876.9918, 1.7076,  -0.0844, -0.2260],
                              [649.8595, 1767.0907, 876.7797, 1.4175,-0.4230, 0.2622],
                              [876.9918, 876.7797, 1488.2037, -2.6521, -0.3766, -0.0616],
                              [1.7076, 1.4175, -2.6521, 368.5446, -0.4371, -0.4127],
                              [-0.0844, -0.4230, -0.3766, -0.4371, 651.1969, -0.1194],
                              [-0.2260, 0.2622, -0.0616, -0.4127, -0.1194, 651.0030]
                              ]
            )
#    ref_modulus_kB =np.array([[1767.2697, 649.8595, 876.9918, 0, 0, 0],
#                              [649.8595, 1767.0907, 876.7797, 0, 0, 0],
#                              [876.9918, 876.7797, 1488.2037, 0, 0, 0],
#                              [ 0, 0, 0, 368.5446, 0, 0],
#                              [ 0, 0, 0, 0, 651.1969, 0],
#                              [ 0, 0, 0, 0, 0, 651.0030]
#                              ]
#            )
    ref_modulus_GPa = 0.1*ref_modulus_kB
    assert (test_modulus == ref_modulus_GPa).all

def test_get_vasp_stress():
    os.system('cp ./test_files/OUTCAR.defect ./test_files/OUTCAR')
    au.get_vasp_stress('./test_files')
    assert os.path.exists('./test_files/residual_stress.outcar')

def test_read_stress_data():
    test_stress = au.read_stress_data('./test_files/residual_stress.outcar')
    ref_stress = np.array([[ 6.48601, -0.00004, -0.00020],
                           [-0.00004, 6.48600, -0.00028],
                           [ -0.00020, -0.00028, 15.53205],
                          ]
            )
    ref_stress_GPa = 0.1*ref_stress
    assert (test_stress ==ref_stress_GPa).all

def test_write_aneto_input():


    test_residual_stress=np.zeros((3,3))
    test_residual_stress[0,0]=-1.239
    test_residual_stress[1,1]=-1.239
    test_residual_stress[2,2]=-0.837

    test_elastic_testor=np.zeros((6,6))
    test_elastic_testor[0,0] = 140.28
    test_elastic_testor[0,1] = 69.72
    test_elastic_testor[0,2] = 65.11
    test_elastic_testor[1,0] = 69.72
    test_elastic_testor[1,1] = 140.28
    test_elastic_testor[1,2] = 65.11
    test_elastic_testor[2,0] = 65.11
    test_elastic_testor[2,1] = 65.11
    test_elastic_testor[2,2] = 168.00
    test_elastic_testor[3,3] = 25.96
    test_elastic_testor[4,4] = 25.96
    test_elastic_testor[5,5] = 35.28

    test_cell=np.zeros((3,3))
    test_cell[0,0] = 4.0
    test_cell[1,0] = -2.0
    test_cell[1,1] = 3.464102
    test_cell[0,2] = 4.803

    test_alat=3.230
    
    test_filename='intput_elast_test'

    au.write_aneto_input(test_residual_stress, test_elastic_testor, test_cell, test_alat, test_filename)

    test_file=open(test_filename, 'r')
    test_line = test_file.readline().strip()
    ref_filename = './test_files/input_elast_eps0'
    ref_file=open(ref_filename, 'r')
    ref_line = ref_file.readline().strip()
    while (len(ref_line)>0):
        if (ref_line!=test_line):
            print(ref_line)
            print(test_line)
            assert False
        ref_line=ref_file.readline().strip()
        test_line = test_file.readline().strip()
    assert True







