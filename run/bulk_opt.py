"""bulk_opt.py"""
import glob,os
from sys import path
path.append('../module')
from pymatgen.io.ase import AseAtomsAdaptor
from surface_cleave_adsorption import gen_gjf_from_mp_id
current_dir = os.getcwd()
fn1 = open('joblist','w')
fn2 = open('BinaryAlloy_selected.txt','r')
for line in fn2:
    line = line.split()
    mp_id = line[0]
    formula_pretty = line[1]
    is_magnetic = line[2]
    atoms = gen_gjf_from_mp_id(mp_id)
    cell_lengths = atoms.cell.lengths()
    fn_k = open('KPOINTS_tmp','w')
    K_POINTS = []
    for cell_length in cell_lengths:
        K_x_a = int(45/float(cell_length))
        if (K_x_a % 2) == 0:
            K_x = K_x_a + 1
        else:
            K_x = K_x_a
        if K_x == 1:
            K_x = 3
        K_POINTS.append(K_x)
    fn_k.write('Automatic generation\n0\nMonhkorst-Pack\n')
    for K_i in K_POINTS:
        fn_k.write(str(K_i)+'  ')
    fn_k.write('\n0.0 0.0 0.0')
    fn_k.close()
    filename_local = mp_id + '_' + formula_pretty
    os.system('gjf2vasp ' + filename_local + '.gjf 0')
    target_dir = current_dir + '/' + filename_local
    os.chdir(target_dir)
    if not os.path.exists('bulk'):
        os.mkdir('bulk')
    fn1.write(target_dir + '/bulk\n')
    os.system('cp ../INCAR_1 ' + target_dir + '/bulk/INCAR_1')
    os.system('cp ../INCAR_2 ' + target_dir + '/bulk/INCAR_2')
    os.system('cp ../INCAR_3 ' + target_dir + '/bulk/INCAR_3')
    os.system('mv INCAR KPO* PO* ./bulk/')
    os.system('mv ../KPOINTS_tmp ' + target_dir + '/bulk/KPOINTS')
    os.chdir(current_dir)
fn1.close()
fn2.close()
