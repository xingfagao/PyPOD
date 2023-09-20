"""slab_build.py"""
import glob,os
from sys import path
path.append('../module')
import surface_cleave_adsorption
from pymatgen.io.ase import AseAtomsAdaptor
from surface_cleave_adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, \
    Lattice, ReconstructionGenerator, center_slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
current_dir = os.getcwd()
fn1 = open('BinaryAlloy_selected.txt','r')
fn2 = open('bulk_opt_info.txt','w')
fn3 = open('bulk_opt_failed.txt','w')
for line in fn1:
    line1 = line.split()
    mp_id = line1[0]
    formula_pretty = line1[1]
    is_magnetic = line1[2]
    target_dir = current_dir + '/' + mp_id + '_' + formula_pretty
    if os.path.exists(target_dir + '/bulk/'):
        termination_info = os.popen('grep \'reached required\' '
                                    + target_dir+'/bulk/OUTCAR').read().split()
        if len(termination_info) > 0:
            is_normal_termination = True
        else:
            is_normal_termination = False
            fn2.write("The bulk structure of " + mp_id + formula_pretty +
                      " is not optimized successfully!!! Please try again!!!\n")
            fn3.write(line)
        if is_normal_termination is True:
            os.chdir(target_dir)
            if not os.path.exists('slab'):
                os.system('mkdir slab')
            poscar = Poscar.from_file(target_dir+'/bulk/CONTCAR')
            struct = poscar.structure
            struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
            slabs = generate_all_slabs(struct, max_index=1,
                                       min_slab_size=8.0,
                                       min_vacuum_size=15.0,
                                       center_slab=False,
                                       symmetrize=False)
            alloy_name = mp_id + '_' + formula_pretty
            slab_io(slabs,alloy_name)
            os.system('mv *cif ./slab/')
            os.system('mv *gjf ./slab/')
        os.chdir(current_dir)
    else:
        fn2.write("The bulk structure of " + mp_id + formula_pretty +
                  " is not optimized successfully!!! Please try again!!!\n")
        fn3.write(line)
fn1.close()
fn2.close()
fn3.close()
