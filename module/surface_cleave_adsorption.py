"""surface_cleave_adsorption.py"""
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, \
    Structure, Lattice, ReconstructionGenerator, center_slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints
from pymatgen.core.structure import Structure, Molecule
import glob, os
import numpy as np
from pymatgen.ext.matproj import MPRester
USER_API_KEY="YOU_API_KEY"
def my_write_gjf(fname, atoms, adsorbate, ads_struc):
    atom_z = []
    for site_i in range(len(ads_struc)-len(adsorbate)):
        atom_z.append(ads_struc[site_i].c)
    atom_z.sort()
    z_value_fix = 0.5*atom_z[0] + 0.5*atom_z[-1]
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom_i in range(len(atoms)):
        atom = atoms[atom_i]
        if ads_struc[atom_i].c > z_value_fix:
            f.write('%s %s %f %f %f\n'%(atom.symbol,
                                        str(0),
                                        atom.position[0],
                                        atom.position[1],
                                        atom.position[2]))
        else:
            f.write('%s %s %f %f %f\n'%(atom.symbol,
                                        str(-1),
                                        atom.position[0],
                                        atom.position[1],
                                        atom.position[2]))
    cell = atoms.cell
    for i in cell:
        f.write('Tv ')
        for j in i:
            f.write('%f '%j)
        f.write('\n')
    f.write('\n')
    f.close()
def my_write_gjf2(fname, atoms, adsorbate, ads_struc):
    atom_z = []
    for site_i in range(len(ads_struc)-len(adsorbate)):
        atom_z.append(ads_struc[site_i].c)
    atom_z.sort()
    z_value_fix = 0.5*atom_z[0] + 0.5*atom_z[-1]
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom_i in range(len(atoms)-len(adsorbate)):
        atom = atoms[atom_i]
        if ads_struc[atom_i].c > z_value_fix:
            f.write('%s %s %f %f %f\n'%(atom.symbol,
                                        str(0),
                                        atom.position[0],
                                        atom.position[1],
                                        atom.position[2]))
        else:
            f.write('%s %s %f %f %f\n'%(atom.symbol,
                                        str(-1),
                                        atom.position[0],
                                        atom.position[1],
                                        atom.position[2]))
    cell = atoms.cell
    for i in cell:
        f.write('Tv ')
        for j in i:
            f.write('%f '%j)
        f.write('\n')
    f.write('\n')
    f.close()
def write_gjf_bulk(fname, atoms):
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom_i in range(len(atoms)):
        atom = atoms[atom_i]
        f.write('%s %s %f %f %f\n'%(atom.symbol,
                                    str(0),
                                    atom.position[0],
                                    atom.position[1],
                                    atom.position[2]))
    cell = atoms.cell
    for i in cell:
        f.write('Tv ')
        for j in i:
            f.write('%f '%j)
        f.write('\n')
    f.write('\n')
    f.close()
def is_slab_symmetrical(slab):
    atom_list = slab.species
    coords_old = slab.frac_coords
    slab_layer = []
    for i in range(len(coords_old)):
        if "%.8f" % (coords_old[i][2]) not in slab_layer:
            slab_layer.append("%.8f" % (coords_old[i][2]))
    slab_layer.sort()
    number_of_slab_layers = len(slab_layer)
    for j in range(int(number_of_slab_layers / 2)):
        slab_layer_i = []
        slab_layer_j = []
        for k in range(len(coords_old)):
            if "%.8f" % (coords_old[k][2]) == slab_layer[j]:
                slab_layer_i.append(atom_list[k])
            elif "%.8f" % (coords_old[k][2]) == slab_layer[-1 - j]:
                slab_layer_j.append(atom_list[k])
            else:
                continue
        slab_layer_i.sort()
        slab_layer_j.sort()
        if slab_layer_i == slab_layer_j:
            continue
        else:
            return False
    return True
def slab_flipping(slab):
    Lattice_new = Lattice.from_parameters(a=slab.lattice.a,
                                          b=slab.lattice.b,
                                          c=slab.lattice.c,
                                          alpha=180 - slab.lattice.alpha,
                                          beta=180 - slab.lattice.beta,
                                          gamma=slab.lattice.gamma)
    flip_array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    atom_list = slab.species
    coords_old = slab.frac_coords
    coords_new = np.dot(coords_old, flip_array)
    slab_new = Structure(Lattice_new, atom_list, coords_new)
    return slab_new
def gen_gjf_from_mp_id(mp_id):
    with MPRester(USER_API_KEY) as m:
        doc = m.summary.get_data_by_id(mp_id)
        structure = doc.structure
        print(structure)
        print(doc.formula_pretty)
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
        atoms = AseAtomsAdaptor.get_atoms(structure)
        fname = mp_id + '_' + doc.formula_pretty + '.gjf'
        write_gjf_bulk(fname,atoms)
    return atoms
H = Molecule("H", [[0, 0, 0]])
O = Molecule("O", [[0, 0, 0]])
OH = Molecule("OH", [[0, 0, 0], [-0.793, 0.384, 0.422]])
OOH = Molecule("OOH", [[0, 0, 0], [-1.067, -0.403, 0.796], [-0.696, -0.272, 1.706]])
def generate_adsorption_structures_user_defined(
    slab,
    molecule,
    site_type,
    repeat=None,
    min_lw=5.0,
    translate=True,
    reorient=True,
    find_args=None,
    ):
    if repeat is None:
        xrep = np.ceil(min_lw / np.linalg.norm(slab.lattice.matrix[0]))
        yrep = np.ceil(min_lw / np.linalg.norm(slab.lattice.matrix[1]))
        repeat = [xrep, yrep, 1]
    structs = []
    asf = AdsorbateSiteFinder(slab)
    find_args = find_args or {}
    for coords in asf.find_adsorption_sites(**find_args)[site_type]:
        structs.append(
            asf.add_adsorbate(
                molecule,
                coords,
                repeat=repeat,
                translate=translate,
                reorient=reorient,
            )
        )
    return structs
def gen_adsorbate_on_struct_site(slab,adsorbate, slab_name,site_type):
    ads_structs_names = []
    slab = center_slab(slab)
    sg = SpacegroupAnalyzer(slab, symprec=0.01, angle_tolerance=5.0)
    cs = sg.get_crystal_system()
    ads_structs = generate_adsorption_structures_user_defined(slab,
                                                              adsorbate,
                                                              site_type,
                                                              min_lw=6.0)
    num = len(ads_structs)
    print('%i adducts in total.'%num)
    for i in range(num):
        atoms = AseAtomsAdaptor.get_atoms(ads_structs[i])
        fname = slab_name + '_%i.gjf'%(i)
        my_write_gjf(fname, atoms, adsorbate, ads_structs[i])
        my_write_gjf2(slab_name+'.gjf', atoms, adsorbate, ads_structs[i])
        if fname not in ads_structs_names:
            ads_structs_names.append(fname)
        if slab_name + '.gjf' not in ads_structs_names:
            ads_structs_names.append(slab_name + '.gjf')
    return ads_structs_names
def slab_io(slabs,alloy_name):
    mi_string_old = 'old'
    slab_count = 1
    for slab in slabs:
        mi_string = "".join([str(i) for i in slab.miller_index])
        if mi_string == mi_string_old:
            slab_count += 1
            slab_name = alloy_name + '_' + mi_string + '_' + '%i' % (slab_count)
        else:
            mi_string_old = mi_string
            slab_count = 1
            slab_name = alloy_name + '_' + mi_string + '_' + '%i' % (slab_count)
        if is_slab_symmetrical(slab):
            center_slab(slab)
            fname1 = slab_name+"u.cif"
            CifWriter(slab).write_file(fname1)
            atoms = AseAtomsAdaptor.get_atoms(slab)
            adsorbate = []
            fname2 = slab_name + "u.gjf"
            my_write_gjf(fname2, atoms, adsorbate, slab)
        else:
            slab_new = slab_flipping(slab)
            center_slab(slab_new)
            center_slab(slab)
            fname1a = slab_name + "u.cif"
            fname1b = slab_name + "d.cif"
            CifWriter(slab).write_file(fname1a)
            CifWriter(slab_new).write_file(fname1b)
            atoms_u = AseAtomsAdaptor.get_atoms(slab)
            atoms_d = AseAtomsAdaptor.get_atoms(slab_new)
            adsorbate = []
            fname2a = slab_name + "u.gjf"
            fname2b = slab_name + "d.gjf"
            my_write_gjf(fname2a, atoms_u, adsorbate, slab)
            my_write_gjf(fname2b, atoms_d, adsorbate, slab_new)
    return
