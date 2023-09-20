"""generate_adsorption_strucs.py"""
import glob,os
from surface_cleave_adsorption import *
from get_most_stable_slabs import * 
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, \
    Structure, Lattice, ReconstructionGenerator, center_slab
def gen_ads_structs(stable_slabs,adsorbate_name):
    slab_adsorbate_name = []
    for item in stable_slabs:
        slab_name = item[0]
        slab_struct = Structure.from_file(slab_name + '.cif')
        if is_slab_symmetrical(slab_struct):
            asf1 = AdsorbateSiteFinder(slab_struct)
            ads_sites1 = asf1.find_adsorption_sites()
            site_types = ads_sites1.keys()
            for site_type in site_types:
                slab_adsorbate_name_tmp = gen_adsorbate_on_struct_site(slab_struct,
                                                                       adsorbate_name,
                                                                       slab_name,
                                                                       site_type)
                for item in slab_adsorbate_name_tmp:
                    if item not in slab_adsorbate_name:
                        slab_adsorbate_name.append(item)
        else:
            asf1 = AdsorbateSiteFinder(slab_struct)
            ads_sites1 = asf1.find_adsorption_sites()
            site_types = ads_sites1.keys()
            for site_type1 in site_types:
                slab_adsorbate_name_tmp1 = gen_adsorbate_on_struct_site(slab_struct,
                                                                        adsorbate_name,
                                                                        slab_name,
                                                                        site_type1)
                for item1 in slab_adsorbate_name_tmp1:
                    if item1 not in slab_adsorbate_name:
                        slab_adsorbate_name.append(item1)
            slab_struct_new = slab_flipping(slab_struct)
            if slab_name[-1]=='u':
                slab_name_new = slab_name[:-1] + 'd'
            else:
                slab_name_new = slab_name[:-1] + 'u'
            slab_struct_new = Structure.from_file(slab_name_new + '.cif')
            asf2 = AdsorbateSiteFinder(slab_struct_new)
            ads_sites2 = asf2.find_adsorption_sites()
            site_types2 = ads_sites2.keys()
            for site_type2 in site_types:
                slab_adsorbate_name_tmp2 = gen_adsorbate_on_struct_site(slab_struct_new,
                                                                        adsorbate_name,
                                                                        slab_name_new,
                                                                        site_type2)
                for item2 in slab_adsorbate_name_tmp2:
                    if item2 not in slab_adsorbate_name:
                        slab_adsorbate_name.append(item2)
        os.system('mv *gjf ../')
    return slab_adsorbate_name
