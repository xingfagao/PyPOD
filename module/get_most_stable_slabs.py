"""get_most_stable_slabs.py"""
import glob,os
def get_most_stable_slabs(filenames):
    slab_structure = []
    slab_mi_index = []
    for i in range(len(filenames)):
        filename = filenames[i].split('.')[0]
        slab_info = filename.split('_')
        mp_id = slab_info[0]
        alloy_name = slab_info[1]
        mi_index = slab_info[2]
        slab_num = slab_info[3]
        slab_energy = os.popen('grep \'energy  without\' ' +
                               filename + '/OUTCAR').read().split()
        print(filename + '    grep \'energy  without\' ' +
              filename + '/OUTCAR')
        slab_energy = slab_energy[-1]
        slab_tmp = (filename,float(slab_energy),mp_id,alloy_name,mi_index,slab_num)
        slab_structure.append(slab_tmp)
        if mi_index not in slab_mi_index:
            slab_mi_index.append(mi_index)
    fn = open('slab_energy.log','w')
    for item in slab_structure:
        fn.write(item[0] + '    ' + str(item[1]) + '\n')
    most_stable_slabs = []
    for mi_index in slab_mi_index:
        energy_ref = 1E10
        for item_2 in slab_structure:
            if item_2[4] == mi_index:
                energy_tmp = item_2[1]
                if float(energy_tmp) < float(energy_ref):
                    energy_ref = energy_tmp
                    most_stable_slab_tmp = item_2
        most_stable_slabs.append(most_stable_slab_tmp)
    fn.write("\n\n\nThe most stable slabs are:\n")
    for stable_slab in most_stable_slabs:
        fn.write(stable_slab[0]+'    '+str(stable_slab[1])+'\n')
    return most_stable_slabs
