"""get_Eads.py"""
import glob,os
def get_Eads(E_adsorbate):
    filenames = glob.glob('*.txt')
    current_dir = os.getcwd()
    slab_structure = []
    slab_Eads = []
    for i in range(len(filenames)):
        filename = filenames[i].split('.')[0]
        filename = filename.split('_')
        mp_id = filename[0]
        alloy_name = filename[1]
        mi_index = filename[2]
        slab_num = filename[3]
        if len(filename) == 5:
            ads_num = filename[4]
        slab_name = mp_id + '_' + alloy_name + '_' + mi_index + '_' + slab_num
        if os.path.exists(slab_name+'/OSZICAR') and len(filename) == 4:
            slab_energy = os.popen('grep E0 ' + slab_name + '/OSZICAR').read().split()
            if len(slab_energy) > 6:
                slab_energy = slab_energy[-6]
                slab_tmp = (slab_name,float(slab_energy))
                slab_structure.append(slab_tmp)
    for i in range(len(slab_structure)):
        slab_tmp_ads = []
        for j in range(len(filenames)):
            filename = filenames[j].split('.')[0]
            filename2 = filename.split('_')
            mp_id = filename2[0]
            alloy_name = filename2[1]
            mi_index = filename2[2]
            slab_num = filename2[3]
            slab_name = mp_id + '_' + alloy_name + '_' + mi_index + '_' + slab_num
            if slab_name == slab_structure[i][0] and len(filename2) == 5:
                if os.path.exists(filename+'/OSZICAR'):
                    slab_ads_energy = os.popen('grep E0 ' + filename +
                                               '/OSZICAR').read().split()
                    if len(slab_ads_energy)>6:
                        slab_ads_energy = slab_ads_energy[-6]
                        slab_tmp = (float(slab_ads_energy),filename)
                        slab_tmp_ads.append(slab_tmp)
        slab_tmp_ads.sort()
        if len(slab_tmp_ads) > 0:
            slab_Eads_tmp = (slab_structure[i][0], slab_structure[i][1],
                             slab_tmp_ads[0][1], slab_tmp_ads[0][0])
            slab_Eads.append(slab_Eads_tmp)
    print(slab_Eads)
    fn = open('ads_energy.log', 'w')
    for item in slab_Eads:
        fn.write(item[0] + '    ' + str(item[1]) + '   ' + item[2] + '    '
                 + str(item[3]) + '   ' + str(item[3]-item[1]-E_adsorbate) + '\n')
    fn.close()
    return slab_Eads
if __name__=='__main__':
    current_dir = os.getcwd()
    E_OH = float(-.77133971E+01)
    fn1 = open('BinaryAlloy_selected.txt', 'r')
    fn2 = open('Eads_OH.txt', 'w')
    fn2.write('%-25s%-15s%-25s%15s%8s\n'%('slab_name',
                                          'slab_energy',
                                          'stable_ads_str',
                                          'stable_ads_energy',
                                          'Eads'))
    for line in fn1:
        line = line.split()
        mp_id = line[0]
        formula_pretty = line[1]
        target_dir = current_dir + '/' + mp_id + '_' + formula_pretty
        if os.path.exists(target_dir):
            os.chdir(target_dir)
            slab_Eads_tmp=get_Eads(E_OH)
            os.chdir(current_dir)
            for item in slab_Eads_tmp:
                fn2.write('%-25s%-15s%-25s%15s%10.2f\n'%(item[0],
                                                         str(item[1]),
                                                         item[2],
                                                         str(item[3]),
                                                         item[3]-item[1]-E_OH))
    fn2.close()
    fn1.close() 
