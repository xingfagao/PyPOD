"""slabads_buildopt.py"""
import glob, os, math
from sys import path
path.append('../module')
from surface_cleave_adsorption import *
from get_most_stable_slabs import get_most_stable_slabs
from generate_adsorption_strucs import gen_ads_structs
def write_pbs(filename,jobpath):
    fn = open(filename + '.sh', 'w')
    fn.write("#!/bin/sh\n\n")
    fn.write("#SBATCH -N 1\n#SBATCH -n 64\n#SBATCH -t 7-00:00\n#SBATCH -J a\n")
    fn.write("#SBATCH -p wzhcnormal\n#SBATCH -o out.%jt\n#SBATCH -e err.%j\n\n\n")
    fn.write("module purge\n\n")
    fn.write("module load compiler/intel/2017.5.239\n")
    fn.write("module load mpi/intelmpi/2017.4.239\n\n")
    fn.write("CDIR=$PWD\n")
    fn.write("joblist=$(cat $CDIR/"+filename+")\n")
    fn.write("for jobi in $joblist\ndo\n")
    fn.write("    cd $jobi\n")
    fn.write("    for (( i=1; i<=3; i++ ));\n")
    fn.write("    do\n")
    fn.write("        cp INCAR_$i INCAR\n")
    fn.write("        cp CONTCAR POSCAR\n")
    fn.write("        mpirun -np 64 /YOUR_PATH_TO_VASP/vasp_std > log_$i\n")
    fn.write("    done\n")
    fn.write("    rm CHG CHGCAR DOSCAR DYNMAT EIGENVAL PCDAT PROCAR")
    fn.write(" REPORT vasprun.xml WAVECAR\n")
    fn.write("done\n")
    fn.close()
    return
current_dir = os.getcwd()
fn2 = open('BinaryAlloy_selected.txt','r')
Alloy_count = 0
numofalloy = 1
os.system("rm *ads*joblist*")
for line in fn2:
    line = line.split()
    mp_id = line[0]
    formula_pretty = line[1]
    is_magnetic = line[2]
    filename_local = mp_id + '_' + formula_pretty
    target_dir1 = current_dir + '/' + filename_local + '/slab'
    os.chdir(target_dir1)
    slab_names = glob.glob('*.txt')
    stable_slabs = get_most_stable_slabs(slab_names)
    slab_adsorbate_names = gen_ads_structs(stable_slabs,OH)
    os.chdir(current_dir)
    ads_i = int(Alloy_count/numofalloy)
    joblist_name = 'Alloy_slabads_joblist'+str(ads_i)
    if os.path.exists(joblist_name):
        fn1 = open(joblist_name,'a')
    else:
        fn1 = open(joblist_name,'w')
    target_dir2 = current_dir + '/' + filename_local
    os.chdir(target_dir2)
    for i in range(len(slab_adsorbate_names)):
        filename = slab_adsorbate_names[i].split('.')[0]
        if os.path.exists(filename):
            os.system('mv ' + filename + '/' + filename +'.gjf ./')
        os.system('gjf2vasp' + filename + '.gjf 0')
        fn1.write(target_dir2 + '/' + filename + '\n')
        os.system('cp ' + current_dir + '/INCAR_1 ' +
                  target_dir2 + '/' + filename + '/INCAR_1')
        os.system('cp ' + current_dir + '/INCAR_2 ' +
                  target_dir2 + '/' + filename + '/INCAR_2')
        os.system('cp ' + current_dir + '/INCAR_3 ' +
                  target_dir2 + '/' + filename + '/INCAR_3')
        os.system('cp ' + current_dir +  '/KPOINTS_1 ' +
                  target_dir2 + '/' + filename + '/KPOINTS')
    fn1.close()
    os.chdir(current_dir)
    write_pbs(joblist_name,current_dir)
    os.system("chmod u+x " + joblist_name + '.sh')    
    Alloy_count += 1
fn2.close()
