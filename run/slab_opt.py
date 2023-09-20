"""slab_opt.py"""
import glob,os,math
from sys import path
path.append('../module')
def write_pbs(filename,jobpath):
    fn = open(filename + '.sh','w')
    fn.write("#!/bin/sh\n\n")
    fn.write("#SBATCH -N 1\n#SBATCH -n 64\n#SBATCH -t 14-00:00\n#SBATCH -J a"
             "\n#SBATCH -p wzhcnormal\n#SBATCH -o out.%jt\n#SBATCH -e err.%j\n\n\n")
    fn.write("module purge\n\n")
    fn.write("module load compiler/intel/2017.5.239\n")
    fn.write("module load mpi/intelmpi/2017.4.239\n\n")
    fn.write("CDIR=$PWD\n")
    fn.write("joblist=$(cat $CDIR/"+filename+")\n")
    fn.write("for jobi in $joblist\ndo\n")
    fn.write("    cd $jobi\n")
    fn.write("    for (( i=1; i<=4; i++ ));\n")
    fn.write("    do\n")
    fn.write("        cp INCAR_$i INCAR\n")
    fn.write("        cp CONTCAR POSCAR\n")
    fn.write("        mpirun -np 64 /work/share/ac0i6ytcwi/apps/"
             "vasp544/vasp.5.4.4/bin/vasp_std > log_$i\n")
    fn.write("    done\n")
    fn.write("    rm CHG CHGCAR DOSCAR DYNMAT EIGENVAL IBZKPT PCDAT"
             " PROCAR REPORT vasprun.xml WAVECAR XDATCAR\n")
    fn.write("done\n")
    fn.close()
    return
current_dir = os.getcwd()
fn2 = open('BinaryAlloy_selected.txt','r')
Alloy_count = 1
numofalloy = 1
os.system("rm *joblist*")
for line in fn2:
    line = line.split()
    mp_id = line[0]
    formula_pretty = line[1]
    is_magnetic = line[2]
    filename_local = mp_id + '_' + formula_pretty
    target_dir = current_dir + '/' + filename_local + '/slab'
    joblist_name = 'Alloy_slab_joblist'+str(Alloy_count)
    if os.path.exists(joblist_name):
        fn1 = open(joblist_name,'a')
    else:
        fn1 = open(joblist_name,'w')
    os.chdir(target_dir)
    filenames = glob.glob('*.cif')
    for i in range(len(filenames)):
        filename = filenames[i].split('.')[0]
        if os.path.exists(filename):
            os.system('mv ' + filename + '/' + filename +'.gjf ./' )
        cell_para = os.popen('grep Tv ' + filename +'.gjf').read().split()
        cell_a = math.sqrt(float(cell_para[1])**2+
                           float(cell_para[2])**2+
                           float(cell_para[3])**2)
        cell_b = math.sqrt(float(cell_para[5])**2+
                           float(cell_para[6])**2+
                           float(cell_para[7])**2)
        cell_lengths = [cell_a,cell_b]
        fn_k = open('KPOINTS_tmp','w')
        K_POINTS = []
        for cell_length in cell_lengths:
            K_x_a = int(45/float(cell_length))
            if (K_x_a % 2) == 0:
                K_x = K_x_a + 1
            else:
                K_x = K_x_a
            K_POINTS.append(K_x)
        K_POINTS.append(1)
        fn_k.write('Automatic generation\n0\nMonhkorst-Pack\n')
        for K_i in K_POINTS:
            fn_k.write(str(K_i)+'  ')
        fn_k.write('\n0.0 0.0 0.0')
        fn_k.close()
        os.system('gjf2vasp ' + filename + '.gjf 0')
        fn1.write(target_dir + '/' + filename + '\n')
        os.system('cp ' + current_dir + '/INCAR_1 ' +
                  target_dir + '/' + filename + '/INCAR_1')
        os.system('cp ' + current_dir + '/INCAR_2 ' +
                  target_dir + '/' + filename + '/INCAR_2')
        os.system('cp ' + current_dir + '/INCAR_3 ' +
                  target_dir + '/' + filename + '/INCAR_3')
        os.system('cp ' + current_dir + '/INCAR_4 ' +
                  target_dir + '/' + filename + '/INCAR_4')
        os.system('mv KPOINTS_tmp ' + target_dir + '/' +
                  filename + '/KPOINTS')
    fn1.close()
    os.chdir(current_dir)
    write_pbs(joblist_name,current_dir)
    os.system("chmod u+x " + joblist_name + '.sh')    
    Alloy_count += 1
fn2.close()
