"""search_cif.py"""
import os.path,os
from pymatgen.io.cif import CifWriter
from mp_api.client import MPRester
USER_API_KEY="JBkm9BAbwqy33jXXfIsRhz8C0R2v9gkI" 
elements_list = ["Au", "Cu"]
with MPRester(api_key=USER_API_KEY) as m:
    results_a = []
    results_b = []
    for element_i in range(len(elements_list)):
        for element_j in range(element_i+1,len(elements_list)):
            results = m.summary.search(chemsys=elements_list[element_i]+"-"
                                               +elements_list[element_j],
                                       fields=["formula_pretty",
                                                "material_id",
                                                "is_stable",
                                                "theoretical",
                                                "nelements",
                                                "is_magnetic",
                                                "nsites",
                                                "structure"])
            for i in range(len(results)):
                result_0 = results[i]
                if result_0.nelements == 2:
                    if result_0.is_stable == True or result_0.theoretical == False:
                        results_a.append(results[i])
    fn = open("BinaryAlloy_all.txt", "w")
    fn2 = open("BinaryAlloy_selected.txt","w")
    for i in range(len(results_a)):
        fn.write(str(results_a[i].material_id) + '   ')
        fn.write(results_a[i].formula_pretty + '   ')
        fn.write(str(results_a[i].is_magnetic) + '    ')
        fn.write(str(results_a[i].theoretical) + '    ')
        fn.write(str(results_a[i].is_stable) + '    ')
        fn.write(str(results_a[i].nsites) + '\n')
        if results_a[i].is_magnetic == False and results_a[i].nsites <= 4:
            fn2.write(str(results_a[i].material_id) + '   ')
            fn2.write(results_a[i].formula_pretty + '   ')
            fn2.write(str(results_a[i].is_magnetic) + '    ')
            fn2.write(str(results_a[i].theoretical) + '    ')
            fn2.write(str(results_a[i].is_stable) + '    ')
            fn2.write(str(results_a[i].nsites) + '\n')
        fname = results_a[i].material_id + "_" + results_a[i].formula_pretty + ".cif"
        CifWriter(results_a[i].structure).write_file(fname)
    fn.close()
    fn2.close()
