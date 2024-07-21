#! /usr/bin/env python3
import os, sys, glob, re
import csv
from optparse import OptionParser

def ParseInput(ArgsIn):
    UseMsg = '''
    Usage: python [script] [options] [target_dir]
    '''
    parser = OptionParser(usage=UseMsg)
    #parser.add_option('--enthalpy',dest='enthalpy',action='store',type='string',default=None,help='keyword for parsing the enthalpy from the output file')
    #parser.add_option('--entropy',dest='entropy',action='store',type='string',default=None,help='keyword for parsing the entropy from the output file')

    options, args = parser.parse_args(ArgsIn)
    if len(args) < 2:
        parser.print_help()
        sys.exit(0)
    return options, args

def parse_enthalpy_entropy(target_dir, options):
    #if options.enthalpy:
    #    start_keyword = "Translational Enthalpy"
    #    end_keyword = "Vibrational Enthalpy"
    #    result_file = "QChem_enthalpy_results.csv"
    #if options.entropy:
    #    start_keyword = "Translational Entropy"
    #    end_keyword = "Vibrational Entropy"
    #    result_file = "QChem_entropy_results.csv"
    #else:
    start_keyword = "Zero point vibrational energy"
    end_keyword = "Vibrational Entropy"
    result_file = "QChem_results.csv"
    
    output_files = glob.glob(os.path.join(target_dir, '*.out'))
    all_results = []

    for outfile in output_files:
        results = {'filename': outfile}
        start_parsing = False

        with open(outfile, 'r') as fr:
            for line in fr:
                if start_keyword in line:
                    start_parsing = True
                if start_parsing:
                    l_sp = line.split()
                    if len(l_sp) < 2:
                        continue
                    if 'Zero' in line:
                        results['zpve'] = float(l_sp[-2])
                    if 'Translational' in line and 'Enthalpy' in line:
                        results['trans_enthalpy'] = float(l_sp[-2])
                    if 'Rotational' in line and 'Enthalpy' in line:
                        results['rot_enthalpy'] = float(l_sp[-2])
                    if 'Vibrational' in line and 'Enthalpy' in line:
                        results['vib_enthalpy'] = float(l_sp[-2])

                    if 'Translational' in line and 'Entropy' in line:
                        results['trans_entropy'] = float(l_sp[-2])
                    if 'Rotational' in line and 'Entropy' in line:
                        results['rot_entropy'] = float(l_sp[-2])
                    if 'Vibrational' in line and 'Entropy' in line:
                        results['vib_entropy'] = float(l_sp[-2])
                    if end_keyword in line:
                        break
        
        all_results.append(results)
    
    if all_results:
        keys = all_results[0].keys()
        output_path = os.path.join(target_dir, result_file)
        with open(output_path, 'w', newline='') as output_file:
            dict_writer = csv.DictWriter(output_file, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(all_results)
    print(f"Q-Chem Results saved to '{os.path.join(target_dir, result_file)}'")

if __name__ == "__main__":
    options, args = ParseInput(sys.argv)
    target_dir = args[1]
    parse_enthalpy_entropy(target_dir,options)
