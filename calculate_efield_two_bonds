#! /usr/bin/env python

import os, glob, sys, re, math
import numpy as np
import xyzgeom
import pandas as pd
from optparse import OptionParser

def ParseInput(ArgsIn):
    UseMsg = '''
    Usage: python [script] [options] [target_dir] [xyz_dir]
    '''
    parser = OptionParser(usage=UseMsg)
    parser.add_option('--cutoff_values', dest='cutoff_values', action='callback', type='string', default=None, callback=string_sp_callback, help='specify multiple cutoff values (default: 0A and 7A)')
    parser.add_option('--skip_efield', dest='skip_efield', action='store_true', default=False, help='set true when .efield files have been generated')
    parser.add_option('--no_subtract', dest='no_subtract', action='store_true', default=False, help='skip subtraction of solute field')
    parser.add_option('--second_only', dest='second_only', action='store_true', default=False, help='only parse the efield produced by the second job')
    parser.add_option('--probe_atoms', dest='probe_atoms', action='callback', type='string', default=None, callback=string_sp_callback, help='the indices of three atoms. Two bonds are 1->2 and 1->3')

    options, args = parser.parse_args(ArgsIn)
    if len(args) < 3:
        parser.print_help()
        sys.exit(0)
    if options.probe_atoms != None and len(options.probe_atoms) != 3:
        print("Has to specify the index of three probe atoms (C, D, O)")
        parser.print_help()
        sys.exit(0)

    return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def parse_efield(efld_file, infile, second_only = False):
    fr = open(infile, 'r')
    start_parsing = False
    begin = -1
    end = -1
    efield = np.array([])
    efield_temp = None
    N = -1
    job_counter = 0
    for line in fr.readlines():
        l_sp = line.split()
        if len(l_sp) == 1 and l_sp[0] == 'EField':
            job_counter += 1
            if not (second_only and job_counter == 1):
                start_parsing = True
                counter = 0
            continue
        if len(l_sp) > 2 and l_sp[1] == 'DONE':
            if not (second_only and job_counter == 1):
                break
        if start_parsing:
            if second_only and job_counter == 1:
                continue
            #print line
            if counter == 0: #atom index
                begin = int(l_sp[0])
                end = int(l_sp[-1])
                N = end - begin + 1
                efield_temp = np.zeros((N, 3))
                counter += 1
            else:
                for idx in range(N):
                    efield_temp[idx][counter-1] = float(l_sp[idx+1])
                if counter == 3:
                    if len(efield) == 0:
                        efield = efield_temp
                    else:
                        efield = np.vstack((efield, efield_temp))
                    counter = 0
                else:
                    counter += 1
    fr.close()
    #dump field to file
    np.savetxt(efld_file, efield, fmt="%.7f")

def generate_all_efld_files(second_only):
   outfile_list = glob.glob('*.out')
   for outfile in outfile_list:
       tmpfile = "tmp.txt"
       command = "sed -n " + "'/ EField /,/DONE ESP/p' " + outfile + " > " + tmpfile
       os.system(command)
       efld_file = outfile + ".efld"
       parse_efield(efld_file, tmpfile, second_only)
   os.system('rm '+tmpfile)
   

def get_bond_direction(xyzfile, atom_idx1, atom_idx2):
    atomlist, coords = xyzgeom.parse_xyz_file(xyzfile)
    direction = coords[atom_idx2 - 1] - coords[atom_idx1 - 1]
    direction /= np.linalg.norm(direction)
    return direction

def get_bond_direction_all_frames(xyz_path, atom_idx_list):
    bond_direction = {}
    xyzfile_list = glob.glob(xyz_path+'/'+'r*.xyz')
    for xyzfile in xyzfile_list:
        frame = int(re.search("frame_([^_]+).xyz", xyzfile).group(1))
        bond_direction[frame] = {}
        bond_direction[frame]["bond1"] = get_bond_direction(xyzfile, atom_idx_list[0], atom_idx_list[1])
        bond_direction[frame]["bond2"] = get_bond_direction(xyzfile, atom_idx_list[0], atom_idx_list[2])
    return bond_direction

def get_efield_single_cutoff(data_efield, target_dir, cutoff, atom_idx_list):
    efield_file_list = glob.glob(target_dir+"/r{:.2f}_frame*.efld".format(cutoff))
    #print "Cutoff = %.1f, processing %d E-field files" %(cutoff, len(efield_file_list))
    data_efield[cutoff] = {}
    for efield_file in efield_file_list:
        #print (efield_file)
        data = np.loadtxt(efield_file)
        frame = int(re.search("frame_([^_]+)_", efield_file).group(1))
        data_efield[cutoff][frame] = {}
        for idx in atom_idx_list:
            data_efield[cutoff][frame][idx] = data[idx-1]
    return data_efield

def calculate_efield_on_atom(data_efield, bond_direction, atom_idx_list):
    au_to_MVcm = 5.142E+3
    data_efield_atom = {}
    for cutoff in data_efield:
        data_efield_atom[cutoff] = {}
        sum_bond1_1 = sum2_bond1_1 = 0.0
        sum_bond1_2 = sum2_bond1_2 = 0.0
        sum_bond2_1 = sum2_bond2_1 = 0.0
        sum_bond2_2 = sum2_bond2_2 = 0.0
        for frame in data_efield[cutoff]:
            dir_bond1 = bond_direction[frame]["bond1"]
            dir_bond2 = bond_direction[frame]["bond2"]
            field_atm1 = data_efield[cutoff][frame][atom_idx_list[0]] * au_to_MVcm
            field_atm2 = data_efield[cutoff][frame][atom_idx_list[1]] * au_to_MVcm
            field_atm3 = data_efield[cutoff][frame][atom_idx_list[2]] * au_to_MVcm
            field_bond1_1 = np.dot(field_atm1, dir_bond1)
            field_bond1_2 = np.dot(field_atm2, dir_bond1)
            field_bond2_1 = np.dot(field_atm1, dir_bond2)
            field_bond2_2 = np.dot(field_atm3, dir_bond2)
            sum_bond1_1 += field_bond1_1
            sum_bond1_2 += field_bond1_2
            sum_bond2_1 += field_bond2_1
            sum_bond2_2 += field_bond2_2
            sum2_bond1_1 += field_bond1_1 ** 2
            sum2_bond1_2 += field_bond1_2 ** 2
            sum2_bond2_1 += field_bond2_1 ** 2
            sum2_bond2_2 += field_bond2_2 ** 2
        n_frames = len(data_efield[cutoff])
        data_efield_atom[cutoff]["avg_bond1_1"] = sum_bond1_1 / float(n_frames)
        data_efield_atom[cutoff]["avg_bond1_2"] = sum_bond1_2 / float(n_frames)
        data_efield_atom[cutoff]["avg_bond2_1"] = sum_bond2_1 / float(n_frames)
        data_efield_atom[cutoff]["avg_bond2_2"] = sum_bond2_2 / float(n_frames)
        data_efield_atom[cutoff]["std_bond1_1"] = math.sqrt(sum2_bond1_1 / float(n_frames) - data_efield_atom[cutoff]["avg_bond1_1"]**2)
        data_efield_atom[cutoff]["std_bond1_2"] = math.sqrt(sum2_bond1_2 / float(n_frames) - data_efield_atom[cutoff]["avg_bond1_2"]**2)
        data_efield_atom[cutoff]["std_bond2_1"] = math.sqrt(sum2_bond2_1 / float(n_frames) - data_efield_atom[cutoff]["avg_bond2_1"]**2)
        data_efield_atom[cutoff]["std_bond2_2"] = math.sqrt(sum2_bond2_2 / float(n_frames) - data_efield_atom[cutoff]["avg_bond2_2"]**2)
        data_efield_atom[cutoff]["se_bond1_1"] = data_efield_atom[cutoff]["std_bond1_1"] / math.sqrt(float(n_frames))
        data_efield_atom[cutoff]["se_bond1_2"] = data_efield_atom[cutoff]["std_bond1_2"] / math.sqrt(float(n_frames))
        data_efield_atom[cutoff]["se_bond2_1"] = data_efield_atom[cutoff]["std_bond2_1"] / math.sqrt(float(n_frames))
        data_efield_atom[cutoff]["se_bond2_2"] = data_efield_atom[cutoff]["std_bond2_2"] / math.sqrt(float(n_frames))
    return data_efield_atom

def calculate_rel_efield(data_efield, bond_direction, target_dir, atom_idx_list, no_subtract = False):
    au_to_MVcm = 5.142E+3
    idx_1, idx_2, idx_3 = atom_idx_list[0], atom_idx_list[1], atom_idx_list[2]
    data_efield_stat = {}
    for cutoff in data_efield:
        efield_data_file = target_dir +'/'+"efield_per_frame_"+str(cutoff)+"A.csv"
        fw = open(efield_data_file, 'w')
        fw.write("frame,E_bond1,E_bond2\n")
        data_efield_stat[cutoff] = {}
        sum_E_bond1 = sum2_E_bond1 = 0.0
        sum_E_bond2 = sum2_E_bond2 = 0.0
        for frame in sorted(data_efield[cutoff]):
            dir_bond1 = bond_direction[frame]["bond1"]
            dir_bond2 = bond_direction[frame]["bond2"]
            if not no_subtract:
                field_atm1 = (data_efield[cutoff][frame][idx_1] - data_efield[0][frame][idx_1]) * au_to_MVcm
                field_atm2 = (data_efield[cutoff][frame][idx_2] - data_efield[0][frame][idx_2]) * au_to_MVcm
                field_atm3 = (data_efield[cutoff][frame][idx_3] - data_efield[0][frame][idx_3]) * au_to_MVcm
            else:
                field_atm1 = data_efield[cutoff][frame][idx_1] * au_to_MVcm
                field_atm2 = data_efield[cutoff][frame][idx_2] * au_to_MVcm
                field_atm3 = data_efield[cutoff][frame][idx_3] * au_to_MVcm
            field_bond1_1 = np.dot(field_atm1, dir_bond1)
            field_bond1_2 = np.dot(field_atm2, dir_bond1)
            field_bond2_1 = np.dot(field_atm1, dir_bond2)
            field_bond2_2 = np.dot(field_atm3, dir_bond2)
            E_bond1 = 0.5 * (field_bond1_1 + field_bond1_2)
            E_bond2 = 0.5 * (field_bond2_1 + field_bond2_2)
            fw.write("%d,%.3f,%.3f\n" %(frame, E_bond1, E_bond2))
            sum_E_bond1 += E_bond1
            sum2_E_bond1 += E_bond1 ** 2
            sum_E_bond2 += E_bond2
            sum2_E_bond2 += E_bond2 ** 2
        n_frames = len(data_efield[cutoff])
        print("Cutoff = %.1f A, Number of frames: %d" %(cutoff, n_frames))
        fw.close()
        data_efield_stat[cutoff]["avg_E_bond1"] = sum_E_bond1 / float(n_frames)
        data_efield_stat[cutoff]["std_E_bond1"] = math.sqrt(sum2_E_bond1 / float(n_frames) - data_efield_stat[cutoff]["avg_E_bond1"]**2)
        data_efield_stat[cutoff]["se_E_bond1"] = data_efield_stat[cutoff]["std_E_bond1"] / math.sqrt(float(n_frames))
        data_efield_stat[cutoff]["avg_E_bond2"] = sum_E_bond2 / float(n_frames)
        data_efield_stat[cutoff]["std_E_bond2"] = math.sqrt(sum2_E_bond2 / float(n_frames) - data_efield_stat[cutoff]["avg_E_bond2"]**2)
        data_efield_stat[cutoff]["se_E_bond2"] = data_efield_stat[cutoff]["std_E_bond2"] / math.sqrt(float(n_frames))
    return data_efield_stat


options, args = ParseInput(sys.argv) 
curdir = os.getcwd()
print("Generating all E-field files from Q-Chem output...")
target_dir = args[1]
if not options.skip_efield:
    os.chdir(target_dir)
    generate_all_efld_files(options.second_only)
    os.chdir(curdir)

#probe atoms
atom_idx_list = []
for iatm in options.probe_atoms:
   atom_idx_list.append(int(iatm))

xyz_path = args[2]
if xyz_path[:-1] == '/':
    xyz_path = xyz_path[:-1]
bond_direction = get_bond_direction_all_frames(xyz_path, atom_idx_list)

print("Getting the E-field vectors on relevant atoms")
data_efield = {}

cutoff_values = []
if options.cutoff_values == None:
    if options.no_subtract:
        cutoff_values = [7]
    else:
        cutoff_values = [0, 7]  #default
else:
    [cutoff_values.append(int(val)) for val in options.cutoff_values]
for cutoff in cutoff_values:
    get_efield_single_cutoff(data_efield, target_dir, cutoff, atom_idx_list)

print("Calculating E-fields along two bond directions")
data_relE_stat = calculate_rel_efield(data_efield, bond_direction, target_dir, atom_idx_list, options.no_subtract)
data_efield_on_atom = calculate_efield_on_atom(data_efield, bond_direction, atom_idx_list)
#print E-field on chemical bonds
df_relE = pd.DataFrame.from_dict(data_relE_stat, orient="index")
df_relE.index.name = 'cutoff'
for cutoff in sorted(cutoff_values):
   if cutoff == 0:
      continue
   print("E-field along bond 1 (atm1->atm2, cutoff = %dA):" %cutoff)
   print("Mean: %.3f MV/cm, Stdev: %.3f MV/cm, StdErr: %.3f MV/cm" %(df_relE["avg_E_bond1"][cutoff], df_relE["std_E_bond1"][cutoff], df_relE["se_E_bond1"][cutoff]))

   print("E-field along bond 2 (atm1->atm3, cutoff = %dA):" %cutoff)
   print("Mean: %.3f MV/cm, Stdev: %.3f MV/cm, StdErr: %.3f MV/cm" %(df_relE["avg_E_bond2"][cutoff], df_relE["std_E_bond2"][cutoff], df_relE["se_E_bond2"][cutoff]))
df_relE = df_relE.reindex(sorted(df_relE.columns), axis=1)
df_relE.to_csv(target_dir + '/efield_bond.csv')

#save E-field on atoms to csvfile
df_efield_atom = pd.DataFrame.from_dict(data_efield_on_atom, orient="index")
df_efield_atom.index.name = "cutoff"
df_efield_atom = df_efield_atom.reindex(sorted(df_efield_atom), axis=1)
df_efield_atom.to_csv(target_dir + '/efield_atom.csv')
