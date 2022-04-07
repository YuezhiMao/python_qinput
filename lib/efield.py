import os, glob, sys, re, math
import numpy as np
import xyzgeom

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

def generate_efield_file(efld_file, outfile, second_only):
    tmpfile = outfile+".tmp"
    command = "sed -n " + "'/ EField /,/DONE ESP/p' " + outfile + " > " + tmpfile
    os.system(command)
    parse_efield(efld_file, tmpfile, second_only)
    os.system('rm '+tmpfile)
   
def get_bond_direction(xyzfile, atom_idx1, atom_idx2):
    atomlist, coords = xyzgeom.parse_xyz_file(xyzfile)
    direction = coords[atom_idx2 - 1] - coords[atom_idx1 - 1]
    direction /= np.linalg.norm(direction)
    return direction

def calculate_efield_on_bond(outfile, xyzfile, bond_idx1, bond_idx2, second_only=False):
   efld_file = outfile + ".efld"
   generate_efield_file(efld_file, outfile, second_only)
   au_to_MVcm = 5.142E+3
   dir_bond = get_bond_direction(xyzfile, bond_idx1, bond_idx2)
   efield_all = np.loadtxt(efld_file)
   E_atom1 = np.dot(efield_all[bond_idx1-1, :], dir_bond) * au_to_MVcm
   E_atom2 = np.dot(efield_all[bond_idx2-1, :], dir_bond) * au_to_MVcm
   E_bond = 0.5 * (E_atom1 + E_atom2)
   return E_bond

#This function is useful when multiple bonds are needed for one output (more efficient)
def calculate_efield_on_bond(efld_file, xyzfile, bond_idx1, bond_idx2):
   au_to_MVcm = 5.142E+3
   dir_bond = get_bond_direction(xyzfile, bond_idx1, bond_idx2)
   efield_all = np.loadtxt(efld_file)
   E_atom1 = np.dot(efield_all[bond_idx1-1, :], dir_bond) * au_to_MVcm
   E_atom2 = np.dot(efield_all[bond_idx2-1, :], dir_bond) * au_to_MVcm
   E_bond = 0.5 * (E_atom1 + E_atom2)
   return E_bond

def parse_esp(esp_file, tmpfile):
   fr = open(tmpfile, 'r')
   raw_data = []
   count = 0
   for line in fr.readlines():
      l_sp = line.split()
      if l_sp[0] == 'TotESP':
         count += 1
      if len(l_sp) == 2 and re.search('(\d+)', l_sp[0]) != None:
         raw_data.append(float(l_sp[1]))
   fr.close()
   if count == 2:
      esp = np.array(raw_data[len(raw_data)/2 : ])
   else:
      esp = np.array(raw_data)
   np.savetxt(esp_file, esp, fmt="%.7f")

def generate_esp_file(esp_file, outfile):
   tmpfile = "tmp.txt"
   command = "sed -n " + "'/TotESP/,/Ele EField/p' " + outfile + " > " + tmpfile
   os.system(command)
   parse_esp(esp_file, tmpfile)
   os.system('rm ' + tmpfile)

def get_bond_length(xyzfile, atom_idx1, atom_idx2):
    atomlist, coords = xyzgeom.parse_xyz_file(xyzfile)
    bond_length = xyzgeom.compute_distance(coords, atom_idx1, atom_idx2)
    return bond_length

def calculate_efield_from_esp(outfile, xyzfile, atom_idx1, atom_idx2):
   esp_file = outfile + ".esp"
   generate_esp_file(esp_file, outfile)
   au_to_MVcm = 5.142E+3
   esp_all = np.loadtxt(esp_file)
   V_1 = esp_all[atom_idx1-1]
   V_2 = esp_all[atom_idx2-1]
   r_bond = get_bond_length(xyzfile, atom_idx1, atom_idx2) * 1.88973 #convert from A to Bohr (au)
   E_bond = (V_1 - V_2) / r_bond * au_to_MVcm
   return E_bond

#This function is useful when multiple bonds are needed for one output (more efficient)
def calculate_efield_from_esp(esp_file, xyzfile, atom_idx1, atom_idx2):
   au_to_MVcm = 5.142E+3
   esp_all = np.loadtxt(esp_file)
   V_1 = esp_all[atom_idx1-1]
   V_2 = esp_all[atom_idx2-1]
   r_bond = get_bond_length(xyzfile, atom_idx1, atom_idx2) * 1.88973 #convert from A to Bohr (au)
   E_bond = (V_1 - V_2) / r_bond * au_to_MVcm
   return E_bond
