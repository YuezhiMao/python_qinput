import os, sys, glob, re
import subprocess as sp
import numpy as np

def get_analytic_forces(output):
    tmpfile = 'tmp.txt'
    command = "sed -n '/Gradient of/,/Max gradient/p' " + output + " > " + tmpfile
    os.system(command)
    coord_x = []
    coord_y = []
    coord_z = []
    fr = open(tmpfile, 'r')
    for line in fr.readlines():
        l_sp = line.split()
        if l_sp[0] == "Gradient" or l_sp[0] == "Max":
            continue
        elif l_sp[0] == "1" and l_sp[1] != "2":
            for i in range(1, len(l_sp)):
                coord_x.append(float(l_sp[i]))
        elif l_sp[0] == "2":
            for i in range(1, len(l_sp)):
                coord_y.append(float(l_sp[i]))
        elif l_sp[0] == "3":
            for i in range(1, len(l_sp)):
                coord_z.append(float(l_sp[i]))
    print "length of x, y, z vectors: %d %d %d" %(len(coord_x), len(coord_y), len(coord_z))
    fr.close()
    sp.call(['rm', tmpfile])
    #dump to csvfile
    force_filename = output[:-4]+'.force.csv'
    fw = open(force_filename, 'w')
    fw.write('X,Y,Z\n')
    for i in range(len(coord_x)):
        fw.write("%.7f,%.7f,%.7f\n" %(coord_x[i], coord_y[i], coord_z[i]))
    fw.close()
    #return numpy array
    all_forces = np.array(zip(coord_x, coord_y, coord_z))
    return all_forces

def get_numerical_forces(output):
    tmpfile = 'tmp.txt'
    command = "sed -n '/Order 1/,/Archival/p' " + output + " > " + tmpfile
    os.system(command)
    forces = []
    fr = open(tmpfile, 'r')
    for line in fr.readlines():
        l_sp = line.split()
        if len(l_sp) == 0: continue
        if re.match('(\d+)', l_sp[0]):
            #print line
            x, y, z = float(l_sp[1]), float(l_sp[2]), float(l_sp[3])
            forces.append((x, y, z))
    fr.close()
    sp.call(['rm', tmpfile])
    #dump to csvfile
    force_filename = output[:-4]+'.force.csv'
    fw = open(force_filename, 'w')
    fw.write('X,Y,Z\n')
    for i in range(len(forces)):
        fw.write("%.7f,%.7f,%.7f\n" %(forces[i][0], forces[i][1], forces[i][2]))
    fw.close()
    #return numpy array
    all_forces = np.array(forces)
    return all_forces
