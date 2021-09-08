import os, sys, re, glob, math
import numpy as np

def parse_pop_section(outfile):
   tmpfile = "pop.txt"
   command = "sed -n " + "'/Ground-State Mulliken Net/,/Sum of atomic/p' " + outfile + " > " + tmpfile
   os.system(command)
   fw = open(tmpfile)
   AtomList = []
   PopData = []
   for line in fw.readlines():
      l_sp = line.split()
      if len(l_sp) == 0:
         continue
      elif l_sp[0] == "Atom":
         AtomList = []
         PopData = []
      elif len(l_sp) == 3:
         AtomList.append(l_sp[1])
         PopData.append((float(l_sp[2]), 0.0))
      elif len(l_sp) == 4:
         AtomList.append(l_sp[1])
         PopData.append( (float(l_sp[2]), float(l_sp[3])) )
   fw.close()
   os.system("rm " + tmpfile)
   PopData = np.array(PopData)
   return AtomList, PopData 

def get_frgm_net_charge(PopData, frgm_begin, frgm_end):
   return sum(PopData[frgm_begin-1 : frgm_end , 0])

def get_frgm_net_spin(PopData, frgm_begin, frgm_end):
   return sum(PopData[frgm_begin-1 : frgm_end , 1])
