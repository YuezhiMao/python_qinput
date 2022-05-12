import os, glob, re, sys
import subprocess as sp
import numpy as np

class XYZ:  #class for XYZ coordinates
   def __init__(self, xyz_file):
      #print (xyz_file)
      self.Name = re.search("([^/]+).xyz$", xyz_file).group(1)
      f = open(xyz_file, 'r')
      self.AtomList = []
      self.CoordList = []
      for line in f.readlines():
         l = re.search('^\s*(\d+)\s*$', line)
         if l!=None:
            self.NAtom = int(l.group(1))
         l = re.search('^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
         if l!=None:
            self.AtomList.append(l.group(1))
            coord = float(l.group(2)), float(l.group(3)), float(l.group(4))
            self.CoordList.append(coord)

      if(len(self.AtomList)!=self.NAtom):
         print("Error in number of atoms: "+xyzfile)
         sys.exit(1)

class FRGM: #class for fragment partition information
   def __init__(self, frgm_file):
      fr = open(frgm_file, 'r')
      line = fr.readline()
      counter = 0
      while line != '':
         counter = counter + 1
         ls = line.split()
         if counter == 1:
            if len(ls) != 1:
               print("Error in line %d of %s" %(counter, frgm_file))
               sys.exit(0)
            else:
               self.total_charge = int(ls[0])
         elif counter == 2:
            if len(ls) != 1:
               print("Error in line %d of %s" %(counter, frgm_file))
               sys.exit(0)
            else:
               self.total_mult = int(ls[0])
         elif counter == 3:
            self.n_frgm = len(ls)
            if self.n_frgm < 2:
               print("Error in line %d of %s: At least two fragments required" %(counter, frgm_file))
               sys.exit(0)
            self.charge_frgm = [int(i) for i in ls]
         elif counter == 4:
            if len(ls) != self.n_frgm:
               print("Error in line %d of %s" %(counter, frgm_file))
               sys.exit(0)
            self.mult_frgm = [int(i) for i in ls]
         elif counter == 5:
            if len(ls) != self.n_frgm:
               print("Error in line %d of %s" %(counter, frgm_file))
               sys.exit(0)
            self.natoms_frgm = [int(i) for i in ls] 
            
         else:
            print("skip the unexpected lines in " + frgm_file) 
         line = fr.readline()

      self.atoms_offset = np.zeros(self.n_frgm, dtype=int)
      offset = 0
      for index in range(0, self.n_frgm):
         self.atoms_offset[index] = offset
         offset += self.natoms_frgm[index] 

def WriteMolecule(fw, XYZ, charge, mult):
   fw.write('$molecule\n')
   #total charge and mult
   fw.write("%d %d\n" %(charge, mult))
   for iatom in range(0, XYZ.NAtom):
      atomic_symbol = XYZ.AtomList[iatom]
      x, y, z = XYZ.CoordList[iatom]
      fw.write("%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z))
   fw.write('$end\n')
   fw.write('\n')

def WriteMolecule_Frgm(fw, XYZ, FRGM):
   #print "Working on "+XYZ.Name
   fw.write('$molecule\n')
   #total charge and mult
   fw.write("%d %d\n" %(FRGM.total_charge, FRGM.total_mult))
   #loop over fragments
   for ifrgm in range(0, FRGM.n_frgm):
      fw.write('--\n')
      fw.write("%d %d\n" %(FRGM.charge_frgm[ifrgm], FRGM.mult_frgm[ifrgm]))
      natom_frgm = FRGM.natoms_frgm[ifrgm]
      index_start = FRGM.atoms_offset[ifrgm]
      for iatom in range(index_start, index_start+natom_frgm):
         atomic_symbol = XYZ.AtomList[iatom]
         x, y, z = XYZ.CoordList[iatom]
         fw.write("%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z))
   fw.write('$end\n')
   fw.write('\n')

def WriteMolecule_supersub(fw, XYZ, FRGM, monomer=0, ghost=False):
   fw.write('$molecule\n')
   if monomer == 0:   #supersystem
      if ghost:
         print("Can't add ghost atoms for the supersystem calculation")
         sys.exit(0)
      else:
         fw.write("%d %d\n" %(FRGM.total_charge, FRGM.total_mult))
         for index in range(0, XYZ.NAtom):
            atomic_symbol = XYZ.AtomList[index]
            x, y, z = XYZ.CoordList[index]
            fw.write("%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z)) 
   else:   #monomer job
      index_frgm = monomer - 1   #start from zero
      index_start = FRGM.atoms_offset[index_frgm]
      index_end = index_start+FRGM.natoms_frgm[index_frgm]
      fw.write("%d %d\n" %(FRGM.charge_frgm[index_frgm], FRGM.mult_frgm[index_frgm]))
      if not ghost:   #isolated monomer
         for index in range(index_start, index_end):
            atomic_symbol = XYZ.AtomList[index]
            x, y, z = XYZ.CoordList[index]
            fw.write("%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z)) 

      else:  #monomer with ghost
         for index in range(0, XYZ.NAtom):
            atomic_symbol = XYZ.AtomList[index]
            x, y, z = XYZ.CoordList[index]
            if index in range(index_start, index_end):
               fw.write("%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z)) 
            else:  #ghost atoms
               fw.write("@%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z)) 
   fw.write("$end\n")
   fw.write("\n")


def WriteMolecule_Read(fw):
   fw.write('$molecule\n')
   fw.write('read\n')
   fw.write('$end\n\n')

def detect_unrestricted_frgm(FRGM):
    for ifrgm in range(0, FRGM.n_frgm):
        if FRGM.mult_frgm[ifrgm] > 1:
            return True
            break
    return False

def ParseMolFile(MoleculeFile):
   MOLECULE = {}
   MOLECULE["name"] = re.search('([^\/]+).mol', MoleculeFile).group(1)
   MOLECULE["nlines"] = 0
   curline = 0
   MOLECULE["the_line"] = {}
   MOLECULE["unrestricted"] = False
   MOLECULE["fragmented"] = False
   f = open(MoleculeFile,'r')
   line = f.readline()
   mult_checked = False
   while line!='':
      MOLECULE["the_line"][str(curline)] = line
      curline += 1
      MOLECULE["nlines"] = curline
      l=re.search("^\s*(\S+)\s+(\S+)\s*$",line) #charge and mult
      if (l!=None and not mult_checked):
         if (int(l.group(2)) != 1):	#multiplicity larger than 1
            MOLECULE["unrestricted"] = True
         mult_checked=True
      l=re.search("(--)",line) #fragment indicator
      if (not l==None):
         MOLECULE["fragmented"] = True
      line = f.readline()
   f.close()
   return MOLECULE
				
def WriteMolSection(fw, MyMOLECULE):
   fw.write('$molecule\n')
   for line in range(MyMOLECULE["nlines"]):
      fw.write(MyMOLECULE["the_line"][str(line)])
   fw.write('\n')
   fw.write('$end\n')
   fw.write('\n')
   return
