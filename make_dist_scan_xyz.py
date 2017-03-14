#A one-shot script to create compressed and stretched geometries

import os, glob, re, sys
import subprocess as sp
from optparse import OptionParser
import numpy as np

def ParseInput(ArgsIn):
   UseMsg = "python make_dist_scan_xyz.py [options] [xyz_seed] [frgm_file]"
   parser = OptionParser(usage=UseMsg)
   parser.add_option('--name_root',dest='name_root',action='store',type='str',default=None,help='The name root for the input files generated. Default name_root.frgm->name_root_coordvalue.in')
   parser.add_option('--dest_path',dest='dest_path',action='store',type='str',default='new_dist_scan_xyz',help='The relative path to the directory where the generated files should be dumped.')
   parser.add_option('--type',dest='type',action='store',type='str',default='absolute',help='Type of scan. Options are "absolute" and "percent". An absolut scan goes from --start_abs to --end_abs with the specified increment. A percent scan goes from --start_pct*dist to --end_pct*dist with the specified increment.')
   parser.add_option('--increment',dest='increment',action='store',type='float',default=0.05,help='Set the increment of the scan. Default 0.05 Angs or 5%')
   parser.add_option('--start_absolute',dest='start_absolute',action='store',type='float',default=1.0,help='Set the smallest value for coord. Relevant for the scan option "absolute". Default 1.0')
   parser.add_option('--end_absolute',dest='end_absolute',action='store',type='float',default=6.0,help='Set the largest value for coord. Relevant for the scan option "absolute". Default 6.0')
   parser.add_option('--start_percent',dest='start_percent',action='store',type='float',default=0.7,help='Set the smallest value for coord to this percent of the initial value. Relevant for the scan option "percent". Default 0.7')
   parser.add_option('--end_percent',dest='end_percent',action='store',type='float',default=3.0,help='Set the largest value for coord to this percent of the initial value. Relevant for the scan option "percent". Default 3.0')
   parser.add_option('--scan_coord',dest='scan_coord',action='callback',type='string', callback=string_sp_callback, default=None, help='specify the intermolecular distance')

   options, args = parser.parse_args(ArgsIn)
   if options.scan_coord == None or len(options.scan_coord)!=2:
      print "specify the distance we are going to scan (given by a pair of atom indices)"
      parser.print_help()
      sys.exit(0)

   if len(args) < 2:
      parser.print_help()
      sys.exit(0)

   return options, args
        
def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

class XYZ:
   def __init__(self, xyz_file):
      #print xyz_file
      self.Name = re.search("/([^/]+).xyz", xyz_file).group(1)
      f = open(xyz_file, 'r')
      self.AtomList = []
      self.CoordList = []
      for line in f.readlines():
         l = re.search('^\s*(\d+)\s*$', line)
         if l!=None:
            self.NAtom = int(l.group(1))
         l = re.search('(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
         if l!=None:
            self.AtomList.append(l.group(1))
            coord = float(l.group(2)), float(l.group(3)), float(l.group(4))
            self.CoordList.append(coord)

      if(len(self.AtomList)!=self.NAtom):
         print "Error in number of atoms"
         sys.exit(1)

def ParseFRGM(FRGMFile, debug = False):
   FRGM = {}
   counter = 0
   f = open(FRGMFile,'r')
   line = f.readline()
   while line!='':
      l=re.search("(\#)",line) 
      if (l==None):
         counter = counter + 1
         if counter == 1:
            l=re.search("(\S+)",line) 
            if (not l==None):
               FRGM["total_charge"] = int(l.group(1)) 
               if debug:
                  print "supersystem charge is "+ str(FRGM["total_charge"])
            else:
               print "error parsing frgm file"
         if counter == 2:
            l=re.search("(\S+)",line) 
            if (not l==None):
               FRGM["total_mult"] = int(l.group(1)) 
               if debug:
                  print "supersystem multiplicity is "+ str(FRGM["total_mult"])
            else:
               print "error parsing frgm file"
         if counter == 3:
            l=re.search("(\S+)",line) 
            if (not l==None):
               FRGM["charge_frgm"] = [int(i) for i in line.split()]
               if debug:
                  print "printing fragment charges"
                  print FRGM["charge_frgm"]
            else:
               print "error parsing frgm file"
         if counter == 4:
            l=re.search("(\S+)",line) 
            if (not l==None):
               FRGM["mult_frgm"] = [int(i) for i in line.split()]
               if debug:
                  print "printing fragment multiplicites"
                  print FRGM["mult_frgm"]
            else:
               print "error parsing frgm file"
         if counter == 5:
            l=re.search("(\S+)",line) 
            if (not l==None):
               FRGM["atoms_frgm"] = [int(i) for i in line.split()]
               if debug:
                  print "printing number of atoms in each fragment"
                  print FRGM["atoms_frgm"]
            else:
               print "error parsing frgm file"
            #make a list for number of atoms prior to each fragment (offset)
            FRGM["NFrgm"] = len(FRGM["atoms_frgm"])
            FRGM["atoms_offset"] = np.zeros(FRGM["NFrgm"],dtype=int)
            offset = 0
            for index in range(0,FRGM["NFrgm"]):
               FRGM["atoms_offset"][index] = offset
               offset += FRGM["atoms_frgm"][index]
      line = f.readline()
   if counter != 5:
      print "did not parse the expected number of lines in "+FRGMFile 
   return FRGM

def make_atom_to_frgm_list(XYZ, FRGM):
   atom_to_frgm_list = np.zeros(XYZ.NAtom, dtype=int)
   for ifrgm in range(0, FRGM["NFrgm"]):
      offset = FRGM["atoms_offset"][ifrgm]
      for iatom in range(0, FRGM["atoms_frgm"][ifrgm]):
         atom_to_frgm_list[offset+iatom] = ifrgm

   return atom_to_frgm_list


#set the distance between the two specified atoms to the given value
#shifted_frgm is atom1's fragment
def fragment_translate_abs(XYZ, FRGM, name_root, atom1_idx, atom2_idx, shifted_frgm, dist_abs, dest_path): 
   #compute the shift vector
   Coords = np.copy(XYZ.CoordList)
   orig_vector = Coords[atom1_idx-1] - Coords[atom2_idx-1]
   trans_vector = orig_vector/np.sqrt(np.dot(orig_vector, orig_vector))*dist_abs - orig_vector
   
   #figure out the atoms to shift
   offset = FRGM["atoms_offset"][shifted_frgm]
   for iatom in range(0, FRGM["atoms_frgm"][shifted_frgm]):
      Coords[offset+iatom] += trans_vector

   outfile = dest_path + '/'+name_root+'_dist_'+format(dist_abs, '.2f')+'.xyz'
   write_xyz_file(outfile, Coords, XYZ.AtomList)


#set the distance between the two specified atoms to a certain ratio of the eql distance
def fragment_translate_rel(XYZ, FRGM, name_root, atom1_idx, atom2_idx, shifted_frgm, ratio, dest_path): 
   #compute the shift vector
   Coords = np.copy(XYZ.CoordList)
   orig_vector = Coords[atom1_idx-1] - Coords[atom2_idx-1]
   trans_vector = orig_vector*(ratio - 1.0)
   
   #figure out the atoms to shift
   offset = FRGM["atoms_offset"][shifted_frgm]
   for iatom in range(0, FRGM["atoms_frgm"][shifted_frgm]):
      Coords[offset+iatom] += trans_vector

   outfile = dest_path + '/' + name_root+'_rel_dist_'+format(dist_abs, '.2f')+'.xyz'
   write_xyz_file(outfile, Coords, XYZ.AtomList)


def write_xyz_file(outfile, Coords, AtomList):
   fw = open(outfile, 'w')
   fw.write("%d\n" %len(AtomList))
   fw.write("\n")
   for iAtom in range(0, len(AtomList)):
      x, y, z = Coords[iAtom][0], Coords[iAtom][1], Coords[iAtom][2]
      fw.write("%-3s %15.10f %15.10f %15.10f\n" %(AtomList[iAtom], x, y, z))
   fw.close()


#the script
curdir = os.getcwd()
options, args = ParseInput(sys.argv)

xyz_seed = args[1]
parse_XYZ = XYZ(xyz_seed)
#print XYZ.AtomList
frgm_file = args[2]
FRGM = ParseFRGM(frgm_file)

name_root = options.name_root
if name_root == None:
   l = re.search('(\S+).frgm', frgm_file)
   if l!=None:
      name_root = l.group(1)
   else:
      print "missing a name root!"
      sys.exit(1) 

dest_path = options.dest_path
if dest_path[-1] == '/':
   dest_path = dest_path[:-1]
if not os.path.exists(dest_path):
   sp.call(['mkdir', dest_path])

atom_to_frgm_list = make_atom_to_frgm_list(parse_XYZ, FRGM)
atom1_idx = int(options.scan_coord[0])
atom2_idx = int(options.scan_coord[1])
shifted_frgm = atom_to_frgm_list[atom1_idx-1] 

#absolute
if options.type == 'absolute':
   start_value = float(options.start_absolute)
   end_value = float(options.end_absolute)
   increment = float(options.increment)
   dist_abs = start_value
   while dist_abs <= end_value + 1.0E-4:
      fragment_translate_abs(parse_XYZ, FRGM, name_root, atom1_idx, atom2_idx, shifted_frgm, dist_abs, dest_path)
      dist_abs += increment

#relative
elif options.type == 'percent':
   start_percent = float(options.start_percent)
   end_percent = float(options.end_percent)
   increment = float(options.increment)
   ratio = start_percent
   while ratio <= end_percent + 1.0E-4:
      fragment_translate_rel(parse_XYZ, FRGM, name_root, atom1_idx, atom2_idx, shifted_frgm, ratio, dest_path)
      ratio += increment

else:
   print "Unsupported type: " +options.type
   sys.exit(0)

