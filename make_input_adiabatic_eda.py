import os, glob, re, sys
import subprocess as sp
import numpy as np
from optparse import OptionParser

def ParseInput(ArgsIn):
   UseMsg = '''
   python make_input_adiabatic_eda.py [options] [xyz_path] [frgm_partition]
   ''' 
   parser = OptionParser(usage=UseMsg)
   parser.add_option('-a','--all',dest='all',action='store_true',default=False,help='run all the xyz files under the xyz_path')
   parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='create inputs for xyz files containing the keyword')
   parser.add_option('-t','--target',dest='target',action='callback',callback=string_sp_callback,type='string',default=None,help='create input for certain xyz files')
   parser.add_option('--onsite',dest='onsite',action='store_true',default=False,help='frgm file and xyz file under the same dir and sharing the same nameroot')
   parser.add_option('-b','--basis',dest='basis',action='callback',callback=string_sp_callback, type='string', default='def2-tzvppd',help='The target basis (default is def2-tzvppd)')
   parser.add_option('-m','--method',dest='method',action='callback',type='string', default='wB97M-V', callback=string_sp_callback,help='The methods (density functionals) to use')
   parser.add_option('-i','--input_path',dest='input_path',action='store',type='string', default='input/',help='The directory storing the generated inputs')
   parser.add_option('--pes',dest='pestype',action='store',type='string',default='integrated',help='Integrated or separated optimizations on FRZ, POL or fully relaxed PESs')
   parser.add_option('--sym_ignore', dest='sym_ignore', action='store_true', default=False, help='ignore point group symmetry during geometry optimization')
   parser.add_option('--do_freq', dest='do_freq', action='store_true', default=False, help='do frequency calculation at the minimum of each PES')
   parser.add_option('--do_ferf', dest='do_ferf', action='store_true', default=False, help='use FERF to construct the polarized PES')
   parser.add_option('--loose', dest='loose', action='store_true', default=False, help='Use looser convergence criterion for geom opt')
    
   options, args = parser.parse_args(ArgsIn)

   if not options.onsite:
      if len(args) < 3:
         parser.print_help()
         sys.exit(1)
   else:
      if len(args) < 2:
         parser.print_help()
         sys.exit(1)

   if not options.all and options.target==None and options.keyword==None:
      print "The target directory must be specified: one or some or all"
      parser.print_help()
      sys.exit(1)
   return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))
   
class XYZ:
   def __init__(self, xyz_file):
      #print xyz_file
      self.Name = re.search("/([^/]+).xyz$", xyz_file).group(1)
      nameroot_parser = re.search('(\S+)_([^_]+)geom', self.Name)
      if nameroot_parser!=None:
         self.NameRoot = nameroot_parser.group(1)
      else:
         self.NameRoot = self.Name
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

def ParseRems(RemFile):
   f = open(RemFile,'r')
   REMS = {}
   REMS["name"] = RemFile[4:]
   REMS["nrems"] = 0
   currem = 0
   REMS["the_rem"] = {}
   line = f.readline()
   while line!='':
      l=re.search("(\S+)\s+(\S+)",line) #rem and value
      if (not l==None):
         REMS["the_rem"][str(currem)] = {}
         REMS["the_rem"][str(currem)]["name"] = l.group(1).upper()
         REMS["the_rem"][str(currem)]["value"] = l.group(2).upper()
         currem += 1
         REMS["nrems"] = currem
      line = f.readline()
   f.close()
   return REMS

def ModRem(r_name,r_value,MyREMS):
	#see if the rem is already set
   REM_NAME = r_name.upper()
   REM_VALUE = r_value.upper()
   done = False
   for rem in range(MyREMS["nrems"]):
      if (MyREMS["the_rem"][str(rem)]["name"] == REM_NAME):
	#found in the rem already set up
         MyREMS["the_rem"][str(rem)]["value"] = REM_VALUE
         done = True
   if (not done):	#not find, add this new rem
      currem = MyREMS["nrems"]
      MyREMS["the_rem"][str(currem)] = {}
      MyREMS["the_rem"][str(currem)]["name"] = REM_NAME
      MyREMS["the_rem"][str(currem)]["value"] = REM_VALUE
      MyREMS["nrems"] = currem + 1
   return

def ModRem_Frgm(r_name, r_value, MyREMS):
   REM_NAME = r_name.upper()
   REM_VALUE = r_value.upper()
   done = False
   for rem_frgm in range(MyREMS["nrems_frgm"]):
      if (MyREMS["rem_frgm"][str(rem_frgm)]["name"] == REM_NAME):
	     #found in the rem already set up
         MyREMS["rem_frgm"][str(rem_frgm)]["value"] = REM_VALUE
         done = True
   if (not done):	#not find, add this new rem
      currem = MyREMS["nrems_frgm"]
      MyREMS["rem_frgm"][str(currem_frgm)] = {}
      MyREMS["rem_frgm"][str(currem_frgm)]["name"] = REM_NAME
      MyREMS["rem_frgm"][str(currem_frgm)]["value"] = REM_VALUE
      MyREMS["nrems_frgm"] = currem_frgm + 1

def AppendRem(fh,MyREMS):
   fh.write('$rem\n')
   for rem in range(MyREMS["nrems"]):
      fh.write(MyREMS["the_rem"][str(rem)]["name"]+'  '+MyREMS["the_rem"][str(rem)]["value"]+'\n')
   fh.write('$end\n')
   return

def WriteMolecule_Frgm(fw, XYZ, FRGM):
   fw.write('$molecule\n')
   #total charge and mult
   fw.write("%d %d\n" %(FRGM["total_charge"], FRGM["total_mult"]))
   #loop over fragments
   for ifrgm in range(0, FRGM["NFrgm"]):
      fw.write('--\n')
      fw.write("%d %d\n" %(FRGM["charge_frgm"][ifrgm], FRGM["mult_frgm"][ifrgm]))
      natom_frgm = FRGM["atoms_frgm"][ifrgm]
      index_start = FRGM["atoms_offset"][ifrgm]
      for iatom in range(index_start, index_start+natom_frgm):
         atomic_symbol = XYZ.AtomList[iatom]
         x, y, z = XYZ.CoordList[iatom]
         fw.write("%-4s %-10.5f %-10.5f %-10.5f\n" %(atomic_symbol, x, y, z))
   fw.write('$end\n')
   fw.write('\n')

def WriteMolecule_Read(fw):
   fw.write('$molecule\n')
   fw.write('read\n')
   fw.write('$end\n\n')

def basis_abbr(BasName):    #generating abbreviated name for basis sets, using for input file names
   fullname = BasName.lower()
   if fullname == 'cc-pvdz':
      return 'dz'
   elif fullname == 'cc-pvtz':
      return 'tz'
   elif fullname == 'cc-pvqz':
      return 'qz'
   elif fullname == 'aug-cc-pvdz':
      return 'adz'
   elif fullname == 'aug-cc-pvtz':
      return 'atz'
   elif fullname == 'aug-cc-pvqz':
      return 'aqz'
   elif fullname == 'def2-svpd' or fullname == 'svpd':
      return 'svpd'
   elif fullname == 'def2-tzvpp' or fullname == 'tzvpp':
      return 'tzvpp'
   elif fullname == 'def2-tzvppd' or fullname == 'tzvppd':
      return 'tzvppd'
   elif fullname == 'def2-qzvpp' or fullname == 'qzvpp':
      return 'qzvpp'
   elif fullname == 'def2-qzvppd' or fullname == 'qzvppd':
      return 'qzvppd'
   elif fullname == 'pc-2':
      return 'pc2'
   elif fullname == 'pc-3':
      return 'pc3'
   elif fullname == 'aug-pc-2':
      return 'apc2'
   elif fullname == 'aug-pc-3':
      return 'apc3'
   elif fullname == 'dp':
      return 'dp'
   elif fullname == 'tp':
      return 'tp' 
   elif fullname == 'lp':
      return 'lp'
   else:
      print "unrecognized basis choice"
      sys.exit(1)

def set_common_rems(curREM, method, basis, options): 
   ModRem('JOBTYPE', 'OPT', curREM)
   ModRem('GEN_SCFMAN','TRUE', curREM)
   ModRem('MEM_TOTAL', '8000', curREM)
   ModRem('MEM_STATIC','2000', curREM)
   ModRem('SCF_GUESS', 'FRAGMO', curREM)
   ModRem('IDERIV', '1', curREM)
   if options.loose:
      ModRem('SCF_CONVERGENCE', '8', curREM)
      ModRem('GEOM_OPT_TOL_GRADIENT', '200', curREM)
      ModRem('GEOM_OPT_TOL_DISPLACEMENT', '200', curREM)
      ModRem('GEOM_OPT_TOL_ENERGY', '20', curREM)
   else:
      ModRem('SCF_CONVERGENCE', '9', curREM)
      ModRem('GEOM_OPT_TOL_GRADIENT', '100', curREM)
      ModRem('GEOM_OPT_TOL_DISPLACEMENT', '100', curREM)
      ModRem('GEOM_OPT_TOL_ENERGY', '10', curREM)
   if (method != 'HF'):
      ModRem ('XC_GRID', '000075000302', curREM)		
      ModRem ('NL_GRID', '1', curREM)
   if options.sym_ignore:
      ModRem('SYM_IGNORE', 'TRUE', curREM)
   else:
      ModRem('SYM_IGNORE', 'FALSE', curREM)
   if basis.upper() == 'SVPD':
      ModRem('BASIS','DEF2-SVPD',curREM)
   elif basis.upper() == 'TZVPP':
      ModRem('BASIS','DEF2-TZVPP',curREM)
   elif basis.upper() == 'TZVPPD':
      ModRem('BASIS','DEF2-TZVPPD',curREM)
   elif basis.upper() == 'QZVPP':
      ModRem('BASIS','DEF2-QZVPP',curREM)
   elif basis.upper() == 'QZVPPD':
      ModRem('BASIS','DEF2-QZVPPD',curREM)
   elif basis.upper() == 'DP': #double-zeta Pople: 6-31+G(d)
      ModRem('BASIS','6-31+G(d)', curREM)
   elif basis.upper() == 'TP': #triple-zeta Pople: 6-311++G(2df,2pd)
      ModRem('BASIS','6-311++G(2df,2pd)',curREM)
   elif basis.upper() == 'LP':
      ModRem('BASIS','6-311++G(3df,3pd)', curREM)
   else: #copy the name of standard basis
      ModRem('BASIS', basis, curREM)

   ModRem('METHOD', method, curREM)	
   #deal with -D3 in the method
   if (method == 'PBE-D3'):
      ModRem('METHOD', 'PBE', curREM)
      ModRem('DFT_D','EMPIRICAL_GRIMME3', curREM)
   if (method == 'REVPBE-D3'):
      ModRem('METHOD', 'REVPBE', curREM)
      ModRem('DFT_D','EMPIRICAL_GRIMME3', curREM)
   if (method == 'B3LYP-D3'):
      ModRem('METHOD', 'B3LYP', curREM)
      ModRem('DFT_D','EMPIRICAL_GRIMME3', curREM)
   if (method == 'BLYP-D3'):
      ModRem('METHOD', 'BLYP', curREM)
      ModRem('DFT_D','EMPIRICAL_GRIMME3', curREM)
   if (method == 'TPSS-D3'):
      ModRem('METHOD', 'TPSS', curREM)
      ModRem('DFT_D','EMPIRICAL_GRIMME3', curREM)

   #default rem_frgm
   #ModRem_Frgm('MEM_TOTAL', '8000', curREM)
   #ModRem_Frgm('MEM_STATIC', '2000', curREM)
   #ModRem_Frgm('SCF_GUESS', 'SAD', curREM)

def set_ferf_rems(curREM):
   ModRem('IDERIV', '0', curREM)
   ModRem('POL_GEOM', 'FALSE', curREM)
   ModRem('SCFMI_MODE', '1', curREM)
   ModRem('CHILD_MP', 'TRUE', curREM)
   ModRem('CHILD_MP_ORDERS', '122', curREM) #nDQ
   ModRem('GUESS_ALWAYS', 'TRUE', curREM)

def XYZ_to_Input(fw, XYZ, FRGM, curREM, read_mol=False, do_freq=False, do_ferf=False):
   #print "generate %s" %inputfile
   #check if unrestricted (any open-shell fragment?)
   unrestricted=False
   for ifrgm in range(0, FRGM["NFrgm"]):
      if FRGM["mult_frgm"][ifrgm] > 1:
         unrestricted=True
         break
   if unrestricted:
      ModRem('UNRESTRICTED','TRUE',curREM)

   if do_ferf:
      set_ferf_rems(curREM)

   if not read_mol:
      WriteMolecule_Frgm(fw, XYZ, FRGM)
   else:
      WriteMolecule_Read(fw) 
   AppendRem(fw, curREM)

   #if frequency is desired, append a second job
   if do_freq:
      fw.write('\n@@@\n\n')
      WriteMolecule_Read(fw)
      #Use the current rem but modify the jobtype to freq
      ModRem('JOBTYPE', 'FREQ', curREM)
      AppendRem(fw, curREM)

#Generate an integrated adiabatic EDA job
def XYZ_to_Input_Integrated(fw, XYZ, FRGM, method, basis, options):
   #Fully relaxed PES first 
   curREM = ParseRems('rems/rem_scfgeom')
   set_common_rems(curREM, method, basis, options)
   XYZ_to_Input(fw, XYZ, FRGM, curREM, read_mol=False, do_freq=options.do_freq)
   fw.write('\n@@@\n\n')
   #Polarized PES
   curREM = ParseRems('rems/rem_polgeom')
   set_common_rems(curREM, method, basis, options)
   XYZ_to_Input(fw, XYZ, FRGM, curREM, read_mol=True, do_freq=options.do_freq, do_ferf=options.do_ferf)
   fw.write('\n@@@\n\n')
   #Frozen PES
   curREM = ParseRems('rems/rem_frzgeom')
   set_common_rems(curREM, method, basis, options)
   XYZ_to_Input(fw, XYZ, FRGM, curREM, read_mol=True, do_freq=options.do_freq)

#the script
options, args = ParseInput(sys.argv)
xyz_path = args[1]
if not options.onsite:
   frgm_partition = args[2]
if xyz_path[-1:] != '/':
   xyz_path += '/'
#determine xyz_file list
#"target" should be the most common mode here
curdir = os.getcwd()
xyzfile_list = []
if options.all:
   xyzfile_list = glob.glob(xyz_path+'*.xyz')
elif options.keyword:
   xyzfile_list = glob.glob(xyz_path+'*'+options.keyword+'*.xyz')
if options.target!=None:
   for xyz_file in options.target:
      if xyz_file not in xyzfile_list:
         xyzfile_list.append(xyz_file)

#check the input_path
input_path = options.input_path
if input_path[-1:] != '/':
   input_path += '/'
if not os.path.exists(input_path):
   sp.call(['mkdir', input_path])

#parse the fragment file
FRGM = ''
if not options.onsite:   #all xyz files share the common frgm file: parse it beforehand
   FRGM = ParseFRGM(frgm_partition)

#create the input file(s), and put it into the input_path
specify_method = False
if len(options.method)*len(options.basis) > 1:
   specify_method = True
for method in options.method:
   for basis in options.basis:
      for xyz_file in xyzfile_list:
         parsed_XYZ = XYZ(xyz_file)
         if options.onsite:
            frgm_file = xyz_file[:-4]+'.frgm'
            FRGM = ParseFRGM(frgm_file)
         if options.pestype.lower()=='integrated' or options.pestype.lower()=='all':  #integrated adaibatic EDA job
            if specify_method:
               inputfile = input_path + parsed_XYZ.NameRoot+'_adiabatic_eda_'+method+'_'+basis_abbr(basis)+'.in' 
            else:
               inputfile = input_path + parsed_XYZ.NameRoot+'_adiabatic_eda.in'
            fw = open(inputfile, 'w')
            XYZ_to_Input_Integrated(fw, parsed_XYZ, FRGM, method, basis, options)
            fw.close()

         elif options.pestype.lower()=='frz' or options.pestype.lower()=='frozen': #separated calculation on the frozen PES
            if specify_method:
               inputfile = input_path + parsed_XYZ.NameRoot+'_frzgeom_'+method+'_'+basis_abbr(basis)+'.in' 
            elif options.do_freq:
               inputfile = input_path + parsed_XYZ.NameRoot+'_frzgeom_freq.in'
            else:
               inputfile = input_path + parsed_XYZ.NameRoot+'_frzgeom.in'
            fw = open(inputfile, 'w')
            curREM = ParseRems('rems/rem_frzgeom')
            set_common_rems(curREM, method, basis, options)
            XYZ_to_Input(fw, parsed_XYZ, FRGM, curREM, do_freq=options.do_freq)
            fw.close()

         elif options.pestype.lower()=='pol' or options.pestype.lower()=='polarized': #separated calculation on the polarized PES
            if specify_method:
               if not options.do_ferf:
                  inputfile = input_path + parsed_XYZ.NameRoot+'_polgeom_'+method+'_'+basis_abbr(basis)+'.in' 
               else:
                  inputfile = input_path + parsed_XYZ.NameRoot+'_polgeom_ferf_'+method+'_'+basis_abbr(basis)+'.in' 

            elif options.do_freq:
               if not options.do_ferf:
                  inputfile = input_path + parsed_XYZ.NameRoot+'_polgeom_freq.in'
               else:
                  inputfile = input_path + parsed_XYZ.NameRoot+'_polgeom_ferf_freq.in'
            else:
               if not options.do_ferf:
                  inputfile = input_path + parsed_XYZ.NameRoot+'_polgeom.in'
               else:
                  inputfile = input_path + parsed_XYZ.NameRoot+'_polgeom_ferf.in'

            fw = open(inputfile, 'w')
            curREM = ParseRems('rems/rem_polgeom')
            set_common_rems(curREM, method, basis, options)
            XYZ_to_Input(fw, parsed_XYZ, FRGM, curREM, do_freq=options.do_freq, do_ferf=options.do_ferf)
            fw.close()
            
         elif options.pestype.lower()=='scf' or options.pestype.lower()=='full': #separated calculation on the fully relaxed PES
            if specify_method:
               inputfile = input_path + parsed_XYZ.NameRoot+'_scfgeom_'+method+'_'+basis_abbr(basis)+'.in' 
            elif options.do_freq:
               inputfile = input_path + parsed_XYZ.NameRoot+'_scfgeom_freq.in'
            else:
               inputfile = input_path + parsed_XYZ.NameRoot+'_scfgeom.in'
            fw = open(inputfile, 'w')
            curREM = ParseRems('rems/rem_scfgeom')
            set_common_rems(curREM, method, basis, options)
            XYZ_to_Input(fw, parsed_XYZ, FRGM, curREM, do_freq=options.do_freq)
            fw.close()

         else:
            print "Unsupported PES type"
            sys.exit(0)
