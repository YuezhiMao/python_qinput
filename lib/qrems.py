import os, glob, re, sys
import subprocess as sp
import numpy as np

def ParseRems(RemFile, do_rem_frgm=False):
   f = open(RemFile,'r')
   REMS = {}
   REMS["name"] = RemFile[4:]
   REMS["nrems"] = 0
   currem = 0
   REMS["the_rem"] = {}
   line = f.readline()
   while not re.search('end', line):
      l=re.search("(\S+)\s+(\S+)",line) #rem and value
      if (not l==None):
         REMS["the_rem"][str(currem)] = {}
         REMS["the_rem"][str(currem)]["name"] = l.group(1).upper()
         REMS["the_rem"][str(currem)]["value"] = l.group(2).upper()
         currem += 1
         REMS["nrems"] = currem
      line = f.readline()
   if do_rem_frgm:
      #then, parse the rem_frgm section
      currem_frgm = 0
      REMS["nrems_frgm"] = 0
      REMS["rem_frgm"] = {}
      line = f.readline()
      while not re.search('end',line):
         l=re.search("(\S+)\s+(\S+)",line) #rem and value
         if (not l==None):
            REMS["rem_frgm"][str(currem_frgm)] = {}
            REMS["rem_frgm"][str(currem_frgm)]["name"] = l.group(1).upper()
            REMS["rem_frgm"][str(currem_frgm)]["value"] = l.group(2).upper()
            currem_frgm += 1
            REMS["nrems_frgm"] = currem_frgm
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
      currem_frgm = MyREMS["nrems_frgm"]
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

def AppendRemFrgm(fh,MyREMS):
   fh.write('\n$rem_frgm\n')
   for rem_frgm in range(MyREMS["nrems_frgm"]):
      fh.write(MyREMS["rem_frgm"][str(rem_frgm)]["name"]+'  '+ MyREMS["rem_frgm"][str(rem_frgm)]["value"]+'\n')
   fh.write('$end\n')
   return

def copy_section_over(fw, filename):
   fr = open(filename, 'r')
   fw.write('\n')
   for line in fr.readlines():
      fw.write(line)
   fr.close()

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
   elif fullname == 'def2-svp' or fullname == 'svp':
      return 'svp'
   elif fullname == 'def2-svpd' or fullname == 'svpd':
      return 'svpd'
   elif fullname == 'def2-tzvp' or fullname == 'tzvp':
      return 'tzvp'
   elif fullname == 'def2-tzvpp' or fullname == 'tzvpp':
      return 'tzvpp'
   elif fullname == 'def2-tzvpd' or fullname == 'tzvpd':
      return 'tzvpd'
   elif fullname == 'def2-tzvppd' or fullname == 'tzvppd':
      return 'tzvppd'
   elif fullname == 'def2-qzvp' or fullname == 'qzvp':
      return 'qzvp'
   elif fullname == 'def2-qzvpp' or fullname == 'qzvpp':
      return 'qzvpp'
   elif fullname == 'def2-qzvpd' or fullname == 'qzvpd':
      return 'qzvpd'
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
   elif fullname == 'mp':
      return 'mp'
   elif fullname == 'tp':
      return 'tp' 
   elif fullname == 'lp':
      return 'lp'
   elif fullname == '3-21g':
      return '321g'
   else:
      print("unrecognized basis choice")
      sys.exit(1)

def set_rems_common(curREM, method, basis, loose=False):
   if basis.upper() == 'SVP':
      ModRem('BASIS','DEF2-SVP',curREM)
   elif basis.upper() == 'SVPD':
      ModRem('BASIS','DEF2-SVPD',curREM)
   elif basis.upper() == 'TZVP':
      ModRem('BASIS', 'DEF2-TZVP', curREM)
   elif basis.upper() == 'TZVPP':
      ModRem('BASIS','DEF2-TZVPP',curREM)
   elif basis.upper() == 'TZVPD':
      ModRem('BASIS','DEF2-TZVPD',curREM)
   elif basis.upper() == 'TZVPPD':
      ModRem('BASIS','DEF2-TZVPPD',curREM)
   elif basis.upper() == 'QZVP':
      ModRem('BASIS','DEF2-QZVP',curREM)
   elif basis.upper() == 'QZVPP':
      ModRem('BASIS','DEF2-QZVPP',curREM)
   elif basis.upper() == 'QZVPD':
      ModRem('BASIS','DEF2-QZVPD',curREM)
   elif basis.upper() == 'QZVPPD':
      ModRem('BASIS','DEF2-QZVPPD',curREM)
   elif basis.upper() == 'DP': #double-zeta Pople: 6-31+G(d)
      ModRem('BASIS','6-31+G(d)', curREM)
   elif basis.upper() == 'MP':
      ModRem('BASIS','6-311G(d,p)',curREM)
   elif basis.upper() == 'TP': #triple-zeta Pople: 6-311++G(2df,2pd)
      ModRem('BASIS','6-311++G(2df,2pd)',curREM)
   elif basis.upper() == 'LP':
      ModRem('BASIS','6-311++G(3df,3pd)', curREM)
   else: #copy the name of standard basis
      ModRem('BASIS', basis, curREM)

   if loose:
      #ModRem('THRESH', '12', curREM)
      ModRem('MEM_TOTAL', '16000', curREM)
      ModRem('MEM_STATIC','8000', curREM)
   else:
      ModRem('THRESH', '14', curREM)
      ModRem('MEM_TOTAL', '8000', curREM)
      ModRem('MEM_STATIC','2000', curREM)
   ModRem('SCF_GUESS', 'SAD', curREM)
   ModRem('METHOD', method, curREM)	
   if (method != 'HF'):
      if loose:
         ModRem('XC_GRID', '000075000302', curREM)
      else:
         ModRem('XC_GRID', '000099000590', curREM)
      ModRem ('NL_GRID', '1', curREM)
   #deal with -D3 in the method
   if "-D3" in method.upper():
      add_d3_tail(method, curREM)

def add_d3_tail(method, curREM):
   if method.upper() == 'WB97X-D3' or method.upper() == 'B97-D3':
      return
   l = re.search('(\S+)-(D3\S*)', method.upper())
   parent_func = l.group(1)
   tail = l.group(2)
   ModRem('METHOD', parent_func, curREM)
   if 'BJ' in tail.upper():
      ModRem('DFT_D', 'D3_BJ', curREM)
   elif 'OP' in tail.upper():
      ModRem('DFT_D', 'D3_OP', curREM)
   elif 'CSO' in tail.upper():
      ModRem('DFT_D', 'D3_CSO', curREM)
   else:
      ModRem('DFT_D', 'EMPIRICAL_GRIMME3', curREM)

def add_aux_basis(basis, curREM):
   if basis.upper() == 'CC-PVTZ':
      ModRem('AUX_BASIS', 'RIMP2-CC-PVTZ', curREM)
   elif basis.upper() == 'AUG-CC-PVTZ':
      ModRem('AUX_BASIS', 'RIMP2-AUG-CC-PVTZ', curREM)
   elif basis.upper() == 'CC-PVDZ':
      ModRem('AUX_BASIS', 'RIMP2-CC-PVDZ', curREM)
   elif basis.upper() == 'AUG-CC-PVDZ':
      ModRem('AUX_BASIS', 'RIMP2-AUG-CC-PVDZ', curREM)
   elif basis.upper() == 'CC-PVQZ':
      ModRem('AUX_BASIS', 'RIMP2-CC-PVQZ', curREM)
   elif basis.upper() == 'AUG-CC-PVQZ':
      ModRem('AUX_BASIS', 'RIMP2-AUG-CC-PVQZ', curREM)
   else:
      print("No corresponding auxiliary basis for "+basis+" yet")

def apply_single_geom_constraint(fw, geom_param, constraint_template):
   fr = open(constraint_template, 'r')
   fw.write('\n$opt\n')
   for line in fr.readlines():
      l_split = line.split()
      if len(l_split) > 2: #constraint line
         l_split[-1] = str(geom_param)
         for idx, item in enumerate(l_split, 1):
            if idx != len(l_split):
               fw.write('%s\t' %item)
            else:
               fw.write('%s\n' %item)
      else:
         fw.write(line)
   fw.write('$end\n')

def AppendSolvationSecs(fw, pcm_epsilon=0.0, smd_solvent=None):
   if pcm_epsilon > 0.0:
      fw.write('\n$pcm\n')
      fw.write('Theory  CPCM\n')
      fw.write('Method  SWIG\n')
      fw.write('Solver  INVERSION\n') 
      fw.write('HPoints     302\n')
      fw.write('HeavyPoints 302\n')
      fw.write('$end\n\n')
      fw.write('$solvent\n')
      fw.write('Dielectric  %.2f\n' %pcm_epsilon)
      fw.write('$end\n')
   elif smd_solvent != None:
      fw.write('\n$smx\n')
      fw.write('solvent  %s\n' %smd_solvent)
      fw.write('$end\n')
   else:
      print("why here!?")
      sys.exit(0)
