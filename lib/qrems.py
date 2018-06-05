import os, glob, re, sys
import subprocess as sp
import numpy as np

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

def AppendRem(fh,MyREMS):
   fh.write('$rem\n')
   for rem in range(MyREMS["nrems"]):
      fh.write(MyREMS["the_rem"][str(rem)]["name"]+'  '+MyREMS["the_rem"][str(rem)]["value"]+'\n')
   fh.write('$end\n')
   return

def copy_section_over(fw, filename):
   fr = open(filename, 'r')
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

def set_rems_common(curREM, method, basis):
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
   ModRem('MEM_TOTAL', '8000', curREM)
   ModRem('MEM_STATIC','2000', curREM)
   ModRem('SCF_GUESS', 'SAD', curREM)
   ModRem('METHOD', method, curREM)	
   if (method != 'HF'):
      ModRem ('XC_GRID', '000099000590', curREM)		
      ModRem ('NL_GRID', '1', curREM)
   #deal with -D3 in the method
   if "-D3" in method.upper():
      add_d3_tail(method, curREM)

def add_d3_tail(method, curREM):
   if method.upper() == 'WB97X-D3' or method.upper() == 'B97-D3':
      return
   l = re.search('(\S+)-(D3\S*)', method)
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
