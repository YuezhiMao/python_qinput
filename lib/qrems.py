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
