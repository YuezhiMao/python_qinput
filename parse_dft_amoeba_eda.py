import os, glob, re, sys
import subprocess as sp
from optparse import OptionParser

def ParseInput(ArgsIn):
   UseMsg = "python parse_binding_energy.py [options] [result_dir]"
   parser = OptionParser(usage=UseMsg)

   parser.add_option('--number_only',dest='numonly',action='store_true',default=False,help='use the serial number as the name of the job')
   parser.add_option('--pes', dest='pes', action='store_true',default=False,help='use a certain geometry parameter as the jobname')
   parser.add_option('-a','--all',dest='all',action='store_true',default=False, help='parse all the output directories under the specified directory')
   parser.add_option('-t','--target',dest='target',action='callback',callback=string_sp_callback,type='string',default=None,help='the output directories to parse (separated by \",\")')
   parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='parse all the output dirs containing the keyword')
   parser.add_option('--blur',dest='blur',action='store_true',default=False,help='For calculations with Gaussian blurring in QM/AMOEBA')
   parser.add_option('-e','--embed',dest='embed',action='store_true',default=False,help='Parse coulomb embedding jobs')

   options, args = parser.parse_args(ArgsIn)
   if len(args) < 2:
      if options.all or options.keyword:
         print "You must specify a result folder for -a or -k mode" 
         parser.print_help()
         sys.exit(0)

   if not options.all and options.target==None and options.keyword==None:
      print "The target directory must be specified: one or some or all"
      parser.print_help()
      sys.exit(1)

   return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def determine_jobname(name, options):
   if options.numonly:
      name_parser = re.search('(\S+)_([^_]+)_([^_]+).out', name)
      jobname = re.search('(\d+)',name_parser.group(1)).group(1)
   elif options.pes:
      placeholder = ''
      if 'dist' in name:
         placeholder = 'dist'
      elif 'angle' in name:
         placeholder = 'angle'
      else:
         print "Can't find a valid placeholder"
         sys.exit(1)
      name_parser = re.search(placeholder+'_([^_]+)_', name)
      jobname = name_parser.group(1) 
      #print "jobname: %s" %jobname
   else:
      name_parser = re.search('(\S+)_([^_]+)_([^_]+).out', name)
      jobname = name_parser.group(1)
   return jobname


def Parse_Single_Dir(options):
   hart_to_kcal = 627.5095

   #parse the energy items
   tmpfile = 'energy.tmp'
   #parse isolated AMOEBA ENERGY
   os.system('grep -B 1000 \"Job 2 of 3\" *.out | grep \"TOTAL EFP ENERGY\" > '+tmpfile)
   fr = open(tmpfile, 'r')
   #create the entries based on this grep
   AllEnergy = {}
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      #print name
      jobname = determine_jobname(name, options)
      if jobname not in AllEnergy:
         AllEnergy[jobname]={}
         AllEnergy[jobname]["nitem"] = 0
      energy = float(l[4]) #amoeba energy in a.u.
      AllEnergy[jobname]["MM_energy"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()

   #parse isolated QM Energy
   os.system('grep -A 10000 \"Job 2 of 3\" *.out | grep -B 300 \"Job 3 of 3\" | grep \"ion met\" > '+tmpfile)
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      jobname = determine_jobname(name, options)
      if jobname not in AllEnergy:
         print "The isolated EFP calculation for job %s has failed horribly" %jobname
         sys.exit(0)
      energy = float(l[2]) #qm energy in a.u.
      AllEnergy[jobname]["QM_energy"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()

   #parse the total energy of QM/MM system
   if options.blur:
      os.system('grep \"SCF energy (excluding charge-charge)\" *.out > '+tmpfile)
   else:
      os.system('grep -A 5000 \"Job 3 of 3\" *.out | grep \"SCF   energy in the final basis set\" > '+tmpfile) 
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      jobname = determine_jobname(name, options)
      if AllEnergy[jobname]["nitem"]!=2:
         print "Job %s has %d items before: something has failed before" %(jobname, AllEnergy[jobname]["nitem"])
         sys.exit(0)
      energy = float(l[-1])
      AllEnergy[jobname]["QMMM_energy"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()

   #parse QM/MM permanent electrostatics
   os.system('grep \"QM/AMOEBA PERM\" *.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      jobname = determine_jobname(name, options)
      if AllEnergy[jobname]["nitem"]!=3:
         print "Job %s has %d items before: something has failed before" %(jobname, AllEnergy[jobname]["nitem"])
         sys.exit(0)
      energy = float(l[4])
      AllEnergy[jobname]["QMMM_ELEC"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()

   #parse QM/MM vdW
   os.system('grep \"QM/AMOEBA VDW\" *.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      jobname = determine_jobname(name, options)
      if AllEnergy[jobname]["nitem"]!=4:
         print "Job %s has %d items before: something has failed before" %(jobname, AllEnergy[jobname]["nitem"])
         sys.exit(0)
      energy = float(l[4])
      AllEnergy[jobname]["QMMM_VDW"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()
   sp.call(['rm', tmpfile])

   #output
   outfile = 'qmmm_interaction.csv'
   fw = open(outfile, 'w')
   if options.numonly or options.pes:
      for jobname in sorted(AllEnergy,key=lambda jobname:float(jobname)):
         #print jobname
         E_int = hart_to_kcal*(AllEnergy[jobname]["QMMM_energy"] - AllEnergy[jobname]["QM_energy"] - AllEnergy[jobname]["MM_energy"])
         E_elec = hart_to_kcal*AllEnergy[jobname]["QMMM_ELEC"]
         E_vdw = hart_to_kcal*AllEnergy[jobname]["QMMM_VDW"]
         E_pol = E_int - E_elec - E_vdw
         fw.write("%s,%.4f,%.4f,%.4f,%.4f\n" %(jobname, E_int, E_elec, E_pol, E_vdw))
   else:
      for jobname in sorted(AllEnergy):
         E_int = hart_to_kcal*(AllEnergy[jobname]["QMMM_energy"] - AllEnergy[jobname]["QM_energy"] - AllEnergy[jobname]["MM_energy"])
         E_elec = hart_to_kcal*AllEnergy[jobname]["QMMM_ELEC"]
         E_vdw = hart_to_kcal*AllEnergy[jobname]["QMMM_VDW"]
         E_pol = E_int - E_elec - E_vdw
         fw.write("%s,%.4f,%.4f,%.4f,%.4f\n" %(jobname, E_int, E_elec, E_pol, E_vdw))
   fw.close()

def Parse_Single_Dir_embedding(options):
   hart_to_kcal = 627.5095
   tmpfile = 'energy.tmp'
   AllEnergy = {}
   #parse permanent electrostatics
   os.system('grep \"EelecEmbed\" *.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = re.search('(\S+.out)', l[0]).group(1)
      jobname = determine_jobname(name, options)
      if jobname not in AllEnergy:
         AllEnergy[jobname] = {}
         AllEnergy[jobname]["nitem"] = 0
      energy_str = l[-2]
      if energy_str[0] == '(':
         energy_str = energy_str[1:]
      energy = float(energy_str)  #the second last entry is the perm_elec energy in kcal/mol 
      AllEnergy[jobname]["PERM_ELEC"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()
   sp.call(['rm', tmpfile])

   #parse polarization (forward only)
   os.system('grep \"Energies,\" *.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = re.search('(\S+.out)', l[0]).group(1)
      jobname = determine_jobname(name, options)
      if jobname not in AllEnergy:
         print "Somehow the perm elec energy of job %s is missing" %jobname
         sys.exit(0)
      energy_str = l[-2]
      if energy_str[0] == '(':
         energy_str = energy_str[1:]
      energy = float(energy_str)
      AllEnergy[jobname]["POL"] = energy
      AllEnergy[jobname]["nitem"] += 1
   fr.close()
   sp.call(['rm', tmpfile])
   
   #output
   if options.embed:
      outfile = 'coulomb_embedding.csv'
   else:
      outfile = 'qmmm_interaction.csv'
   fw = open(outfile, 'w')
   if options.numonly or options.pes:
      for jobname in sorted(AllEnergy,key=lambda jobname:float(jobname)):
         fw.write("%s,%.4f,%.4f\n" %(jobname, AllEnergy[jobname]["PERM_ELEC"], AllEnergy[jobname]["POL"]))
   else:
      for jobname in sorted(AllEnergy):        
         fw.write("%s,%.4f,%.4f\n" %(jobname, AllEnergy[jobname]["PERM_ELEC"], AllEnergy[jobname]["POL"]))
   fw.close()


#the main script
options, args = ParseInput(sys.argv)
curdir = os.getcwd()
result_dir = None
if len(args) > 1:
   result_dir = args[1]
   if result_dir[-1:]!='/':
      result_dir += '/'
#get the outdir list
outdir_list = []
if options.all:
   outdir_list = glob.glob(result_dir+'*/')
elif options.keyword!=None:
   outdir_list = glob.glob(result_dir+'*'+options.keyword+'*/')
if options.target!=None:
   for target_dir in options.target:
      print target_dir
      if target_dir not in outdir_list:
         outdir_list.append(target_dir)
print outdir_list
#parse them
for outdir in outdir_list:
   os.chdir(outdir)
   if options.embed:
      Parse_Single_Dir_embedding(options)
   else:
      Parse_Single_Dir(options)
   os.chdir(curdir)
