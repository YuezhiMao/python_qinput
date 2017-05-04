import os, glob, re, sys
import subprocess as sp
from optparse import OptionParser

def ParseInput(ArgsIn):
   UseMsg = "python parse_binding_energy.py [options] [result_dir]"
   parser = OptionParser(usage=UseMsg)

   parser.add_option('--number_only',dest='numonly',action='store_true',default=False,help='use the serial number as the name of the job')
   parser.add_option('-a','--all',dest='all',action='store_true',default=False, help='parse all the output directories under the specified directory')
   parser.add_option('-t','--target',dest='target',action='callback',callback=string_sp_callback,type='string',default=None,help='the output directories to parse (separated by \",\")')
   parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='parse all the output dirs containing the keyword')
   parser.add_option('-n','--nfrgm',dest='nfrgm',action='store',type='int',default=2,help='Number of fragments')
   parser.add_option('--mm_only',dest='mm_only',action='callback',callback=string_sp_callback,type='string',default=None,help='flag MM-only fragments')
   options, args = parser.parse_args(ArgsIn)
   if len(args) < 2:
      if options.all or options.keyword:
         print "You must specify a result folder for -a and -k mode"
         parser.print_help()
         sys.exit(0)

   if not options.all and options.target==None and options.keyword==None:
      print "The target directory must be specified: one or some or all"
      parser.print_help()
      sys.exit(1)

   return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def Parse_Single_Dir(options, nfrag):
   tmpfile = 'energy.tmp'
   #Assume all fragments are with QM first
   os.system('grep \"ion met\" *.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   AllEnergy = {}
   jobname = ''
   count = 0 #count for each job's energy items
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      #print name
      energy = float(l[2])
      name_parser = re.search('(\S+)_([^_]+)_([^_]+).out', name)
      if not options.numonly:
         jobname = name_parser.group(1)
      else:
         jobname = re.search('(\d+)',name_parser.group(1)).group(1)
      if jobname not in AllEnergy:
         AllEnergy[jobname]={}
         #count = 0
      method = name_parser.group(2),name_parser.group(3)
      if method not in AllEnergy[jobname]:
         AllEnergy[jobname][method] = {}
         count = 0
      AllEnergy[jobname][method][count] = energy
      count += 1
      if count > (nfrag + 1):
         print "I got too many (%d) energies for job %s" %(count, jobname)
         sys.exit(0)
   fr.close()
   sp.call(['rm', tmpfile])

   #for the pure MM fragments, we need to redo their energies
   if options.mm_only!=None:
      for mm_only_frgm in options.mm_only: 
         njobs = nfrag + 1
         placeholder1 = '\"Job ' + mm_only_frgm + ' of ' + str(njobs)+'\"' 
         placeholder2 = '\"Job ' + str(int(mm_only_frgm)+1) + ' of ' + str(njobs)+'\"' 
         os.system('grep -A 10000 '+placeholder1+' *.out | grep -B 200 '+placeholder2+' | grep \"TOTAL EFP ENERGY\" > '+tmpfile)
         fr = open(tmpfile, 'r')
         for line in fr.readlines():
            l = line.split()
            name = l[0][:-1]
            mm_energy = float(l[4])
            name_parser = re.search('(\S+)_([^_]+)_([^_]+).out', name)
            if not options.numonly:
               jobname = name_parser.group(1)
            else:
               jobname = re.search('(\d+)',name_parser.group(1)).group(1)
            if jobname not in AllEnergy:
               print "Jobname %s hasn't been created" %jobname
               sys.exit(0)
            method = name_parser.group(2),name_parser.group(3)
            if method not in AllEnergy[jobname]:
               print "Method %s/%s hasn't been created for jobname %s" %(method, jobname)
               sys.exit(0)
            AllEnergy[jobname][method][int(mm_only_frgm)-1] = mm_energy
         fr.close()
         sp.call(['rm', tmpfile])
   
   outfile = 'binding_energy.csv'
   fw = open(outfile, 'w')
   if options.numonly:
      for jobname in sorted(AllEnergy,key=lambda jobname:float(jobname)):
         compute_binding_energy(fw, AllEnergy, jobname, nfrag)
         #for method in AllEnergy[jobname]:
         #   functional, basis = method
         #   E_super = AllEnergy[jobname][method][2*nfrag]
         #   E_monomers_nocp = 0.0
         #   E_monomers_cp = 0.0
         #   for i in range(0, nfrag):
         #      E_monomers_nocp += AllEnergy[jobname][method][2*i]
         #      E_monomers_cp += AllEnergy[jobname][method][2*i+1]
         #   E_bind_nocp = E_super - E_monomers_nocp
         #   E_bind_cp = E_super - E_monomers_cp
         #   BSSE = E_bind_nocp - E_bind_cp #define BSSE to be positive
         #   fw.write("%s,%s,%s,%.4f,%.4f,%.4f\n" %(jobname,functional,basis,hart_to_kcal*E_bind_nocp, hart_to_kcal*E_bind_cp, hart_to_kcal*BSSE))
   else:
      for jobname in sorted(AllEnergy):
         compute_binding_energy(fw, AllEnergy, jobname, nfrag)
         #for method in sorted(AllEnergy[jobname]):
         #   functional, basis = method
         #   #print jobname,method
         #   E_super = AllEnergy[jobname][method][2*nfrag]
         #   E_monomers_nocp = 0.0
         #   E_monomers_cp = 0.0
         #   for i in range(0, nfrag):
         #      E_monomers_nocp += AllEnergy[jobname][method][2*i]
         #      E_monomers_cp += AllEnergy[jobname][method][2*i+1]
         #   E_bind_nocp = E_super - E_monomers_nocp
         #   E_bind_cp = E_super - E_monomers_cp
         #   BSSE = E_bind_nocp - E_bind_cp #define BSSE to be positive
         #   fw.write("%s,%s,%s,%.4f,%.4f,%.4f\n" %(jobname,functional,basis,hart_to_kcal*E_bind_nocp, hart_to_kcal*E_bind_cp, hart_to_kcal*BSSE))

   fw.close()

def compute_binding_energy(fw, AllEnergy, jobname, nfrag, doCP=False):
   hart_to_kcal = 627.5095
   for method in AllEnergy[jobname]:
      functional, basis = method
      if doCP:
         E_super = AllEnergy[jobname][method][2*nfrag]
         E_monomers_nocp = 0.0
         E_monomers_cp = 0.0
         for i in range(0, nfrag):
            E_monomers_nocp += AllEnergy[jobname][method][2*i]
            E_monomers_cp += AllEnergy[jobname][method][2*i+1]
         E_bind_nocp = E_super - E_monomers_nocp
         E_bind_cp = E_super - E_monomers_cp
         BSSE = E_bind_nocp - E_bind_cp #define BSSE to be positive
         fw.write("%s,%s,%s,%.4f,%.4f,%.4f\n" %(jobname,functional,basis,hart_to_kcal*E_bind_nocp, hart_to_kcal*E_bind_cp, hart_to_kcal*BSSE))
      else:
         E_super = AllEnergy[jobname][method][nfrag] 
         E_monomers = 0.0
         for i in range(0, nfrag):
            E_monomers += AllEnergy[jobname][method][i]
         E_bind = E_super - E_monomers
         fw.write("%s,%s,%s,%.4f\n" %(jobname,functional,basis,hart_to_kcal*E_bind))

#the main script
options, args = ParseInput(sys.argv)
nfrag = options.nfrgm
result_dir = None
if len(args) > 1:
   result_dir = args[1]
   if result_dir[-1:]!='/':
      result_dir += '/'
outdir_list = []
if options.all:
   outdir_list = glob.glob(result_dir+'*/')
elif options.keyword!=None:
   outdir_list = glob.glob(result_dir+'*'+options.keyword+'*/')
if options.target!=None:
   for target_dir in options.target:
      if target_dir not in outdir_list:
         outdir_list.append(target_dir)
print outdir_list
curdir = os.getcwd()
for outdir in outdir_list:
   os.chdir(outdir)
   Parse_Single_Dir(options,nfrag)
   os.chdir(curdir)




    


   
       
