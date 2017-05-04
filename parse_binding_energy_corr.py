import os, glob, re, sys
import subprocess as sp
from optparse import OptionParser

def ParseInput(ArgsIn):
   UseMsg = "python parse_binding_energy_corr.py [options] [result_dir]"
   parser = OptionParser(usage=UseMsg)
   parser.add_option('-a','--all',dest='all',action='store_true',default=False,help='parsing all directories under the result_dir')
   parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='parsing directories containing the keyword')
   parser.add_option('-t','--target',dest='target',action='store',type='string',default=None, help='parsing a single directory of outputs')
   parser.add_option('-n','--nfrgm',dest='nfrgm',action='store',type='int',default=2,help='Number of fragments')
   parser.add_option('-v','--verbose',dest='verbose',action='store_true',default=False,help='Generate verbose output')
   parser.add_option('--pes',dest='pes',action='store_true',default=False,help='do PES scan, save geom_param as the jobname')
   parser.add_option('--placeholder',dest='placeholder',action='store',type='string',default='dist', help='placeholder for parsing geom parameters')
   parser.add_option('--kcal',dest='kcal',action='store_true',default=False,help='use kcal/mol (default is kJ/mol)')
   parser.add_option('--order', dest='order',action='store', type='int', default=0, help='Extrapolation order:\n TQ: 0, Q5: 1, 56: 2')
   options, args = parser.parse_args(ArgsIn)

   if len(args) < 2:
      parser.print_help()
      sys.exit(0)
   if not options.all and options.keyword==None and options.target==None:
      #print "-t is the only supported mode for now. Specify the folder"
      print "select one mode: -t, -k, -a"
      parser.print_help()
      sys.exit(0)

   return options, args


def Parse_MP2(AllEnergy, options):
   nfrag = options.nfrgm
   placeholder = options.placeholder
   hart_to_kcal = 627.5095
   hart_to_kj = hart_to_kcal * 4.184

   tmpfile = 'energy.tmp'
   #MP2 total (QZ)
   os.system('grep \"MP2         total energy\" *MP2.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   count = 0
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      energy = float(l[-2])
      jobname = re.search('(\S+)_MP2.out', name).group(1)
      if options.pes:
         jobname = re.search(placeholder+'_(\S+)_MP2.out', name).group(1)
      if jobname not in AllEnergy:  #find a new system 
         AllEnergy[jobname] = {}
      if "MP2TOT_QZ" not in AllEnergy[jobname]:
         AllEnergy[jobname]["MP2TOT_QZ"] = {}
         count = 0
      AllEnergy[jobname]["MP2TOT_QZ"][count] = energy
      count += 1
      if count > (2*nfrag+1):
         print "Too many (%d) MP2 total energies (QZ) for job %s" %(count, jobname)
         sys.exit(0)
   fr.close()

   #MP2 corr (QZ)
   os.system('grep \"Total  MP2   correlation\" *MP2.out > ' + tmpfile)
   fr = open(tmpfile, 'r')
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      energy = float(l[-2])
      jobname = re.search('(\S+)_MP2.out', name).group(1)
      if options.pes:
         jobname = re.search(placeholder+'_(\S+)_MP2.out', name).group(1)
      if jobname not in AllEnergy:  #find a new system 
         AllEnergy[jobname] = {}
      if "MP2CORR_QZ" not in AllEnergy[jobname]:
         AllEnergy[jobname]["MP2CORR_QZ"] = {}
         count = 0
      AllEnergy[jobname]["MP2CORR_QZ"][count] = energy
      count += 1
      if count > (2*nfrag+1):
         print "Too many (%d) MP2 correlation energies (QZ) for job %s" %(count, jobname)
         sys.exit(0)
   fr.close()

   #SCF energy (QZ)
   for jobname in AllEnergy:
      AllEnergy[jobname]["SCF_QZ"] = {}
      for count in range(0, 2*nfrag+1):
         AllEnergy[jobname]["SCF_QZ"][count] = AllEnergy[jobname]["MP2TOT_QZ"][count] - AllEnergy[jobname]["MP2CORR_QZ"][count]

   sp.call(['rm', tmpfile])
   #return AllEnergy


def Parse_ParenT(AllEnergy, options):
   nfrag = options.nfrgm
   placeholder = options.placeholder
   hart_to_kcal = 627.5095
   hart_to_kj = hart_to_kcal * 4.184

   tmpfile = 'energy.tmp'
   #SCF energy (TZ)
   os.system('grep \"SCF energy\" *parenT.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   count = 0
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      energy = float(l[-1])
      jobname = re.search('(\S+)_parenT.out', name).group(1)
      if options.pes:
         jobname = re.search(placeholder+'_(\S+)_parenT.out', name).group(1)
      if "SCF_TZ" not in AllEnergy[jobname]:
         AllEnergy[jobname]["SCF_TZ"] = {}
         count = 0
      AllEnergy[jobname]["SCF_TZ"][count] = energy
      count += 1
      if count > (2*nfrag+1):
         print "Too many (%d) SCF energies (TZ) for job %s" %(count, jobname)
         sys.exit(0)
   fr.close() 

   #MP2 total (TZ)
   os.system('grep \"MP2 energy\" *parenT.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   count = 0
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      energy = float(l[-1])
      jobname = re.search('(\S+)_parenT.out', name).group(1)
      if options.pes:
         jobname = re.search(placeholder+'_(\S+)_parenT.out', name).group(1)
      if "MP2TOT_TZ" not in AllEnergy[jobname]:
         AllEnergy[jobname]["MP2TOT_TZ"] = {}
         count = 0
      AllEnergy[jobname]["MP2TOT_TZ"][count] = energy
      count += 1
      if count > (2*nfrag+1):
         print "Too many (%d) MP2 total energies (TZ) for job %s" %(count, jobname)
         sys.exit(0)
   fr.close() 

   #CCSD(T) total (TZ)
   os.system('grep \"CCSD(T) total\" *parenT.out > '+tmpfile)
   fr = open(tmpfile, 'r')
   count = 0
   for line in fr.readlines():
      l = line.split()
      name = l[0][:-1]
      energy = float(l[-1])
      jobname = re.search('(\S+)_parenT.out', name).group(1)
      if options.pes:
         jobname = re.search(placeholder+'_(\S+)_parenT.out', name).group(1)
      if "ParenT_TOT_TZ" not in AllEnergy[jobname]:
         AllEnergy[jobname]["ParenT_TOT_TZ"] = {}
         count = 0
      AllEnergy[jobname]["ParenT_TOT_TZ"][count] = energy
      count += 1
      if count > (2*nfrag+1):
         print "Too many (%d) CCSD(T) total energies (TZ) for job %s" %(count, jobname)
         sys.exit(0)
   fr.close() 
   
   #MP2 and CCSD(T) corr (TZ)
   for jobname in AllEnergy:
      AllEnergy[jobname]["MP2CORR_TZ"] = {}
      AllEnergy[jobname]["ParenT_CORR_TZ"] = {}
      for count in range(0, 2*nfrag+1):
         AllEnergy[jobname]["MP2CORR_TZ"][count] = AllEnergy[jobname]["MP2TOT_TZ"][count] - AllEnergy[jobname]["SCF_TZ"][count]
         AllEnergy[jobname]["ParenT_CORR_TZ"][count] = AllEnergy[jobname]["ParenT_TOT_TZ"][count] - AllEnergy[jobname]["SCF_TZ"][count]
   sp.call(['rm', tmpfile])
   #return AllEnergy


#Do the T,Q extrapolation for MP2 corr, and then add the delta correction
def evaluate_interaction(AllEnergy, options):
   nfrag = options.nfrgm
   extrap_order = options.order
   for jobname in AllEnergy:
      #HF(Q) binding energy
      E_sup_hf = AllEnergy[jobname]["SCF_QZ"][2*nfrag]
      E_mon_nocp_hf = 0.0
      E_mon_cp_hf = 0.0
      for i in range(0, nfrag):
         E_mon_nocp_hf += AllEnergy[jobname]["SCF_QZ"][2*i]
         E_mon_cp_hf += AllEnergy[jobname]["SCF_QZ"][2*i+1]
      AllEnergy[jobname]["BIND_HF_QZ(noCP)"] = E_sup_hf - E_mon_nocp_hf
      AllEnergy[jobname]["BIND_HF_QZ(CP)"] = E_sup_hf - E_mon_cp_hf

      #MP2(T, Q) correlation's contribution to binding
      E_sup_mp2corr_tz = AllEnergy[jobname]["MP2CORR_TZ"][2*nfrag]
      E_sup_mp2corr_qz = AllEnergy[jobname]["MP2CORR_QZ"][2*nfrag]
      E_mon_mp2corr_tz_nocp = 0.0
      E_mon_mp2corr_tz_cp = 0.0
      E_mon_mp2corr_qz_nocp = 0.0
      E_mon_mp2corr_qz_cp = 0.0

      E_sup_mp2corr = two_points_extrap(AllEnergy[jobname]["MP2CORR_TZ"][2*nfrag], AllEnergy[jobname]["MP2CORR_QZ"][2*nfrag], extrap_order) 
      E_mon_mp2corr_nocp = 0.0
      E_mon_mp2corr_cp = 0.0
      for i in range(0, nfrag):
         E_mon_mp2corr_nocp += two_points_extrap(AllEnergy[jobname]["MP2CORR_TZ"][2*i], AllEnergy[jobname]["MP2CORR_QZ"][2*i], extrap_order)
         E_mon_mp2corr_cp += two_points_extrap(AllEnergy[jobname]["MP2CORR_TZ"][2*i+1], AllEnergy[jobname]["MP2CORR_QZ"][2*i+1], extrap_order)
         E_mon_mp2corr_tz_nocp += AllEnergy[jobname]["MP2CORR_TZ"][2*i]
         E_mon_mp2corr_tz_cp += AllEnergy[jobname]["MP2CORR_TZ"][2*i+1]
         E_mon_mp2corr_qz_nocp += AllEnergy[jobname]["MP2CORR_QZ"][2*i]
         E_mon_mp2corr_qz_cp += AllEnergy[jobname]["MP2CORR_QZ"][2*i+1]
      AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(noCP)"] = E_sup_mp2corr - E_mon_mp2corr_nocp
      AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(CP)"] = E_sup_mp2corr - E_mon_mp2corr_cp
      AllEnergy[jobname]["BIND_MP2CORR_TZ(noCP)"] = E_sup_mp2corr_tz - E_mon_mp2corr_tz_nocp
      AllEnergy[jobname]["BIND_MP2CORR_TZ(CP)"] = E_sup_mp2corr_tz - E_mon_mp2corr_tz_cp
      AllEnergy[jobname]["BIND_MP2CORR_QZ(noCP)"] = E_sup_mp2corr_qz - E_mon_mp2corr_qz_nocp
      AllEnergy[jobname]["BIND_MP2CORR_QZ(CP)"] = E_sup_mp2corr_qz - E_mon_mp2corr_qz_cp
      
      #Extrapolation
      #if not options.Q5:
      #   AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(noCP)"] = (27.0*AllEnergy[jobname]["BIND_MP2CORR_TZ(noCP)"] - 64.0*AllEnergy[jobname]["BIND_MP2CORR_QZ(noCP)"])/(27.0-64.0)
      #   AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(CP)"] = (27.0*AllEnergy[jobname]["BIND_MP2CORR_TZ(CP)"] - 64.0*AllEnergy[jobname]["BIND_MP2CORR_QZ(CP)"])/(27.0-64.0)
      #else:
      #   AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(noCP)"] = (64.0*AllEnergy[jobname]["BIND_MP2CORR_TZ(noCP)"] - 125.0*AllEnergy[jobname]["BIND_MP2CORR_QZ(noCP)"])/(64.0-125.0)
      #   AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(CP)"] = (64.0*AllEnergy[jobname]["BIND_MP2CORR_TZ(CP)"] - 125.0*AllEnergy[jobname]["BIND_MP2CORR_QZ(CP)"])/(64.0-125.0)

      AllEnergy[jobname]["BIND_MP2(noCP)"] = AllEnergy[jobname]["BIND_HF_QZ(noCP)"] + AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(noCP)"]
      AllEnergy[jobname]["BIND_MP2(CP)"] = AllEnergy[jobname]["BIND_HF_QZ(CP)"] + AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(CP)"]
      AllEnergy[jobname]["BSSE_MP2"] = AllEnergy[jobname]["BIND_MP2(noCP)"] - AllEnergy[jobname]["BIND_MP2(CP)"]


      #CCSD(T) (TZ) correlation's contribution to binding
      E_sup_parenT_corr_tz = AllEnergy[jobname]["ParenT_CORR_TZ"][2*nfrag]
      E_mon_parenT_corr_tz_nocp = 0.0
      E_mon_parenT_corr_tz_cp = 0.0
      for i in range(0, nfrag):
         E_mon_parenT_corr_tz_nocp += AllEnergy[jobname]["ParenT_CORR_TZ"][2*i]
         E_mon_parenT_corr_tz_cp += AllEnergy[jobname]["ParenT_CORR_TZ"][2*i+1]
      AllEnergy[jobname]["BIND_ParenT_CORR_TZ(noCP)"] = E_sup_parenT_corr_tz - E_mon_parenT_corr_tz_nocp
      AllEnergy[jobname]["BIND_ParenT_CORR_TZ(CP)"] = E_sup_parenT_corr_tz - E_mon_parenT_corr_tz_cp

      #delta CCSD(T) Binding energy
      AllEnergy[jobname]["E_BIND(noCP)"] = AllEnergy[jobname]["BIND_HF_QZ(noCP)"] + AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(noCP)"] + (AllEnergy[jobname]["BIND_ParenT_CORR_TZ(noCP)"] - AllEnergy[jobname]["BIND_MP2CORR_TZ(noCP)"])
      AllEnergy[jobname]["E_BIND(CP)"] = AllEnergy[jobname]["BIND_HF_QZ(CP)"] + AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(CP)"] + (AllEnergy[jobname]["BIND_ParenT_CORR_TZ(CP)"] - AllEnergy[jobname]["BIND_MP2CORR_TZ(CP)"])
      AllEnergy[jobname]["E_BIND(avg)"] = AllEnergy[jobname]["BIND_HF_QZ(CP)"] + 0.5*(AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(CP)"]+AllEnergy[jobname]["BIND_MP2CORR_EXTRAP(noCP)"]) + (AllEnergy[jobname]["BIND_ParenT_CORR_TZ(CP)"] - AllEnergy[jobname]["BIND_MP2CORR_TZ(CP)"])
      AllEnergy[jobname]["BSSE"] = AllEnergy[jobname]["E_BIND(noCP)"] - AllEnergy[jobname]["E_BIND(CP)"] 

def Generate_Output(AllEnergy, options):
   hart_to_kcal = 627.5095
   hart_to_kj = hart_to_kcal * 4.184
   ratio = hart_to_kj
   if options.kcal:
      ratio = hart_to_kcal
   outfile = 'binding_energy.csv'   
   fw = open(outfile, 'w')
   if not options.verbose:
      fw.write("System,E_bind(noCP),E_bind(CP),E_bind(avg)\n")
      if options.pes:
         for jobname in sorted(AllEnergy, key=lambda jobname:float(jobname)):
            fw.write("%.2f,%.2f,%.2f,%.2f\n" %(float(jobname), AllEnergy[jobname]["E_BIND(noCP)"]*ratio, AllEnergy[jobname]["E_BIND(CP)"]*ratio, AllEnergy[jobname]["E_BIND(avg)"]*ratio))
      else:
         for jobname in sorted(AllEnergy):
            fw.write("%s,%.2f,%.2f,%.2f\n" %(jobname, AllEnergy[jobname]["E_BIND(noCP)"]*ratio, AllEnergy[jobname]["E_BIND(CP)"]*ratio, AllEnergy[jobname]["E_BIND(avg)"]*ratio))
   else: #also include the MP2 results
      fw.write("System,MP2(noCP),MP2(CP),BSSE(MP2),delta_T(noCP),delta_T(CP),delta_T(avg),BSSE(delta_T)\n")
      if options.pes:
         for jobname in sorted(AllEnergy, key=lambda jobname:float(jobname)):
            fw.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" %(float(jobname), AllEnergy[jobname]["BIND_MP2(noCP)"]*ratio, AllEnergy[jobname]["BIND_MP2(CP)"]*ratio,AllEnergy[jobname]["BSSE_MP2"]*ratio,AllEnergy[jobname]["E_BIND(noCP)"]*ratio,AllEnergy[jobname]["E_BIND(CP)"]*ratio,AllEnergy[jobname]["E_BIND(avg)"]*ratio,AllEnergy[jobname]["BSSE"]*ratio))
      else:
         for jobname in sorted(AllEnergy):
            fw.write("%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" %(jobname, AllEnergy[jobname]["BIND_MP2(noCP)"]*ratio, AllEnergy[jobname]["BIND_MP2(CP)"]*ratio,AllEnergy[jobname]["BSSE_MP2"]*ratio,AllEnergy[jobname]["E_BIND(noCP)"]*ratio,AllEnergy[jobname]["E_BIND(CP)"]*ratio,AllEnergy[jobname]["E_BIND(avg)"]*ratio,AllEnergy[jobname]["BSSE"]*ratio))
   fw.close()

def two_points_extrap(value_small, value_big, extrap_order=0):
   if extrap_order == 0:  #do TQ extrap
      return (27.0*value_small-64.0*value_big)/(27.0-64.0)
   elif extrap_order == 1: #do Q5 extrap
      return (64.0*value_small-125.0*value_big)/(64.0 - 125.0)
   elif extrap_order == 2: #do 56 extrap
      return (125.0*value_small - 216.0*value_big)/(125.0 - 216.0)
   else:
      print "Invalid extrap_order: %d" %extrap_order
      sys.exit(0)

def Parse_Single_Dir(options):
   AllEnergy = {}
   Parse_MP2(AllEnergy, options)
   Parse_ParenT(AllEnergy, options)
   evaluate_interaction(AllEnergy, options)
   Generate_Output(AllEnergy, options) 

#the script
options, args = ParseInput(sys.argv)
result_dir = args[1]
nfrag = options.nfrgm
if result_dir[-1:]!='/':
   result_dir += '/'
outdir_list = []
if options.all:
   outdir_list = glob.glob(result_dir+'*/')
elif options.keyword!=None:
   outdir_list = glob.glob(result_dir+'*'+options.keyword+'*/')
if options.target!=None:
   outdir_list.append(options.target)
print outdir_list
curdir = os.getcwd()
for outdir in outdir_list:
   os.chdir(outdir)
   Parse_Single_Dir(options)
   os.chdir(curdir)

 
