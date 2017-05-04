import os, glob, re, sys
import subprocess as sp
from optparse import OptionParser

def ParseInput(ArgsIn):
   UseMsg = "python edajob_parser.py [options] [result_dir]\nExample for a result dir: jobs_new_eda_opt3_atz/result/"
   parser = OptionParser(usage=UseMsg)
   parser.add_option('-a','--all',dest='all',action='store_true',default=False, help='parse all the output directories under the specified directory')
   parser.add_option('--skip',dest='skip',action='store',type='string',default=None,help='The output directory we want to skip')
   parser.add_option('--skip_word',dest='skip_word',action='store',type='string',default=None,help='skip the output dir(s) containing the specified keyword')
   parser.add_option('-t','--target',dest='target',action='callback',callback=string_sp_callback,type='string',default=None,help='the output directories to parse (separated by \",\")')
   parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='parse all the output dirs containing the keyword')
   parser.add_option('--noCT',dest='noCT',action='store_true',default=False,help='Do not parse the CT related results (required by the EDA_noCT jobs)')
   parser.add_option('--fullavail',dest='fullavail',action='store_true',default=False,help='full SCF data for evaluating CT (for noCT jobs only)')
   parser.add_option('--mp2',dest='mp2',action='store_true',default=False,help='Parse the total binding energy (only) evaluated by MP2 or double hybrid DFT (T-Q extrapolation presumed') #this should ignore a lot of other options
   parser.add_option('--single_atom',dest='single_atom',action='store',default=None, help="When parsing the total binding curve, enable the special case in which one monomer is a single atom so one job covers everything. Recognizable options: mon1(monA,1), mon2(monB,2), both")
   parser.add_option('--old',dest='oldeda',action='store_true',default=False,help='Parsing EDA jobs using the old ALMO-EDA code')
   parser.add_option('--modified',dest='modified',action='store_true',default=False,help='Modified frozen decomposition: E_cls_elec, E_mod_pauli, E_disp')
   parser.add_option('--ff',dest='ff',action='store_true',default=False,help='Three terms only: cls_elec, pol, vdw(pauli+disp+ct)')
   parser.add_option('--pes',dest='pes',action='store_true',default=False,help='scanning a PES (searching for \"dist\" or \"angle\"')
   parser.add_option('--snapshot',dest='snapshot',action='store_true',default=False,help='doing EDA for snapshots instead of scanning a potential energy curve')
   parser.add_option('--snap_interval',dest='snap_interval',action='store',type='int',default=1, help='Only collecting data for snapshot numbers that are multiples of this number')
   parser.add_option('--mbe',dest='mbe',action='store_true',default=False,help='Doing many-body expansion (3-body only)')
   parser.add_option('--placeholder',dest='placeholder',action='store',default='dist', type='string',help='placeholder for pes scan (the word in front of the geom parameter)')
   parser.add_option('--kcal',dest='kcal',action='store_true',default=False,help='use kcal/mol for printed interactions (only applied to new EDA jobs)')
   parser.add_option('--modelchem',dest='modelchem',action='store_true',default=False,help='parse the model chemistry used for the job')
   parser.add_option('--bsse',dest='bsse',action='store_true',default=False,help='including the CP-based BSSE manually')
   
   options, args=parser.parse_args(ArgsIn)
   if options.mp2:     #when parsing MP2, we should prevent other options being triggered
      options.neweda = False
      options.noCT = False
      options.fullavail = False
   if options.single_atom!=None:
      option_list = ["1","2","mon1","mon2","monA","monB", "both"]
      if options.single_atom not in option_list:
         print "Single-atom option unrecognized"
         parser.print_help()
         sys.exit(1)
      #regularize the input: mon1, mon2, both
      elif options.single_atom == "1" or options.single_atom=="mon1" or options.single_atom=="monA":
         options.single_atom = "mon1"
      elif options.single_atom == "2" or options.single_atom=="mon2" or options.single_atom=="monA":
         options.single_atom = "mon2"
      
   if len(args) < 2:
      if options.all or options.keyword:  #in these two scenarios
         print "You must specify a result folder for -a and -k mode"
         parser.print_help()
         sys.exit(1)
   if not options.all and options.target==None and options.keyword==None:
      print "The target directory must be specified: one or some or all"
      parser.print_help()
      sys.exit(1)
   
   else:
      return options,args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def parse_oldeda(options, output_dir):
#currently parsed: E_frz, E_pol, E_vct; computed quantities: E_bind  
#for convenience, it has simply assumed all the jobs are finished successfully
   hart_to_kj = 2625.5311584660003
   curdir = os.getcwd()
   os.chdir(output_dir)
   #shortname is how the directory is called under "result/"
   shortname = re.search('/([^/]+)/$',output_dir).group(1)
   placeholder = ''   #little trick for parsing
   if re.search('dist',output_dir):
      placeholder = 'dist'
   elif re.search('angle',output_dir):
      placeholder = 'angle'
   else:
      print "Neither \"dist\" nor \"angle\", the placeholder got confused"
      sys.exit(1)
   
   DataPoints = {}
   outlist = glob.glob('*.out')
   #outlist=sorted(outlist, key=str)
   #print outlist
   for out in outlist:
      #print out
      name_parser = re.search(placeholder+'_([^_]+)_', out)
      geom_param = name_parser.group(1)
      #print geom_param
      DataPoints[geom_param] ={}
   tmp_filename = 'grep.tmp'
   
   #parse E_frz
   os.system('grep \"E_frz\" *.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   for line in f.readlines():
      l = re.search(placeholder+'_([^_]+)_\S+:(\S+).+\s+(\S+)$', line)
      geom_param = l.group(1)
      item = l.group(2)
      value = float(l.group(3))
      if geom_param not in DataPoints:
         print "geom_param: %s is not in DataPoints dictionary" %geom_param
         sys.exit(1)
      DataPoints[geom_param][item] = value
   f.close()

   #parse E_pol
   os.system('grep \"E_pol\" *.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   for line in f.readlines():
      l = re.search(placeholder+'_([^_]+)_\S+:(\S+).+\s+(\S+)$', line)
      geom_param = l.group(1)
      item = l.group(2)
      value = float(l.group(3))
      if geom_param not in DataPoints:
         print "geom_param: %s is not in DataPoints dictionary" %geom_param
         sys.exit(1)
      DataPoints[geom_param][item] = value
   f.close()
   
   #parse E(pfrz)
   os.system('grep \"E(Pfrz)\" *.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   for line in f.readlines():
      l = re.search(placeholder+'_([^_]+)_\S+:(\S+).+\s+(\S+)$', line)
      geom_param = l.group(1)
      item = l.group(2)
      value = float(l.group(3))
      if geom_param not in DataPoints:
         print "geom_param: %s is not in DataPoints dictionary" %geom_param
         sys.exit(1)
      DataPoints[geom_param][item] = value
   f.close()

   #parse fragment energies
   os.system('grep -A 2 -i \"Fragment Energies\" *.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   frag_count = 0
   for line in f.readlines():
      if not re.search(placeholder, line):
         continue
      new_point = re.search(placeholder+'_([^_]+)_\S+:Fragment Energies', line)
      if new_point!=None:
         frag_count = 0
         continue
      
      l = re.search(placeholder+'_([^_]+)_\S+\s+(\S+)$',line)
      if l!=None:
         frag_count += 1
         if frag_count > 2:
            print "Warning: this script only prints the energies of the first two fragments"
         item = "E_frag"+str(frag_count)
         geom_param = l.group(1)
         value = float(l.group(2))
         if geom_param not in DataPoints:
            print "geom_param: %s is not in DataPoints dictionary" %geom_param
            sys.exit(1)
         DataPoints[geom_param][item] = value
   f.close() 
   
      
   if not options.noCT: 
      #parse E_vct
      os.system('grep \"E_vct\" *.out > '+tmp_filename)
      f = open(tmp_filename, 'r')
      for line in f.readlines():
         l = re.search(placeholder+'_([^_]+)_\S+:(\S+).+\s+(\S+)$', line)
         geom_param = l.group(1)
         item = l.group(2)
         value = float(l.group(3))
         if geom_param not in DataPoints:
            print "geom_param: %s is not in DataPoints dictionary" %geom_param
            sys.exit(1)
         DataPoints[geom_param][item] = value
      f.close()
      #compute E_bind
      for geom_param in DataPoints:
         DataPoints[geom_param]["E_bind"] = DataPoints[geom_param]["E_frz"]+DataPoints[geom_param]["E_pol"]+DataPoints[geom_param]["E_vct"]
   
   elif options.fullavail:
      # we need two numbers: pol SCF energy (from the EDA job) and full SCF energy (from the normal scf job)
      os.system('grep \"ion met\" *eda*.out > '+tmp_filename)
      f = open(tmp_filename, 'r')
      for line in f.readlines():
         l = re.search(placeholder+'_([^_]+)_\S+:\s+\S+\s+(\S+)',line) 
         geom_param = l.group(1)
         item = "pol_scf"
         value = float(l.group(2))
         if geom_param not in DataPoints:
            print "geom_param: %s is not in DataPoints dictionary" %geom_param
            sys.exit(1)
         DataPoints[geom_param][item] = value
      f.close()
      
      os.system('grep \"ion met\" *normalscf*.out > '+tmp_filename)
      f = open(tmp_filename, 'r')
      for line in f.readlines():
         l = re.search(placeholder+'_([^_]+)_\S+:\s+\S+\s+(\S+)',line) 
         geom_param = l.group(1)
         item = "full_scf"
         value = float(l.group(2))
         if geom_param not in DataPoints:
            print "geom_param: %s is not in DataPoints dictionary" %geom_param
            sys.exit(1)
         DataPoints[geom_param][item] = value
      f.close()

      #E_ct = E(full_scf) - E(pol_scf)
      for geom_param in DataPoints:
         DataPoints[geom_param]["E_vct"] = hart_to_kj*(DataPoints[geom_param]["full_scf"] - DataPoints[geom_param]["pol_scf"])
         DataPoints[geom_param]["E_bind"] = DataPoints[geom_param]["E_frz"]+DataPoints[geom_param]["E_pol"]+DataPoints[geom_param]["E_vct"]
          
   #generate output
   csvfile = shortname+'.csv'
   if not options.noCT:   #with CT
      f = open(csvfile,'w')
      placeholder1 = placeholder
      placeholder2 = "E_frz (kJ/mol)"
      placeholder3 = "E_pol (kJ/mol)"
      placeholder4 = "E_ct (kJ/mol)"
      placeholder5 = "E_bind (kJ/mol)"
      placeholder6 = "E(Pfrz) (Ha)"
      placeholder7 = "E_frag1 (Ha)"
      placeholder8 = "E_frag2 (Ha)"
      f.write("%15s %15s %15s %15s %15s %15s %15s %15s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8))
      for geom_param in sorted(DataPoints,key=lambda geom_param:float(geom_param)):
         #print geom_param
         f.write("%15s %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"],DataPoints[geom_param]["E_vct"],DataPoints[geom_param]["E_bind"],DataPoints[geom_param]["E(Pfrz)"],DataPoints[geom_param]["E_frag1"],DataPoints[geom_param]["E_frag2"]))
      f.close()
   elif options.fullavail:
      f = open(csvfile,'w')
      placeholder1 = placeholder
      placeholder2 = "E_frz (kJ/mol)"
      placeholder3 = "E_pol (kJ/mol)"
      placeholder4 = "E_ct (kJ/mol)"
      placeholder5 = "E_bind (kJ/mol)"
      f.write("%15s %15s %15s %15s %15s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5))
      for geom_param in sorted(DataPoints,key=lambda geom_param:float(geom_param)):
         #print geom_param
         f.write("%15s %15.8f %15.8f %15.8f %15.8f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"],DataPoints[geom_param]["E_vct"],DataPoints[geom_param]["E_bind"]))
      f.close()

   else:  #no CT and full SCF not available
      f = open(csvfile,'w')
      placeholder1 = placeholder
      placeholder2 = "E_frz (kJ/mol)"
      placeholder3 = "E_pol (kJ/mol)"
      f.write("%15s %15s %15s\n" %(placeholder1, placeholder2, placeholder3))
      for geom_param in sorted(DataPoints,key=lambda geom_param:float(geom_param)):
         #print geom_param
         f.write("%15s %15.8f %15.8f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"]))
      f.close()
   sp.call(['rm',tmp_filename])
   os.chdir(curdir)
   return


#for parsing a PES scan
def parse_neweda_pes(options, output_dir, keyword=None): 
   if options.noCT:
      print "parse_disp_decomp doesn't support the noCT option for now"
      sys.exit(1)

   hart_to_kj = 2625.5311584660003
   hart_to_kcal = 627.5095
   curdir = os.getcwd()
   os.chdir(output_dir)
   #shortname = re.search('/([^/]+).*$',output_dir).group(1)
   placeholder = ''  #little trick to get the geom_param
   if re.search('dist',output_dir):
      placeholder = 'dist'
   elif re.search('angle',output_dir):
      placeholder = 'angle'
   elif options.placeholder!=None:
      placeholder = options.placeholder
   else:
      print "Neither \"dist\" nor \"angle\", the placeholder got confused"
      sys.exit(1)

   DataPoints = {} 
   outlist = glob.glob('*.out')
   for out in outlist:
      name_parser = re.search(placeholder+'_([^_]+)_', out)
      geom_param = name_parser.group(1)
      DataPoints[geom_param] ={}
   tmp_filename = 'grep.tmp'

   if not keyword:
      os.system('grep \"kJ/mol\" *.out > '+tmp_filename)
   else:
      os.system('grep \"kJ/mol\" *'+keyword+'*.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   for line in f.readlines():
      l = re.search(placeholder+'_([^_]+)_\S+:\s+(\S+).+=\s+(\S+)',line)
      if l!=None:
         geom_param = l.group(1)
         item = l.group(2)
         value = float(l.group(3))
         if geom_param not in DataPoints:
            print "geom_param: %s is not in DataPoints dictionary. Something is very wrong" %geom_param
            sys.exit(1)
         DataPoints[geom_param][item] = value
         if options.kcal:
            DataPoints[geom_param][item] *= 1.0/hart_to_kj * hart_to_kcal #final results in kcal/mol
   f.close()

   #including the BSSE if desired
   if options.bsse:
      for geom_param in DataPoints:
         DataPoints[geom_param]["E_int"] += DataPoints[system]["BSSE"]
         DataPoints[geom_param]["E_vct"] += DataPoints[system]["BSSE"]

   #generate output
   if not keyword:
      csvfile_basic = 'EDA2.csv'
   else:
      csvfile_basic = 'EDA2_'+keyword+'.csv'

   f = open(csvfile_basic,'w')
   placeholder1 = 'distance'
   if placeholder == 'angle':
      placeholder1 = 'angle'
   placeholder2 = "E_frz"
   placeholder3 = "E_pol"
   placeholder4 = "E_vct"
   placeholder5 = "E_int"
   if not options.modified:
      placeholder6 = "E_elec"
      placeholder7 = "E_pauli"
      placeholder8 = "E_disp"
      placeholder9 = "E_cls_elec"
      f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8, placeholder9))
   else:
      placeholder6 = "E_cls_elec"
      placeholder7 = "E_mod_pauli"
      placeholder8 = "E_disp" 
      f.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8))

   for geom_param in sorted(DataPoints,key=lambda geom_param:float(geom_param)):
      if not options.modified:
         f.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"],DataPoints[geom_param]["E_vct"],DataPoints[geom_param]["E_int"],DataPoints[geom_param]["E_elec"],DataPoints[geom_param]["E_pauli"],DataPoints[geom_param]["E_disp"],DataPoints[geom_param]["E_cls_elec"]))
      else:
         f.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"],DataPoints[geom_param]["E_vct"],DataPoints[geom_param]["E_int"],DataPoints[geom_param]["E_cls_elec"],DataPoints[geom_param]["E_frz"]-DataPoints[geom_param]["E_cls_elec"]-DataPoints[geom_param]["E_disp"],DataPoints[geom_param]["E_disp"]))

   f.close()

   sp.call(['rm', tmp_filename])
   os.chdir(curdir)
   return
 

#parse generic EDA jobs (using full jobname to specify each data point)
def parse_neweda_generic(options, output_dir, keyword=None): #TODO: make this the update-to-date parser for new EDA jobs
   if options.noCT:
      print "parse_disp_decomp doesn't support the noCT option for now"
      sys.exit(1)

   hart_to_kj = 2625.5311584660003
   hart_to_kcal = 627.5095
   curdir = os.getcwd()
   os.chdir(output_dir)
   #shortname = re.search('/([^/]+).*$',output_dir).group(1)

   DataPoints = {} 
   outlist = glob.glob('*.out')
   tmp_filename = 'grep.tmp'

   #grep all the items at one time
   if not keyword:
      os.system('grep \"kJ/mol\" *.out > '+tmp_filename)
   else:
      os.system('grep \"kJ/mol\" *'+keyword+'*.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   for line in f.readlines():
      l = re.search('(\S+)_eda2_op(\d)_(\S+).out:\s+(\S+).+=\s+(\S+)',line)
      if l!=None:
         system = l.group(1)
         eda_option = 'op'+l.group(2)
         method = l.group(3)
         item = l.group(4)
         value = float(l.group(5))
         if system not in DataPoints:
            DataPoints[system]={}
         if not options.modelchem:
            DataPoints[system][item] = value
            if options.kcal:
               DataPoints[system][item] *= 1.0/hart_to_kj * hart_to_kcal
         else:
            if method not in DataPoints[system]:
               DataPoints[system][method] = {}
            DataPoints[system][method][item] = value
            if options.kcal:
               DataPoints[system][method][item] *= 1.0/hart_to_kj * hart_to_kcal
   f.close()
   #including the BSSE if desired
   if options.bsse:
      for system in DataPoints:
         DataPoints[system]["E_int"] += DataPoints[system]["BSSE"]
         DataPoints[system]["E_vct"] += DataPoints[system]["BSSE"]

   #generate output
   if not keyword:
      csvfile_basic = 'EDA2.csv'
   else:
      csvfile_basic = 'EDA2_'+keyword+'.csv'
   f = open(csvfile_basic,'w')
   if options.ff:
      for system in sorted(DataPoints):
         if not options.modelchem:
            E_int = DataPoints[system]["E_int"]
            E_perm = DataPoints[system]["E_cls_elec"]
            E_pol = DataPoints[system]["E_pol"]
            E_vdw = DataPoints[system]["E_cls_pauli"] + DataPoints[system]["E_vct"] 
            f.write("%s,%.4f,%.4f,%.4f,%.4f\n" %(system, E_int, E_perm, E_pol, E_vdw))
         else:
            for method in DataPoints[system]:
               E_int = DataPoints[system][method]["E_int"]
               E_perm = DataPoints[system][method]["E_cls_elec"]
               E_pol = DataPoints[system][method]["E_pol"]
               E_vdw = DataPoints[system][method]["E_cls_pauli"] + DataPoints[system][method]["E_vct"] 
               f.write("%s,%s,%.4f,%.4f,%.4f,%.4f\n" %(system, method, E_int, E_perm, E_pol, E_vdw))

   else:
      placeholder1 = 'system'
      placeholder2 = "E_frz"
      placeholder3 = "E_pol"
      placeholder4 = "E_vct"
      placeholder5 = "E_int"
      if not options.modified:
         placeholder6 = "E_elec"
         placeholder7 = "E_pauli"
         placeholder8 = "E_disp"
         placeholder9 = "E_cls_elec"
         f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8, placeholder9))
      else:
         placeholder6 = "E_cls_elec"
         placeholder7 = "E_mod_pauli"
         placeholder8 = "E_disp" 
         f.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8))
      for system in sorted(DataPoints):
         if not options.modelchem:
            E_frz = DataPoints[system]["E_frz"] 
            E_pol = DataPoints[system]["E_pol"]
            E_vct = DataPoints[system]["E_vct"]
            E_int = DataPoints[system]["E_int"]
            E_elec = DataPoints[system]["E_elec"]
            E_pauli = DataPoints[system]["E_pauli"]
            E_disp = DataPoints[system]["E_disp"]
            E_cls_elec = DataPoints[system]["E_cls_elec"]
            if not options.modified:
               f.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(system,E_frz,E_pol,E_vct,E_int,E_elec,E_pauli,E_disp,E_cls_elec))
            else:
               E_mod_pauli = E_frz - E_disp - E_cls_elec
               f.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(system,E_frz,E_pol,E_vct,E_int,E_cls_elec,E_mod_pauli,E_disp))
         else:
            for method in DataPoints[system]:
               E_frz = DataPoints[system][method]["E_frz"] 
               E_pol = DataPoints[system][method]["E_pol"]
               E_vct = DataPoints[system][method]["E_vct"]
               E_int = DataPoints[system][method]["E_int"]
               E_elec = DataPoints[system][method]["E_elec"]
               E_pauli = DataPoints[system][method]["E_pauli"]
               E_disp = DataPoints[system][method]["E_disp"]
               E_cls_elec = DataPoints[system][method]["E_cls_elec"]
               if not options.modified:
                  f.write("%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(system,method,E_frz,E_pol,E_vct,E_int,E_elec,E_pauli,E_disp,E_cls_elec))
               else:
                  E_mod_pauli = E_frz - E_disp - E_cls_elec
                  f.write("%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(system,method,E_frz,E_pol,E_vct,E_int,E_cls_elec,E_mod_pauli,E_disp))

   f.close()

   sp.call(['rm', tmp_filename])
   os.chdir(curdir)
   return


#parse EDA results for snapshots from MD simulations
def parse_neweda_snapshot(options, output_dir, keyword=None): 
   if options.noCT:
      print "parse_disp_decomp doesn't support the noCT option for now"
      sys.exit(1)

   hart_to_kj = 2625.5311584660003
   hart_to_kcal = 627.5095
   curdir = os.getcwd()
   os.chdir(output_dir)

   DataPoints = {} 

   #grep all the items at one time
   tmp_filename = 'grep.tmp'
   if not keyword:
      os.system('grep \"kJ/mol\" *.out > '+tmp_filename)
   else:
      os.system('grep \"kJ/mol\" *'+keyword+'*.out > '+tmp_filename)
   f = open(tmp_filename, 'r')
   for line in f.readlines():
      l = line.split()
      jobname = l[0][:-1]
      snapID = re.search('(\d+)_eda2', jobname).group(1)
      if snapID not in DataPoints:
         DataPoints[snapID]={}
      #make sure it is an energy item
      if re.search('\d+', l[-1])!=None:
         item = l[1]
         energy = float(l[-1])
         DataPoints[snapID][item] = energy
         if options.kcal:
            DataPoints[snapID][item] *= 1.0/hart_to_kj * hart_to_kcal
   f.close()
   #including the BSSE if desired
   if options.bsse:
      for sanpID in DataPoints:
         DataPoints[snapID]["E_int"] += DataPoints[system]["BSSE"]
         DataPoints[snapID]["E_vct"] += DataPoints[system]["BSSE"]

   #generate output
   if not keyword:
      csvfile_basic = 'EDA2.csv'
   else:
      csvfile_basic = 'EDA2_'+keyword+'.csv'

   #the "basic" file
   f = open(csvfile_basic,'w')
   placeholder1 = 'ID'
   placeholder2 = "E_frz"
   placeholder3 = "E_pol"
   placeholder4 = "E_ct"
   placeholder5 = "E_int"
   if not options.modified:
      placeholder6 = "E_elec"
      placeholder7 = "E_pauli"
      placeholder8 = "E_disp"
      placeholder9 = "E_cls_elec"
      f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8, placeholder9))
   else:
      placeholder6 = "E_cls_elec"
      placeholder7 = "E_mod_pauli"
      placeholder8 = "E_disp" 
      f.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %(placeholder1, placeholder2, placeholder3, placeholder4, placeholder5, placeholder6, placeholder7, placeholder8))

   for geom_param in sorted(DataPoints,key=lambda geom_param:float(geom_param)):
      if not options.modified:
         f.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"],DataPoints[geom_param]["E_vct"],DataPoints[geom_param]["E_int"],DataPoints[geom_param]["E_elec"],DataPoints[geom_param]["E_pauli"],DataPoints[geom_param]["E_disp"],DataPoints[geom_param]["E_cls_elec"]))
      else:
         f.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(geom_param,DataPoints[geom_param]["E_frz"],DataPoints[geom_param]["E_pol"],DataPoints[geom_param]["E_vct"],DataPoints[geom_param]["E_int"],DataPoints[geom_param]["E_cls_elec"],DataPoints[geom_param]["E_frz"]-DataPoints[geom_param]["E_cls_elec"]-DataPoints[geom_param]["E_disp"],DataPoints[geom_param]["E_disp"]))

   f.close()

   sp.call(['rm', tmp_filename])
   os.chdir(curdir)
   return

def compute_3b_terms(outdir):
   curdir = os.getcwd()
   os.chdir(outdir)
   data_tot = np.genfromtxt('EDA2.csv', dtype=None, delimiter=',', skip_header=1)
   data_mbe12 = np.genfromtxt('EDA2_mbe12.csv', dtype=None, delimiter=',', skip_header=1)
   data_mbe13 = np.genfromtxt('EDA2_mbe13.csv', dtype=None, delimiter=',', skip_header=1)
   data_mbe23 = np.genfromtxt('EDA2_mbe23.csv', dtype=None, delimiter=',', skip_header=1)
   data_3b = np.zeros(data_tot.shape)
   #copy the first column
   data_3b[:,0] = data_tot[:,0]
   data_3b[:,1:] = data_tot[:,1:] - data_mbe12[:,1:] - data_mbe13[:,1:] - data_mbe23[:,1:]
   np.savetxt("3b_term.csv", data_3b, delimiter=',', fmt='%.2f')
   os.chdir(curdir)

def parse_mp2(options, output_dir):
   hart_to_kj = 2625.5311584660003
   curdir = os.getcwd()
   os.chdir(output_dir)
   #shortname is how the directory is called under "result/"
   shortname = re.search('/([^/]+)/$',output_dir).group(1)
   placeholder = ''   #little trick for parsing
   if re.search('dist',output_dir):
      placeholder = 'dist'
   elif re.search('angle',output_dir):
      placeholder = 'angle'
   else:
      print "Neither \"dist\" nor \"angle\", the placeholder got confused"
      sys.exit(1)
   
   DataPoints = {}
   outlist = glob.glob('*.out')
   #print outlist
   for out in outlist:
      #Apparently we don't want double initialization, so skip all the monomer outputs
      if re.search("mon", out):
         continue
      #print out
      name_parser = re.search(placeholder+'_([^_]+)_', out)
      geom_param = name_parser.group(1)
      #print geom_param
      DataPoints[geom_param] ={}
      #initialize all the three pieces
      DataPoints[geom_param]["complex"]={}
      DataPoints[geom_param]["mon1"]={}
      DataPoints[geom_param]["mon2"]={}
   tmp_filename = 'grep.tmp'

   #parse TZ and QZ SCF energies
   #complex first
   os.system('grep \"ion met\" *.out | grep -v \"mon\" > '+ tmp_filename)
   f = open(tmp_filename, 'r')
   geom_param = None
   scfE_counts = 0
   for line in f.readlines():
      line_parser = line.split()
      l = re.search(placeholder+'_([^_]+)_', line_parser[0])
      value = float(line_parser[2])
      if l.group(1)!=geom_param: # a new geom_param, reset scfE_counts and geom_param
         geom_param = l.group(1)
         scfE_counts = 0
      if scfE_counts == 0: #Triple zeta goes first
         DataPoints[geom_param]["complex"]["SCF_TZ"] = value
      elif scfE_counts == 1:
         DataPoints[geom_param]["complex"]["SCF_QZ"] = value
      else:
         print "There are more than two SCF energies around for the same job. Something is wrong!"
         sys.exit(1)
      scfE_counts += 1
   f.close()
   #then two monomers' SCF energies
   #collect data first, then decide whether to use the same data for all the datapoints
   same_mon1 = False
   same_mon2 = False
   if options.single_atom == 'mon1' or options.single_atom == 'both': 
      same_mon1 = True
   elif options.single_atom == 'mon2' or options.single_atom == 'both':
      same_mon2 = True
   
   #mon1 and mon2
   os.system('grep \"ion met\" *.out | grep \"mon\" > '+ tmp_filename)
   f = open(tmp_filename, 'r')
   geom_param = None
   monomer_tag = None
   scfE_counts = 0
   # for the "same_mon" case
   cur_E_TZ_mon1 = 0.0
   cur_E_QZ_mon1 = 0.0
   cur_E_TZ_mon2 = 0.0
   cur_E_QZ_mon2 = 0.0
   for line in f.readlines():
      line_parser = line.split()
      l = re.search(placeholder+'_([^_]+)_(mon\d+)', line_parser[0])
      value = float(line_parser[2])
      if l.group(1)!=geom_param or l.group(2)!=monomer_tag: # a new geom_param, reset scfE_counts and geom_param
         geom_param = l.group(1)
         monomer_tag = l.group(2)
         scfE_counts = 0
      if scfE_counts == 0: #Triple zeta goes first
         DataPoints[geom_param][monomer_tag]["SCF_TZ"] = value
         if monomer_tag == "mon1":
            cur_E_TZ_mon1 = value
         else:
            cur_E_TZ_mon2 = value
      elif scfE_counts == 1:
         DataPoints[geom_param][monomer_tag]["SCF_QZ"] = value
         if monomer_tag == "mon1":
            cur_E_QZ_mon1 = value
         else:
            cur_E_QZ_mon2 = value
      else:
         print "There are more than two SCF energies around for the same job. Something is wrong!"
         sys.exit(1)
      scfE_counts += 1
   f.close()

   if same_mon1:
      for geom_param in DataPoints:
         DataPoints[geom_param]["mon1"]["SCF_TZ"] = cur_E_TZ_mon1
         DataPoints[geom_param]["mon1"]["SCF_QZ"] = cur_E_QZ_mon1
   elif same_mon2:
      for geom_param in DataPoints:
         DataPoints[geom_param]["mon2"]["SCF_TZ"] = cur_E_TZ_mon2
         DataPoints[geom_param]["mon2"]["SCF_QZ"] = cur_E_QZ_mon2
   
   #parse TZ and QZ MP2 correlation energies
   #complex first
   os.system('grep \"Total  MP2   correlation energy\" *.out | grep -v \"mon\" > '+tmp_filename)
   f = open(tmp_filename, 'r')
   geom_param = None
   corrE_counts = 0
   for line in f.readlines():
      line_parser = line.split()
      l = re.search(placeholder+'_([^_]+)_', line_parser[0])
      value = float(line_parser[6])  #MP2 correlation energy
      if l.group(1)!=geom_param:
         geom_param = l.group(1)  
         corrE_counts = 0
      if corrE_counts == 0:  #TZ
         DataPoints[geom_param]["complex"]["corr_TZ"] = value
      elif corrE_counts == 1:
         DataPoints[geom_param]["complex"]["corr_QZ"] = value
      else:
         print "There are more than two correlation energies around for the same job. Something is wrong!"
         sys.exit(1)
      corrE_counts += 1
   f.close()
   #then two monomers
   os.system('grep \"Total  MP2   correlation energy\"  *.out | grep \"mon\" > '+ tmp_filename)
   f = open(tmp_filename, 'r')
   geom_param = None
   monomer_tag = None
   corrE_counts = 0
   # for the "same_mon" case
   cur_corrE_TZ_mon1 = 0.0
   cur_corrE_QZ_mon1 = 0.0
   cur_corrE_TZ_mon2 = 0.0
   cur_corrE_QZ_mon2 = 0.0
   for line in f.readlines():
      line_parser = line.split()
      l = re.search(placeholder+'_([^_]+)_(mon\d+)', line_parser[0])
      value = float(line_parser[6])
      if l.group(1)!=geom_param or l.group(2)!=monomer_tag: # a new geom_param, reset scfE_counts and geom_param
         geom_param = l.group(1)
         monomer_tag = l.group(2)
         corrE_counts = 0
      if corrE_counts == 0: #Triple zeta goes first
         DataPoints[geom_param][monomer_tag]["corr_TZ"] = value
         if monomer_tag == "mon1":
            cur_corrE_TZ_mon1 = value
         else:
            cur_corrE_TZ_mon2 = value
      elif corrE_counts == 1:
         DataPoints[geom_param][monomer_tag]["corr_QZ"] = value
         if monomer_tag == "mon1":
            cur_corrE_QZ_mon1 = value
         else:
            cur_corrE_QZ_mon2 = value
      else:
         print "There are more than two correlation energies around for the same job. Something is wrong!"
         sys.exit(1)
      corrE_counts += 1
   f.close()
   if same_mon1:
      for geom_param in DataPoints:
         DataPoints[geom_param]["mon1"]["corr_TZ"] = cur_corrE_TZ_mon1
         DataPoints[geom_param]["mon1"]["corr_QZ"] = cur_corrE_QZ_mon1
   elif same_mon2:
      for geom_param in DataPoints:
         DataPoints[geom_param]["mon2"]["corr_TZ"] = cur_corrE_TZ_mon2
         DataPoints[geom_param]["mon2"]["corr_QZ"] = cur_corrE_QZ_mon2
   #do the extrapolation, and compute the binding energy
   for geom_param in DataPoints:
      for item in DataPoints[geom_param]:
         DataPoints[geom_param][item]["corr_extrap"] = (3.0**3*DataPoints[geom_param][item]["corr_TZ"]-4.0**3*DataPoints[geom_param][item]["corr_QZ"])/(3.0**3-4.0**3)
      DataPoints[geom_param]["E_bind"] = (DataPoints[geom_param]["complex"]["SCF_QZ"]+DataPoints[geom_param]["complex"]["corr_extrap"])-(DataPoints[geom_param]["mon1"]["SCF_QZ"]+DataPoints[geom_param]["mon1"]["corr_extrap"])-(DataPoints[geom_param]["mon2"]["SCF_QZ"]+DataPoints[geom_param]["mon2"]["corr_extrap"])
      DataPoints[geom_param]["E_bind"] *= hart_to_kj 
   #generate output
   csvfile =shortname+'_binding.csv'
   f = open(csvfile,'w')
   placeholder1 = placeholder
   placeholder2 = "E_bind(KJ/mol)"
   f.write("%15s %15s\n" %(placeholder1, placeholder2))
   for geom_param in sorted(DataPoints,key=lambda geom_param:float(geom_param)):
      f.write("%15s %15.8f\n" %(geom_param, DataPoints[geom_param]["E_bind"]))
   f.close()
   sp.call(['rm', tmp_filename])
   os.chdir(curdir)
   return 


#The script
options, args = ParseInput(sys.argv)
outdir_list = []
cur_dir = os.getcwd()
result_dir = None
if len(args) > 1:
   result_dir = args[1]
   if result_dir[-1:]!='/':
      result_dir += '/'

#get the outdir_list
if options.all:
   outdir_list = glob.glob(result_dir+'*/')
elif options.keyword!=None:
   outdir_list = glob.glob(result_dir+'*'+options.keyword+'*/')
if options.target!=None:
   for target_dir in options.target:
      if target_dir not in outdir_list:
         outdir_list.append(target_dir)
 
if options.skip!=None:
   skip_dir = options.skip
   if skip_dir[-1:]!='/':
      skip_dir += '/'
   outdir_list.remove(skip_dir)
if options.skip_word!=None:
   skip_list = glob.glob(result_dir+'*'+options.skip_word+'*/')
   for skip_dir in skip_list:
      outdir_list.remove(skip_dir)
print outdir_list

#parse them
for outdir in outdir_list:
   if options.mp2:
      parse_mp2(options, outdir)
   elif options.oldeda:
      parse_oldeda(options, outdir)
   else:
      if options.pes:
         parse_neweda_pes(options, outdir)
      elif options.snapshot:
         parse_neweda_snapshot(options, outdir)
      else:
         parse_neweda_generic(options, outdir)
      #3-body MBE
      if options.mbe:
         mbe_list = ['mbe12','mbe23','mbe13']
         for keyword in mbe_list:
            if options.pes:
               parse_neweda_pes(options, outdir, keyword)
            elif options.snapshot:
               parse_neweda_snapshot(options, outdir, keyword)
            else:
               parse_neweda_generic(options, outdir, keyword)
         compute_3b_terms(outdir)
