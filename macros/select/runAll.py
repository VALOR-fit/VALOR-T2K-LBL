#! /usr/bin/env python
"""
Script to run all "official" SKMC selections for 10a and 11c together quickly.
Just calls superK_selector.C a lot of times.
Puts output into structured folders (i.e. GENIE and NEUT)
To run, just do:

./runAll.py

You can modify some variables in the script to do what you want.
"basecut" contains the cut sets to process (e.g. FV, 1RmuT, 1ReT, None)
"cutsandfiles" contains the extra cuts to apply
  There are 5 cuts to apply, and they are : separated
     cut1: e (electron neutrinos) or m (muon neutrinos)
     cut2: the actual event cut you're applying
     cut3: used as part of the output filename
     cut4: genie or neut - used to specify 27 or 0 MeV for the reconstructed energy. also puts the output in the relevant folder
     cut5: SKMC version. Used to decide which cuts to apply, based on the variables contained in the SKMC files. also puts output into relevant folder
     cut6: Write expanded tree, that can be used to reweight? t or f
     NOTE: Use a unique combination of cut3, cut4, cut5 for each call to superK_selector.C you use, otherwise overwriting will occur
  See superK_selctor.C for more details
"cutsandfiles"  also contains the files to process the cuts on.

You can turn on/off processing a given production with booleans like:
process10a = True
process11c = False
(note the capital letter!)
"""

import sys
import optparse
import os
import math # cos, log, pi etc

# will create directories relative to the current working directory, if they don't already exist
base = os.environ["PWD"]

# input path
file_loc = "/Users/davidesgalaberna/Desktop/ETH_phd/T2K/ANA/VALOR/skfiles/13a/t2ksk_13a_ntuple"

if os.environ["USER"] == "phrhbq": file_loc = "/storage/epp2/t2k/data/sk/13a/neut/t2ksk_13a_ntuple"

# output file names
outfile_numuCCQE  = "sk13a_1RnumuCCQE.root"
outfile_nueCCQE   = "sk13a_1RnueCCQE.root"
outfile_numuCC    = "sk13a_nRnumuCC.root"
outfile_nueCC     = "sk13a_nRnueCC.root"
outfile_2pi0ring  = "sk13a_2pi0ring.root"

#specify the production(s) you want to run over here
process10a = False
process11c = False
process13a = True
process14a = False

processings = []
if process10a: processings.append("10a")
if process11c: processings.append("11c")
if process13a: processings.append("13a")
if process14a: processings.append("14a")

if not os.path.isdir(base+"/selected"):
   try:
      os.makedirs(base+"/selected")
   except:
      sys.exit("Could not make selected directory")
for gen in ["/genie", "/neut"]:
   if not os.path.isdir(base+"/selected"+gen):
      try:
         os.makedirs(base+"/selected"+gen)
      except:
         sys.exit("Could not make selected"+gen+" directory")
   for mc in processings:
      if not os.path.isdir(base+"/selected"+gen+"/"+mc):
         try:
            os.makedirs(base+"/selected"+gen+"/"+mc)
         except:
            sys.exit("Could not make selected "+gen+"/"+mc+" directory")

rwght = "t"

cutsandfiles = []

if process13a:
   cutsandfiles.append(["tnumu", file_loc+"/numu.h2o.sk.flux11a.*.root","/13a/"])
   cutsandfiles.append(["tnumubar", file_loc+"/numubar.h2o.sk.flux11a.*.root","/13a/"])
   cutsandfiles.append(["tnue", file_loc+"/nue.h2o.sk.flux11a.*.root","/13a/"])
   cutsandfiles.append(["tnuebar", file_loc+"/nuebar.h2o.sk.flux11a.*.root","/13a/"])
   cutsandfiles.append(["tnuesig", file_loc+"/nue_x_numuflx.h2o.sk.flux11a.*.root","/13a/"])

#if process14a:
   #cutsandfiles.append(["e:"+basecut[b]+":nue:neut:14a:"    +rwght, file_loc+"/nue_x_nue/nue_x_nue.h2o.sk.flux13a.*.root"])
   #cutsandfiles.append(["e:"+basecut[b]+":nuebar:neut:14a:" +rwght, file_loc+"/nuebar_x_nuebar/nuebar_x_nuebar.h2o.sk.flux13a.*.root"])
   #cutsandfiles.append(["e:"+basecut[b]+":nuesig:neut:14a:" +rwght, file_loc+"/numu_x_nue/numu_x_nue.h2o.sk.flux13a.*.root"])
   #cutsandfiles.append(["m:"+basecut[b]+":numu:neut:14a:"   +rwght, file_loc+"/numu_x_numu/numu_x_numu.h2o.sk.flux13a.*.root"])
   #cutsandfiles.append(["m:"+basecut[b]+":numubar:neut:14a:"+rwght, file_loc+"/numubar_x_numubar/numubar_x_numubar.h2o.sk.flux13a.*.root"])


# Set path
path = base+"/selected"+gen+cutsandfiles[0][2]

# remove existing files
if os.path.isfile(path+outfile_nueCCQE) :
   print "Remove existing ", path+outfile_nueCCQE
   os.remove(path+outfile_nueCC)
if os.path.isfile(path+outfile_numuCCQE) :
   print "Remove existing ", path+outfile_numuCCQE
   os.remove(path+outfile_numuCC)
if os.path.isfile(path+outfile_nueCC) :
   print "Remove existing ", path+outfile_nueCC
   os.remove(path+outfile_nueCC)
if os.path.isfile(path+outfile_numuCC) :
   print "Remove existing ", path+outfile_numuCC
   os.remove(path+outfile_numuCC)
if os.path.isfile(path+outfile_2pi0ring) :
   print "Remove existing ", path+outfile_2pi0ring
   os.remove(path+outfile_2pi0ring)


#sys.exit()

for i in xrange(len(cutsandfiles)):
   args = [cutsandfiles[i][1], cutsandfiles[i][0],outfile_numuCCQE,outfile_nueCCQE,outfile_numuCC,outfile_nueCC,outfile_2pi0ring,path]
   command = './mymultiringsel.py "'+('" "'.join(args))+'" '
   print " "
   print command
   os.system(command)
