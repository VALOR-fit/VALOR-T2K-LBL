#!/usr/bin/env python

import ROOT
import math
import os
import glob
import numpy as n
import getopt
import sys

ROOT.gROOT.SetBatch(1)

from sys import exit
from math import sqrt, sin
from ROOT import TTree
from ROOT import *
from glob import glob
from ROOT import TLorentzVector
from ROOT import TVector3

dir_nu  = [0.669764, -0.742179, 0.024223]; #t2k beam direction
neutronMass  = 939.56536
protonMass   = 938.27203
pionMass     = 139.570 # charged
deltaMass    = 1232.
electronMass = 0.510999
muonMass     = 105.66

binding_energy = 27.

kGAMMA    = 0 
kELECTRON = 1 
kMUON     = 2 
kPION     = 3 
kKAON     = 4 
kPROTON   = 5 
kCONE     = 6

pNames = dict( [ (kGAMMA,"photon"),
                 (kELECTRON,"electron"),
                 (kMUON,"muon"),
                 (kPION,"pion"),
                 (kKAON,"kaon"),
                 (kPROTON,"proton"),
                 (kCONE,"cone")
               ])

# correct for neutrino mode only
neut_codes = dict( [ (1 , "CCQE"           ),
                     (2 , "MEC"            ),
                     (11, "CC1Pi+"         ),
                     (12, "CC1Pi0"         ),
                     (13, "CC1Pi-"         ),
                     (16, "CCCOH1Pi+"      ),
                     (17, "CCgamma"        ),
                     (21, "CCnPi"          ),
                     (22, "CC1Eta0"        ),
                     (23, "CC1K+"          ),
                     (26, "CCDIS"          ),
                     (31, "NC1PiO"         ),
                     (32, "NC1Pi0"         ),
                     (33, "NC1Pi-"         ),
                     (34, "NC1Pi+"         ),
                     (36, "NCCOH1Pi0"      ),
                     (38, "NCgamma"        ),
                     (39, "NCgamma"        ),
                     (41, "NCnPi"          ),
                     (42, "NCEta0"         ),
                     (43, "NCEta0"         ),
                     (44, "NCK0"         ),
                     (45, "NCK+"         ),
                     (46, "NCDIS"          ),
            ])
            
try:
  ROOT.t2kana
except:
  ROOT.gROOT.ProcessLine(".L $T2KOSCANA/root/root2k.C");
  pass

#~ class FitQun1RFit:
  #~ def __init__(self,tree):
    #~ self.fqnse = tree.fnqse
    #~ 
    #~ self.status = []
    
    
class FitQunMRFit:
  def __init__(self, tree, fit_index):
    self.fit_index  = fit_index
    self.fit_id     = tree.fqmrifit  [fit_index]
    self.nring      = tree.fqmrnring [fit_index]
    self.status     = tree.fqmrpcflg [fit_index]
    self.nll        = tree.fqmrnll   [fit_index]
    
    self.ring_pid       = []
    self.ring_momentum  = []
    self.ring_vertex    = []
    self.ring_direction = []
    
    for iring in xrange(self.nring):
       self.ring_pid      .append(tree.fqmrpid[6*fit_index+iring])
       self.ring_momentum .append(tree.fqmrmom[6*fit_index+iring])
       self.ring_vertex   .append(TVector3( tree.fqmrpos[18*fit_index+3*iring+0],
                                            tree.fqmrpos[18*fit_index+3*iring+1],
                                            tree.fqmrpos[18*fit_index+3*iring+2]) )
       self.ring_direction.append(TVector3( tree.fqmrdir[18*fit_index+3*iring+0],
                                            tree.fqmrdir[18*fit_index+3*iring+1],
                                            tree.fqmrdir[18*fit_index+3*iring+2]) )
                                            
    return;
    
  def GetNMuons(self):
    return len([pid for pid in self.ring_pid if pid==kMUON])
        
  def GetNPions(self):
    return len([pid for pid in self.ring_pid if pid==kPION])
    
  def GetNElectrons(self):
    return len([pid for pid in self.ring_pid if pid==kELECTRON])
    
  def AsString(self):
      s = `self.fit_id` + "\t; "
      s+=", ".join(pNames[pid] for pid in self.ring_pid)
      s+=" ; nll="+`self.nll`
      return s
        
#_______________________________________________________________________
def main(argv):
  
  files = glob(argv[0])
  treename = argv[1]
  outfile_numuCCQE = argv[2]
  outfile_nueCCQE  = argv[3]
  outfile_numuCC   = argv[4]
  outfile_nueCC    = argv[5]
  outfile_2pi0ring = argv[6]
  path = argv[7]+"/"
  
  # set path
  outfile_numuCCQE = path+outfile_numuCCQE
  outfile_nueCCQE  = path+outfile_nueCCQE
  outfile_numuCC   = path+outfile_numuCC
  outfile_nueCC    = path+outfile_nueCC
  outfile_2pi0ring = path+outfile_2pi0ring
  
  # output tree (share same variables since they are mututally exclusive)
  tree_numuccqe = ROOT.TTree(treename, treename)
  tree_nueccqe  = ROOT.TTree(treename, treename)
  tree_numucc   = ROOT.TTree(treename, treename)
  tree_nuecc    = ROOT.TTree(treename, treename)
  tree_2rpi0    = ROOT.TTree(treename, treename)

  global tmode
  tmode= n.zeros(1, dtype=n.int)
  global terec
  terec = n.zeros(1, dtype=n.float32)
  global tenu
  tenu = n.zeros(1, dtype=n.float32)
  global twgt
  twgt = n.zeros(1, dtype=n.float32)
  global tplep
  tplep = n.zeros(1, dtype=n.float32)
  global tcostheta
  tcostheta = n.zeros(1, dtype=n.float32)
  global tQ2
  tQ2 = n.zeros(1, dtype=n.float32)
  global tevis
  tevis = n.zeros(1, dtype=n.float32)
  global tnring
  tnring = n.zeros(1, dtype=n.int)

  global terec_delta
  terec_delta = n.zeros(1,dtype=n.float32)
  global terec_pi
  terec_pi = n.zeros(1,dtype=n.float32)
  
  output_trees = [tree_numuccqe,tree_nueccqe,tree_numucc,tree_nuecc,tree_2rpi0]
  
  for tree in output_trees:
      tree.Branch('mode', tmode, 'mode/I') 
      tree.Branch('erec', terec, 'erec/F')
      tree.Branch('enu', tenu, 'enu/F')
      tree.Branch('wgt', twgt, 'wgt/F')
      tree.Branch('plep', tplep, 'plep/F')
      tree.Branch('costheta', tcostheta, 'costheta/F')
      tree.Branch('Q2', tQ2, 'Q2/F')
      tree.Branch('evis', tevis, 'evis/F')
      tree.Branch('nring', tnring, 'nring/I')
      tree.Branch('erec_delta', terec_delta, 'erec_delta/F')
      tree.Branch('erec_pi', terec_pi, 'erec_pi/F')
  
  tree = TChain("h1")
  for f in files: tree.Add(f)

  h_nueCCQEname  = "h_nueCCQE" +treename
  h_numuCCQEname = "h_numuCCQE"+treename
  h_nueCCname    = "h_nueCC" +treename
  h_numuCCname   = "h_numuCC"+treename
  h_2rpi0name    = "h_2rpi0" +treename

  h_nueCCQE  = ROOT.TH1F(h_nueCCQEname ,"h_nueCCQE;;# of events" ,10,0,10)
  h_numuCCQE = ROOT.TH1F(h_numuCCQEname,"h_numuCCQE;;# of events",10,0,10)
  h_nueCC    = ROOT.TH1F(h_nueCCname   ,"h_nueCC;;# of events"   ,10,0,10)
  h_numuCC   = ROOT.TH1F(h_numuCCname  ,"h_numuCC;;# of events"  ,10,0,10)
  h_2rpi0    = ROOT.TH1F(h_2rpi0name   ,"h_2rpi0;;# of events"   ,10,0,10)
  
  hists = [h_nueCCQE , h_numuCCQE , h_nueCC , h_numuCC , h_2rpi0 ]
  
  nev = tree.GetEntries()
  
  if not nev:
    print "No entries found (",treename,")"
    return
    
  fv_nev   = 0
  fcfv_nev = 0
    
  for iev,event in enumerate(tree):
    if iev%10000==0 : print "Event", iev," of ",nev," ("+treename+")"
    
    # apply cuts (all samples are mutually exclusive)
    if select_fv(tree)==False : continue # fiducial volume events
    fv_nev+=1
    if select_fc(tree)==False : continue # fiducial volume events
    fcfv_nev+=1
    
    #~ mr_fits = [FitQunMRFit(tree,i) for i in xrange(tree.fqnmrfit) if tree.fqmrifit[i]/100000000>0]
    #~ mr_fits.sort(key=lambda x: x.nll)
    #~ 
    #~ try:
      #~ print neut_codes[tree.mode]
    #~ except KeyError:
      #~ print tree.mode
    #~ 
    #~ for fit in mr_fits:
      #~ print fit.AsString()
      #~ 
    #~ print "\n===\n"
    
    
    if select_numuccqe(tree):
      isnumu = True
      FillVarTree(tree,isnumu)
      tree_numuccqe.Fill()
      h_numuCCQE.Fill(1)
      continue
        
    if select_nueccqe(tree):
      isnumu = False
      FillVarTree(tree,isnumu)
      tree_nueccqe.Fill()
      h_nueCCQE.Fill(1)
      continue
    
    #~ if reject1ring(tree)==False : continue # reject 1ring events

    if select_pi0(tree,False)==True : # 2-ring pi0
      isnumu = False
      if ErecoCut(tree,isnumu)==False : continue 
      FillVarTree(tree,isnumu)
      tree_2rpi0.Fill()
      h_2rpi0.Fill(1)
      continue
    
    if select_mr_numucc(tree)==True : # multi-ring numucc (highest mom ring)
      isnumu = True
      if ErecoCut(tree,isnumu)==False : continue 
      FillVarTree(tree,isnumu)
      tree_numucc.Fill()
      h_numuCC.Fill(1)
      continue

    if select_mr_nuecc(tree)==True : # multi-ring nuecc (highest mom ring)
      isnumu = False
      if ErecoCut(tree,isnumu)==False : continue 
      FillVarTree(tree,isnumu)
      tree_nuecc.Fill()
      h_nueCC.Fill(1)
      continue


  # Fill scaling histograms
  try:
    scale_name = tree_names[treename]
    scale = sk_event_counts[scale_name]
  except KeyError:
    print "Cannot find scaling for ",treename," scalings will be wrong"
    scale = 1.
    
  for h in hists:
    h.Fill(1.5,1E21)
    h.Fill(0.5, scale / fv_nev)
    h.Fill(2.5, scale / fcfv_nev)
    h.Fill(3.5, scale / nev)
  # write trees in files
  
  fout_numuccqe = ROOT.TFile(outfile_numuCCQE, "update")
  tree_numuccqe.Write()
  h_numuCCQE.Write()
  fout_numuccqe.Close()

  fout_nueccqe = ROOT.TFile(outfile_nueCCQE, "update")
  tree_nueccqe.Write()
  h_nueCCQE.Write()
  fout_nueccqe.Close()
  
  fout_numucc = ROOT.TFile(outfile_numuCC, "update")
  tree_numucc.Write()
  h_numuCC.Write()
  fout_numucc.Close()

  fout_nuecc = ROOT.TFile(outfile_nueCC, "update")
  tree_nuecc.Write()
  h_nueCC.Write()
  fout_nuecc.Close()
  
  fout_2rpi0 = ROOT.TFile(outfile_2pi0ring, "update")
  tree_2rpi0.Write()
  h_2rpi0.Write()
  fout_2rpi0.Close()
  
  
  print "Done! (",treename,"), selected:"
  print "NuMu CCQE: " , tree_numuccqe.GetEntries()
  print "NuE CCQE:  " , tree_nueccqe.GetEntries()
  print "NuMu CC:   " , tree_numucc.GetEntries()
  print "NuE CC:    " , tree_nuecc.GetEntries()
  print "2Rpi0:     " , tree_2rpi0.GetEntries()
  
###################################### Selection ###############################

def select_fv(tree):
  if tree.evis <= 30.: return False
  if tree.wall <= 200.: return False
  #record_cut("fv")
  return True
  
def select_fc(tree):
  if tree.nhitac >= 16: return False
  if tree.evclass != 1: return False
  return True
  
def true_fv(tree):
  if tree.mode == 0: return False
  return (tree.wallv > 200)
  
def distanceToWall(position):
  s = 1690.-math.sqrt(position[0]**2 + position[1]**2)
  t = 1810.-abs(position[2])
  
  return min(s,t)
  
def fitQun1RingPID(tree, subevent=0):
  t = tree
  ind = 7*subevent
  if -t.fq1rnll[ind+kELECTRON]+t.fq1rnll[ind+kMUON] > t.fq1rmom[ind+kELECTRON]*0.2:
    return kELECTRON
  elif (-t.fq1rnll[ind+kPION] + -t.fq1rnll[ind+kMUON] > t.fq1reloss[ind+kPION]*0.35-3.):
    return kPION
  
  return kMUON
  
def select_numuccqe(tree):
  t = tree
  
  if not select_fc(t) : return False     # FC cut
  
  if t.fqnse > 2: return False # 1 ring events need 2 or fewer subevents
  
  subevent=0
  mu_ind = subevent*7+kMUON
  mu_vec_ind = subevent*7*3 + kMUON*3
  
  if t.fq1rpcflg[mu_ind] < 0: return False #MINUIT convergence check
  
  vertex_position = [t.fq1rpos[mu_vec_ind+i] for i in xrange(3)]
  if not distanceToWall(vertex_position) > 200.: return False # FV cut
   
  pid_1Ring = fitQun1RingPID(t)
  if not (pid_1Ring==kMUON):        return False # muon PID cut
  if not (t.fq1rmom[kMUON]>200.):   return False # muon momentum cut
  
  
  mr_fits = [FitQunMRFit(tree,i) for i in xrange(tree.fqnmrfit) if abs(tree.fqmrifit[i]/100000000)==2]
  mr_fits.sort(key=lambda x: x.nll)
  
  if len(mr_fits):
    best_mr_nll = mr_fits[0].nll
    if not (-best_mr_nll+t.fq1rnll[mu_ind] < t.fq1rmom[mu_ind]*0.25 + 12.5): return False #multiring likelihood cut
  
  return True
    
def select_nueccqe(tree):
  t = tree
  if not select_fc(t) : return False     # FC cut
  
  if t.fqnse > 2: return False # 1 ring events need 2 or fewer subevents
  
  if not tree.evis > 100.      : return False # evis cut, should be translated to FQ?
  
  subevent=0
  elec_ind = subevent*7+kELECTRON
  elec_vec_ind = subevent*7*3 + kELECTRON*3
  
  if t.fq1rpcflg[elec_ind] < 0: return False #MINUIT convergence check
  
  vertex_position = [t.fq1rpos[elec_vec_ind+i] for i in xrange(3)]
  if not distanceToWall(vertex_position) > 200.: return False # FV cut
   
  pid_1Ring = fitQun1RingPID(t)
  if not (pid_1Ring==kELECTRON):        return False # electron PID cut
  
  
  mr_fits = [FitQunMRFit(tree,i) for i in xrange(tree.fqnmrfit) if abs(tree.fqmrifit[i]/100000000)==2]
  mr_fits.sort(key=lambda x: x.nll)
  
  if len(mr_fits):
    best_mr_nll = mr_fits[0].nll
    if not (-best_mr_nll+t.fq1rnll[elec_ind] < t.fq1rmom[elec_ind]*0.25 + 12.5): return False #multiring likelihood cut
  
  if not EnuReco(tree,False) < 1250.: return False
  
  #~ fitqun_electron_likelihood = tree.fq1rnll[1]
  #~ fitqun_pizero_likelihood   = tree.fq1rnll[3]
  #~ fitqun_pizero_mass         = tree.fqpi0mass[0]
  #~ 
  #~ if not (fitqun_electron_likelihood-fitqun_pizero_likelihood) < (175 - 0.875*fitqun_pizero_mass): return False
  
  return True
  

def select_mr_numucc(tree):
  if not select_fc(tree)           : return False
  if not select_fv(tree)           : return False
  
  mr_fits = [FitQunMRFit(tree,i) for i in xrange(tree.fqnmrfit) if abs(tree.fqmrifit[i]/100000000)==2]
  mr_fits.sort(key=lambda x: x.nll)
  
  mr_fits = [f for f in mr_fits if f.status>=0]
  
  if not len(mr_fits): return False
  
  bf_nll = mr_fits[0].nll
  fits_with_ok_nll = [f for f in mr_fits]
  
  accepted_fit = None
  
  #~ if True:
    #~ try:
      #~ print neut_codes[tree.mode]
    #~ except KeyError:
      #~ print tree.mode
    #~ 
    #~ for fit in mr_fits:
      #~ print fit.AsString(),"\t",fit.GetNPions()
      #~ 
    #~ print "\n===\n"
    
  
  #~ for fit in fits_with_ok_nll:
    #~ if fit.GetNMuons()==1 and fit.GetNPions()==1 and fit.nring<4:
      #~ accepted_fit = fit
      #~ break
  for fit in fits_with_ok_nll:
    if fit.GetNMuons()==1 and fit.GetNElectrons()==2 and fit.nring<4:
      accepted_fit = fit
      break
  
  if not accepted_fit: return False
    
  for ring,pid in enumerate(accepted_fit.ring_pid):
    if pid==kMUON: break
    
  if not accepted_fit.ring_momentum[ring] > 200. : return False
  
  return True
  
def select_NCpiPlus(tree):
  ind = best_fq_mr_index(tree)
  nring = tree.fqmrnring[ind]
  if tree.fqmrpcflg[ind]<0.         : return False
  if not tree.fqmrnring[ind] == 1   : return False
  muon_inds     = [i for i in xrange(nring) if tree.fqmrpid[ind*6+i]==kMUON]
  pion_inds     = [i for i in xrange(nring) if tree.fqmrpid[ind*6+i]==kPION]
  electron_inds = [i for i in xrange(nring) if tree.fqmrpid[ind*6+i]==kELECTRON]
  
  if len(pion_inds)!=1:
      return False     # require pion like ring
      
  return True

def select_mr_nuecc(tree):
  for i in range(tree.nring):
    if tree.ip[i] == 2 and tree.amome[i] > 1000: return True
  return False

def get_pmax_ring(tree):
  if tree.nring == 1:
    return 0
  
  # return ith in NRING of the highest momentum in event
  imax = 0
  PMAX = -1
  
  for i in range(tree.nring):
    if   tree.ip[i] == 2: P = tree.amome[i] # e-like ring momentum
    elif tree.ip[i] == 3: P = tree.amomm[i] # mu-like ring momentum
    else: P = 0.
    
    if PMAX < P:
      PMAX = P
      imax = i
        
  return imax
    
def select_pi0(tree, nomass = False):
    # 2 rings
    if tree.nring != 2: return False
    #record_cut("2ring")
    # 1st and 2nd ring are e-like
    if  tree.ip[0] != 2 or tree.ip[1] != 2: return False
    #record_cut("erings")
    # no decay electrons
    if tree.nmue != 0: return False
    #record_cut("dcye")
    if nomass:
        return True
    mass, mom, costh =  get_pi0_mmd(tree)
    if not select_pi0_mass(mass): return False
    return True

# Calculate Pi0 Properties
def ix(a, b, c):
    return 3*6*a + 6*b + c

def get_pi0_mmd(tree):
    # Reference : 
    #  /home/atmpd/skrep/07d/atmpd/src/analysis/official_macro_sk3/official_pi0.kumac
    #  $ATMPD_ROOT/inc/apmsfit.h
    if tree.nring != 2:
        print " OfficialEventParser::get_pi0mom   Don't call me if you don't have 2 rings! "
        exit(-1)

    path = 0.
    tmpdir = [0, 0, 0]

    msdir = tree.msdir
    p0 = tree.amome[0]
    p1 = tree.amome[1]
    

    # msdir[nring][0:x, 1:y, 2:z][0:gamma, 1:e, 2:mu, 3:pi, 4:kaon, 5:proton]
    costh = msdir[ix(0,0,1)] * msdir[ix(1,0,1)] + msdir[ix(0,1,1)] * msdir[ix(1,1,1)] + msdir[ix(0,2,1)] * msdir[ix(1,2,1)]
    mass = sqrt(2. * p0 * p1 * (1.-costh))

    for i in range(3):
        tmpdir[i] = p0*msdir[ix(0,i,1)] + p1*msdir[ix(1,i,1)];
        path += tmpdir[i]**2

    pi0mom = sqrt(path)
    #pidir = tmpdir[2]/pi0mom

    return mass, pi0mom, costh


def get_inv_mass(tree, rings = []):
    # Assume all e-like rings are photons
    mumass = 105.7

    if tree.nring == 1:
        return 0
    
    msdir = tree.msdir
    Esum = 0
    psum = [ 0, 0, 0 ]

    if not rings: 
        rings = range(tree.nring)

    for r in rings:
        # muon
        if tree.ip[r] == 3: 
            ptot = tree.amomm[r]
            E = sqrt(ptot**2 + mumass**2)
            p = [ msdir[ix(r,xx,2)] * ptot for xx in range(3) ]
        # electron
        else:
            ptot = tree.amome[r]
            E = ptot
            p = [ msdir[ix(r,xx,1)] * ptot for xx in range(3) ]
        
        Esum += E
        psum = [ a+b for a, b in zip(psum, p) ]

    try:
        mass = sqrt(Esum**2 - sum([ pi**2 for pi in psum ]))
    except ValueError:
        mass = 0
    return mass


# it returns True if pi0mass is in 85 < mass < 185 .
def select_pi0_mass( pi0mass):
    if pi0mass <= 85. or pi0mass >= 185.: return False
    #record_cut("pi0mass")
    return True


# Flux/Osc Weight
def weight(tree, fileType, osc = 1, banff=False):
    if fileType == "data": return 1

    # Set sensible defaults
    weight = 1.0

    if fileType == "nue_x_numuflx": weight = 0.0 # avoid double counting

    if osc == 0:
        pass
    elif is_nc(tree):
        pass
    elif osc < 3:
        pos = { 'nue':0, 'nuebar':0,   # nue to nue
                'numu':4, 'numubar':4, # numu to numu
                'nue_x_numuflx':3 }    # numu to nue

        index = (osc - 1)*6 + pos[fileType ]
        weight = tree.t2kposc[index]
    elif osc == 3:
        if fileType != 'numu' and fileType != 'numubar':
            pass
        else:
            L = 295. #km
            sinth23 = 1.0
            dm23 = 2.4e-3
            posc = 1. - sinth23 * sin( 1.267 * dm23 * L / tree.pnu[0] )**2
            weight = posc
    
    if banff:  weight *= banff.fWeight  # BANFF central value
    else:      weight *= 1.0            # ND Scale

    return weight


############################################

def best_fq_mr_index(tree):
    fits = [(tree.fqmrnll[i],i,tree.fqmrifit[i]) for i in xrange(tree.fqnmrfit) if tree.fqmrifit[i]/100000000==2]
    #~ fits = [(tree.fqmrnll[i],i,tree.fqmrifit[i]) for i in xrange(tree.fqnmrfit)]
    nll,index,fit_id = min(fits)
    return index


def CosThetaNuLepton(tree) :
    ipmax = get_pmax_ring(tree)
    # it would be dir[0][0], dir[0][1] and dir[0][2] in the OLD tree c++ variable
    costnl = tree.dir[0+3*ipmax]*dir_nu[0] + tree.dir[1+3*ipmax]*dir_nu[1] + tree.dir[2+3*ipmax]*dir_nu[2]

    return costnl


def EnuReco(tree,isnumu) :
    V = binding_energy
    
    ind = best_fq_mr_index(tree)
    
    if isnumu==True : 
        leptonMass = muonMass
        muon_inds     = [i for i in xrange(6) if tree.fqmrpid[ind*6+i]==kMUON]
        if not muon_inds: return -999.
        leptonMom = tree.fqmrmom[ind*6+muon_inds[0]]
        #leptonMom = tree.amomm[0]
    else : 
        leptonMass = electronMass
        electron_inds = [i for i in xrange(6) if tree.fqmrpid[ind*6+i]==kELECTRON]
        if not electron_inds: return -999.
        leptonMom = tree.fqmrmom[ind*6+electron_inds[0]]

    leptonEn = math.sqrt(leptonMom**2 + leptonMass**2)
    protonMass2 = protonMass**2
    neutronMass2 = neutronMass**2
    costh = CosThetaNuLepton(tree)

    enu = ((neutronMass-V)*leptonEn - leptonMass*leptonMass/2.0 + neutronMass*V - V*V/2.0 + (protonMass2 - neutronMass2)/2.0 ) / (neutronMass - V - leptonEn + leptonMom*costh )
    
    #if enu<=0 :
    #    print "tree.amomm[ipmax]", tree.amomm[ipmax],"tree.amome[ipmax]", tree.amome[ipmax], "leptonMom", leptonMom, "leptonEn", leptonEn, "protonMass2", protonMass2, "neutronMass2", neutronMass2, "costh", costh, "enu", enu

    #print "tree.amomm[ipmax]", tree.amomm[ipmax],"tree.amome[ipmax]", tree.amome[ipmax], "leptonMom", leptonMom, "leptonEn", leptonEn, "protonMass2", protonMass2, "neutronMass2", neutronMass2, "costh", costh, "enu", enu
    
    return enu

# Calculate  Ereco using a recoil delta, for CC1Pi
def EnuRecoDelta(tree,isnumu) :
    V=binding_energy
    ipmax = get_pmax_ring(tree)
    ind = best_fq_mr_index(tree)
    if isnumu==True : 
        leptonMass = muonMass
        muon_inds     = [i for i in xrange(6) if tree.fqmrpid[ind*6+i]==kMUON]
        if not muon_inds: return -999
        leptonMom = tree.fqmrmom[ind*6+muon_inds[0]]
        #leptonMom = tree.amomm[0]
    else : 
        leptonMass = electronMass
        electron_inds = [i for i in xrange(6) if tree.fqmrpid[ind*6+i]==kELECTRON]
        if not electron_inds: return -999
        leptonMom = tree.fqmrmom[ind*6+electron_inds[0]]
    
    leptonMass2 = leptonMass**2
    leptonEn2 = (leptonMom**2 + leptonMass2);
    leptonEn = math.sqrt(leptonEn2)
    protonMass2 = protonMass**2
    neutronMass2 = neutronMass**2
    deltaMass2 = deltaMass**2
    costh = CosThetaNuLepton(tree)
    
    enu = ((neutronMass-V)*leptonEn - leptonMass2/2.0 + neutronMass*V - V*V/2.0 + (deltaMass2 - neutronMass2)/2.0 ) / (neutronMass - V - leptonEn + leptonMom*costh )
    
    return enu

# Calculate Ereco assuming CC1Pi and using pion ring
def EnuRecoCC1Pi(tree,isnumu) :
    return -1
    return enu


def ErecoCut(tree,isnumu) :
    ereco = EnuReco(tree,isnumu)
    if ereco<0.0 : return False #EnuReco(tree)<1250.0
    return True


def Q2(tree,isnumu) :
    V=binding_energy
    ipmax = get_pmax_ring(tree)
    
    if isnumu : 
        leptonMass = muonMass
        leptonMom = tree.amomm[ipmax]
        #leptonMom = tree.amomm[0]
    else : 
        leptonMass = electronMass
        leptonMom = tree.amome[ipmax]
        #leptonMom = tree.amome[0]

    leptonEn = math.sqrt(leptonMom**2 + leptonMass**2)
    protonMass2 = protonMass**2
    neutronMass2 = neutronMass**2
    costh = CosThetaNuLepton(tree)

    enu = ((neutronMass-V)*leptonEn - leptonMass*leptonMass/2.0 + neutronMass*V - V*V/2.0 + (protonMass2 - neutronMass2)/2.0 ) / (neutronMass - V - leptonEn + leptonMom*costh )
    #enu = EnuReco(tree,isnumu)
    
    p_nu    = TLorentzVector (enu*dir_nu[0], enu*dir_nu[1], enu*dir_nu[2], enu);
    dir_lep = [ tree.dir[0+3*ipmax], tree.dir[1+3*ipmax], tree.dir[2+3*ipmax] ];
    p_lep = TLorentzVector(dir_lep[0]*leptonMom, dir_lep[1]*leptonMom, dir_lep[2]*leptonMom, leptonEn);
    q2 = (p_nu - p_lep).Mag2();

    return -q2 * 1E-6


def FillVarTree(tree,isnumu) :
    ipmax = get_pmax_ring(tree)
    tmode[0] = tree.mode
    terec[0] = EnuReco(tree,isnumu)  # in MeV/c in the tree, then converted in the t2k analysis code
    tenu[0] = tree.pnu[0] # momentum at neutrino interaction
    twgt[0] = 1 
    if isnumu :
        #tplep[0] = tree.amomm[0]/1000.0; # given in MeV/c
        tplep[0] = tree.amomm[ipmax]/1000.0; # given in MeV/c
    else :
        tplep[0] = tree.amome[ipmax]/1000.0; # given in MeV/c
    tcostheta[0] = CosThetaNuLepton(tree)        
    tQ2[0] = Q2(tree,isnumu)
    #print tQ2[0] 
    tevis[0] = tree.evis/1000.0 # given in MeV/c
    tnring[0] = tree.nring 
    #print tnring[0], terec[0], tenu[0], tplep[0], tevis[0]
    
    terec_delta[0] = EnuRecoDelta(tree,isnumu)  # in MeV/c in the tree, then converted in the t2k analysis code
    terec_pi   [0] = EnuRecoCC1Pi(tree,isnumu)  # in MeV/c in the tree, then converted in the t2k analysis code


######################## Scaling factors ####################


nominalPOT = 1e21
#~ scales = {
    #~ "nue_x_numuflx":2.34490339e+00, 
    #~ "nuebar"       :4.66653529e-03, 
    #~ "numu"         :2.23465271e+00, 
    #~ "nue"          :4.46572986e-02, 
    #~ "numubar"      :7.66312101e-02
#~ }
#~ scales = {          # Nexp    N truefid  Nfiles     
    #~ "numu"         :1418.68 / 634998.0 * 1000.,
    #~ "numubar"      :48.6039 / 319640.0 * 504.,
    #~ "nue"          :28.312  / 633715.0 * 1000.,
    #~ "nuebar"       :2.98885 / 319689.0 * 504.,
    #~ "nue_x_numuflx":1490.15 / 320252.0 * 504.,
    #~ "ncel"         :1.0
#~ }

#Normalizing: 1418.68/634998
#Normalizing: 28.312/633715
#Normalizing: 48.6039/319640
#Normalizing: 1490.15/320252
#Normalizing: 2.98885/319689

#######################
sk_event_counts = { 
  "numu_x_numu"         : 1396.245 ,
  "nue_x_nue"           :   26.663 ,
  "numu_bar_x_numu_bar" :   50.964 ,
  "nue_bar_x_nue_bar"   :    2.884 ,
  "numu_x_nue"          : 1483.994 ,
  "numu_bar_x_nue_bar"  :   52.939 
}

tree_names = { "tnumu"        : "numu_x_numu"         ,
               "tnue"         : "nue_x_nue"           ,
               "tnumubar"     : "numu_bar_x_numu_bar" ,
               "tnuebar"      : "nue_bar_x_nue_bar"   ,
               "tnuesig"      : "numu_x_nue"          ,
               "tnuebarsig"   : "numu_bar_x_nue_bar"  
             }
               
####################### 





class T2KAnaError(Exception):
  pass

if __name__=="__main__":
  #import optparse
  #parser = optparse.OptionParser(usage="%prog [options]..",description="Make plots of size of systematic/statistical effects")
  #parser.add_option('-o', dest='osc'  ,default=files  ,help='The unoscillated input ROOT file')
  #options,other = parser.parse_args()
  #main(options.osc) 

  main(sys.argv[1:]) 


