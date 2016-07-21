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
    
try:
    ROOT.t2kana
except:
    ROOT.gROOT.ProcessLine(".L $T2KOSCANA/root/root2k.C");
    pass

def main(argv):

    files = glob(argv[0])
    treename = argv[1]
    outfile_numuCCQE = argv[2]
    outfile_nueCCQE  = argv[3]
    outfile_numuCC   = argv[4]
    outfile_nueCC    = argv[5]
    outfile_2pi0ring = argv[6]
    useereco    = argv[7]
    path = argv[8]+"/"
    
    if useereco=="erec_delta" : print "EREC_DELTA"
    elif useereco=="erec" : print "EREC"
    elif useereco=="erec_pi" : print "EREC_PI"
    else : 
        print "WRONG EREC"
        sys.exit()

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
    global tnmue
    tnmue = n.zeros(1, dtype=n.int)

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
        tree.Branch('nmue', tnmue, 'nmue/I')

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

    nev = tree.GetEntries()
    for iev in range(0,nev) :
        tree.GetEntry(iev)
        if iev%10000==0 : print "Event", iev
        
        #if iev==10 : break
        #if iev==10000 : break
        #print iev
        #print "Mom = ", tree.amomm[0], " - Ereco = ", EnuReco(tree)
        #isnumu = True
        #EnuReco(tree,isnumu)
        #Q2(tree,isnumu)
        #if iev==10 : return 
        
        # apply cuts (all samples are mutually exclusive)
        if select_fv(tree)==False : continue # fiducial volume events
        if select_fc(tree)==False : continue # fully contained events
        
        if select_numuccqe(tree):
            isnumu = True
            FillVarTree(tree,isnumu,useereco)
            tree_numuccqe.Fill()
            h_numuCCQE.Fill(1)
            continue
            
        if select_nueccqe(tree):
            isnumu = False
            FillVarTree(tree,isnumu,useereco)
            tree_nueccqe.Fill()
            h_nueCCQE.Fill(1)
            continue
        
        if reject1ring(tree)==False : continue # reject 1ring events

        if select_pi0(tree,False)==True : # 2-ring pi0
        #if select_pi0_fitqun(tree,False)==True : # 2-ring pi0
            isnumu = False
            if ErecoCut(tree,isnumu,useereco)==False : continue 
            FillVarTree(tree,isnumu,useereco)
            tree_2rpi0.Fill()
            h_2rpi0.Fill(1)
            continue
        
        #if select_mr_numucc(tree)==True : # multi-ring numucc
        if select_mr_numucc_highmom(tree)==True : # multi-ring numucc (highest mom ring)
        #if select_2r_numucc_highmom_fitqun(tree)==True : # multi-ring numucc (highest mom ring)
            isnumu = True
            if ErecoCut(tree,isnumu,useereco)==False : continue 
            FillVarTree(tree,isnumu,useereco)
            tree_numucc.Fill()
            h_numuCC.Fill(1)
            continue

        #if select_mr_nuecc(tree)==True : # multi-ring nuecc
        #if select_mr_nuecc_highmom(tree)==True : # multi-ring nuecc (highest mom ring)
        if select_2r_nuecc_highmom_fitqun(tree)==True : # multi-ring nuecc (highest mom ring)
            isnumu = False
            if ErecoCut(tree,isnumu,useereco)==False : continue 
            FillVarTree(tree,isnumu,useereco)
            tree_nuecc.Fill()
            h_nueCC.Fill(1)
            continue
        
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


###################################### Selection ###############################

def select_fv(tree):
    if tree.evis <= 30.: return False
    if tree.wall <= 200.: return False
    #record_cut("fv")
    return True

def select_fc(tree):
    #if tree.mode != 0 and tree.evis <= 9.: return False -- obsolete
    #if tree.mode != 0: return False
    #record_cut("evis")
    if tree.nhitac >= 16: return False
    #record_cut("nhitac")
    if tree.evclass != 1: return False
    return True
    
def true_fv(tree):
    if tree.mode == 0: return False
    return (tree.wallv > 200)

def select_elike(tree):
    return (ring_type(tree) == 2)

def select_mulike(tree):
    return (ring_type(tree) == 3)
    
def select_numuccqe(tree):
    if not select_fc(tree)       : return False
    if not select_fv(tree)       : return False
    if not tree.nring == 1       : return False
    if not tree.ip[0] == 3       : return False
    if not tree.amomm[0] > 200.  : return False
    if not tree.nmue <= 1        : return False
    if not tree.evis > 30.       : return False
    
    return True
    
def select_nueccqe(tree):
    if not select_fc(tree)       : return False
    if not select_fv(tree)       : return False
    if not tree.evis > 100.      : return False
    if not tree.nring == 1       : return False
    if not tree.ip[0] == 2       : return False
    if not tree.nmue == 0        : return False
    
    if not EnuReco(tree,False) < 1250.: return False
    
    fitqun_electron_likelihood = tree.fq1rnll[1]
    fitqun_pizero_likelihood   = tree.fq1rnll[3]
    fitqun_pizero_mass         = tree.fqpi0mass[0]
    
    if not (fitqun_electron_likelihood-fitqun_pizero_likelihood) < (175 - 0.875*fitqun_pizero_mass): return False
    
    return True
    

def select_mr_numucc(tree):
    for i in range(tree.nring):
        if tree.ip[i] == 3 and tree.amomm[i] > 200: return True
    for inmue in range(tree.nmue):
        evect = [ (tree.pos[xx] - tree.epos[3*inmue+xx])**2 for xx in range(3) ]
        edist = sqrt(sum(evect))
        if edist > 100: return True
    return False

def select_mr_nuecc(tree):
    for i in range(tree.nring):
        if tree.ip[i] == 2 and tree.amome[i] > 1000: return True
        return False

def prmslg_ip(tree, rnum = -1):
    if rnum == -1: rnum = get_pmax_ring(tree)
    eprob  = tree.prmslg[1+6*rnum]
    muprob = tree.prmslg[2+6*rnum]
    #print "prmslg: muprob=",muprob, "eprob=",eprob
    if eprob < muprob: return 2
    if muprob < eprob: return 3
    print "Cannot id with prmslg"
    return 0

def probms_ip(tree, rnum = -1):
    if rnum == -1: rnum = get_pmax_ring(tree)
    eprob  = tree.probms[1+6*rnum]
    muprob = tree.probms[2+6*rnum]
    #print "probms: muprob=",muprob, "eprob=",eprob
    if eprob > muprob: return 2
    if muprob > eprob: return 3
    print "Cannot id with probms"
    return 0

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


def get_pid(tree, rnum = -1):
    if tree.nring == 1:
        eprob  = tree.probms[1]
        muprob = tree.probms[2]
        return sqrt(-eprob) - sqrt(-muprob)
    else:
        if rnum == -1: rnum = get_pmax_ring(tree)
        eprob  = tree.prmslg[1+6*rnum]
        muprob = tree.prmslg[2+6*rnum]
        return (eprob - muprob)/100.



def ring_type(tree):
    maxp = 0
    minp = 0
    evtype = 0
    if tree.nring == 1:
        return tree.ip[0]

    for r in range(tree.nring):
        ip = tree.ip[r]
        #ip = prmslg_ip(tree, r)
        #ip = probms_ip(tree, r)
        #print "ip=", tree.ip[r], "prmslg_ip=", prm_ip, "probms_ip=", prob_ip
        
        if ip == 2:
            p = tree.amome[r]
            #minp = 100
        elif ip == 3:
            p = tree.amomm[r]
            #minp = 200
        else:
            p = tree.amom[r]
        if p > maxp and p > minp:
            maxp = p
            evtype = ip
        #if tree.nring > 1:
        #    print "  ==========="
        #    print "       r =",r
        #    print "      ip =",ip
        #    print "       p =",p
        #    print "    maxp = ", maxp
        #    print "  evtype =", evtype
    #if tree.nring > 1:
    #    print "++++++++++++++"
    #    print "evtype = ",evtype
    #    print "++++++++++++++"
    return evtype
    
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






def select_pi0_fitqun(tree, nomass = False):
    # 2 rings
    if tree.nring != 2: return False
    # fitqun
    #
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


# Truth
def is_nutau(tree):
    if abs( tree.ipnu[0]) == 16: return True
def is_numu(tree):
    if abs( tree.ipnu[0]) == 14: return True
def is_nue(tree):
    if abs( tree.ipnu[0]) == 12: return True


def original_pi0(tree):
    if tree.mode == 0: return False

    #Single pi0, coherent pi0 from neut
    if tree.mode == 31 or tree.mode == 32 or tree.mode == 36: return True

    # DIS with primary pi0
    if tree.mode == 46: 
        for i in range(2, tree.npar):
            if  tree.ipv[i] == 7 : 
                return True

    return False



def has_pi0(tree):
    if tree.mode == 0: return False

    for i in range(2, tree.npar):
        if  tree.ipv[i] == 7 : 
            return True
    for j in range(0, tree.npar2):
        if  tree.ipv2[j] == 7 : 
            return True
    return False

# Is the secondary visiblet
#invisible = [0, 4, 10, 13, 18, 20]
threshold = { 1:0.573, 2:0.573, 3:0.573,
              5:118.4, 6:118.4, 7:0.573,
              8:156.5, 9:156.5, 
              11:553.3, 12:553.3,
              14:1039.0, 
              19:877.7, 21:877.7}

def is_visible(ip, mom):
    # Neutrinos, neutrons, K0, geantinos are invisible
    if ip not in threshold: return False
    if mom < threshold[ip]: return False
    return True


# Simple pi0 event where we 
# only have 2 gammas in the final state
def simple_pi0(tree):
    _has_pi0 = False
    for i in range(2, tree.npar):
        if  tree.ipv[i] == 7 : 
            _has_pi0 = True
    
    if not _has_pi0:
        return False
    if tree.npar2 != 2:
        return False
    for j in range(0, tree.npar2):
        if tree.ipv2[j] != 1:
            return False


# Nu interaction produced a pi0 that was then lost
# through FSI
def lost_fsi_pi0(tree):
    # Nu interaction produces a pi0
    if not (tree.mode == 31 or tree.mode == 32 or tree.mode == 36):
        return False

    _has_pi0 = False
    for i in range(2, tree.npar):
        if  tree.ipv[i] == 7 : 
            return False

    return True
    

# Only secondary above threshold is a pi0 and gammas
def only_pi0(tree):
    _has_pi0 = False
    _only_pi0 = True

    for i in range(2, tree.npar):
        if  tree.ipv[i] == 7 : 
            _has_pi0 = True
    for j in range(0, tree.npar2):
        if  tree.ipv2[j] == 7 : 
            _has_pi0 = True
        elif tree.ipv2[j] == 1:
            pass # Don't check photons, expect those from pi0
        elif is_visible(tree.ipv2[j], tree.pmomv2[j]):
            _only_pi0 = False

    return (_has_pi0 and _only_pi0)

        

def is_nc(tree):
    if tree.mode == 0: return False
    return ( abs(tree.mode) > 30 )

def is_ncel(tree):
    if tree.mode == 0: return False
    return ( abs(tree.mode) == 51 or abs(tree.mode) == 52 )

def is_cc(tree):
    if tree.mode == 0: return False
    return not is_nc(tree)

def is_ccqe(tree):
    if tree.mode == 0: return False
    return ( abs(tree.mode) == 1 )


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


def CosThetaNuLepton(tree) :
    ipmax = get_pmax_ring(tree)
    # it would be dir[0][0], dir[0][1] and dir[0][2] in the OLD tree c++ variable
    costnl = tree.dir[0+3*ipmax]*dir_nu[0] + tree.dir[1+3*ipmax]*dir_nu[1] + tree.dir[2+3*ipmax]*dir_nu[2]

    return costnl


def EnuReco(tree,isnumu) :
    V=27.0
    ipmax = get_pmax_ring(tree)
    
    if isnumu==True : 
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
    
    #if enu<=0 :
    #    print "tree.amomm[ipmax]", tree.amomm[ipmax],"tree.amome[ipmax]", tree.amome[ipmax], "leptonMom", leptonMom, "leptonEn", leptonEn, "protonMass2", protonMass2, "neutronMass2", neutronMass2, "costh", costh, "enu", enu

    #print "tree.amomm[ipmax]", tree.amomm[ipmax],"tree.amome[ipmax]", tree.amome[ipmax], "leptonMom", leptonMom, "leptonEn", leptonEn, "protonMass2", protonMass2, "neutronMass2", neutronMass2, "costh", costh, "enu", enu
    
    return enu

# Calculate  Ereco using a recoil delta, for CC1Pi
def EnuRecoDelta(tree,isnumu) :
    V=27.0
    ipmax = get_pmax_ring(tree)
    
    if isnumu==True : 
        leptonMass = muonMass
        leptonMom = tree.amomm[ipmax]
        #leptonMom = tree.amomm[0]
    else : 
        leptonMass = electronMass
        leptonMom = tree.amome[ipmax]
        #leptonMom = tree.amome[0]
    
    leptonMass2 = leptonMass**2
    leptonEn2 = (leptonMom**2 + leptonMass2);
    leptonEn = math.sqrt(leptonEn2)
    protonMass2 = protonMass**2
    neutronMass2 = neutronMass**2
    deltaMass2 = deltaMass**2
    costh = CosThetaNuLepton(tree)
    
    enu = ((neutronMass-V)*leptonEn - leptonMass2/2.0 + neutronMass*V - V*V/2.0 + (deltaMass2 - neutronMass2)/2.0 ) / (neutronMass - V - leptonEn + leptonMom*costh )
    
    # Equation from Morgan Wascko http://wng.ift.uni.wroc.pl/karp45/presentations/Wascko-CC1pi.pdf
    #~ enu =  ( 0.5 * ( deltaMass2 - protonMass2 + 2*protonMass*leptonEn - leptonMass2) / 
             #~ (protonMass - leptonEn + math.sqrt(leptonEn2 - leptonMass2) * costh ) )
    
    return enu

# Calculate Ereco assuming CC1Pi and using pion ring
def EnuRecoCC1Pi(tree,isnumu) :
    if tree.nring < 2: return EnuReco(tree,isnumu)
    V=27.0
    
    rings_ordered_in_descending_energy = sorted( range(tree.nring), key=lambda i: tree.amomm[i], reverse=True )
    
    lepton_ringid = rings_ordered_in_descending_energy[0]
    pion_ringid   = rings_ordered_in_descending_energy[1]
    
    #~ print [(rings_ordered_in_descending_energy[i],tree.amomm[rings_ordered_in_descending_energy[i]]) for i in xrange(tree.nring)]
    #~ if (lepton_ringid != get_pmax_ring(tree)):
        #~ print get_pmax_ring(tree)
    
    leptonMom = 0
    if isnumu==True : 
        leptonMass = muonMass
        leptonMom = tree.amomm[lepton_ringid]
        #leptonMom = tree.amomm[0]
    else : 
        leptonMass = electronMass
        leptonMom = tree.amome[lepton_ringid]
        #leptonMom = tree.amome[0]

    leptonEn = math.sqrt(leptonMom**2 + leptonMass**2)
    protonMass2  = protonMass **2
    neutronMass2 = neutronMass**2
    leptonMass2  = leptonMass **2
    pionMass2    = pionMass   **2
    
    costhLepton = CosThetaNuLepton(tree)
    
    pionMom = tree.amomp[pion_ringid];
    pionEn  = math.sqrt(pionMom**2 + pionMass**2)
    
    costhPion = sum(( tree.dir[i+3*pion_ringid]*dir_nu[i] for i in xrange(3) ))

    pMuDotpPi = leptonMom * pionMom * sum(( tree.dir[i+3*pion_ringid]*tree.dir[i+3*lepton_ringid] for i in xrange(3) ))
    
    # Equation from Morgan Wascko http://wng.ift.uni.wroc.pl/karp45/presentations/Wascko-CC1pi.pdf
    enu = ( (leptonMass2 + pionMass2 - 2.*neutronMass*(leptonEn + pionEn) + 2.*pMuDotpPi) /
             (2.*(leptonEn + pionEn - leptonMom*costhLepton - pionMom*costhPion - neutronMass) ) )
    
    return enu


def ErecoCut(tree,isnumu,useereco) :
    
    #ereco = EnuReco(tree,isnumu)

    if useereco=="erec" :
        ereco = EnuReco(tree,isnumu)
    elif useereco=="erec_delta" :
        ereco = EnuRecoDelta(tree,isnumu) 
    elif useereco=="erec_pi" :
        ereco = EnuRecoCC1Pi(tree,isnumu) 
    else : 
        sys.exit()
    
    if ereco<0.0 : return False #EnuReco(tree)<1250.0
    return True


def Q2(tree,isnumu) :
    V=27.0
    ipmax = get_pmax_ring(tree)
    
    if is_numu(tree)==True : 
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


def FillVarTree(tree,isnumu,useereco) :
    ipmax = get_pmax_ring(tree)
    tmode[0] = tree.mode

    if useereco=="erec" :
        terec[0] = EnuReco(tree,isnumu)
    elif useereco=="erec_delta" :
        terec[0] = EnuRecoDelta(tree,isnumu) 
    elif useereco=="erec_pi" :
        terec[0] = EnuRecoCC1Pi(tree,isnumu) 
    else : 
        sys.exit()
    
    #terec[0] = EnuReco(tree,isnumu)  # in MeV/c in the tree, then converted in the t2k analysis code
    #terec[0] = (EnuReco(tree,isnumu) - tree.pnu[0])/tree.pnu[0] # for plotting
 

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
    tnmue[0] = tree.nmue
    #print tnring[0], terec[0], tenu[0], tplep[0], tevis[0]
    
    terec_delta[0] = EnuRecoDelta(tree,isnumu)  # in MeV/c in the tree, then converted in the t2k analysis code
    terec_pi   [0] = EnuRecoCC1Pi(tree,isnumu)  # in MeV/c in the tree, then converted in the t2k analysis code
    
    
def select_mr_numucc_highmom(tree): # highest momentum ring
    rtype = ring_type(tree) 
    ipmax = get_pmax_ring(tree)

    #if tree.nmue > tree.nring : return False # new cut --> # mu->e <= # rings  
    #if tree.nmue > 1 : return False # new cut --> # mu->e <= # rings  

    if rtype == 3 and tree.amomm[ipmax] > 200 : return True
    for inmue in range(tree.nmue):
        evect = [ (tree.pos[xx] - tree.epos[3*inmue+xx])**2 for xx in range(3) ]
        edist = sqrt(sum(evect))
        if edist > 100: return True
    return False





def select_2r_numucc_highmom_fitqun(tree): # highest momentum ring
    rtype = ring_type(tree) 
    ipmax = get_pmax_ring(tree)

    #if tree.nmue > tree.nring : return False # new cut --> # mu->e <= # rings  
    #if tree.nmue > 1 : return False # new cut --> # mu->e <= # rings  

    # fiTQun cut for highest mom ring (first seeded ring is the more energetic)
    #if tree.nring==2 : 
    if tree.fq2rpcflg[0] < 0 : return False # whether fiTQun believes the particle is exiting hte ID --> #MINUIT convergence check
    
    indMom = 0
    if tree.fq2rmom2[0]>tree.fq2rmom1[0] : indMom = 1
    if not tree.fq2rmom1[indMom]>200. : return False

    # selected_fitQun_16_9_14
    #fitqun_electron_likelihood = tree.fq2rnll[0] # first seeded ring
    #fitqun_pizero_likelihood = tree.fq2rnll[1] # first seeded ring
    #fitqun_pizero_mass         = tree.fqpi0mass[0]

    # selected_fitQun_16_9_14 II
    fitqun_electron_likelihood = min(tree.fq2rnll[0],tree.fq2rnll[2]) # first seeded ring
    fitqun_pizero_likelihood = min(tree.fq2rnll[1],tree.fq2rnll[3]) # first seeded ring
    fitqun_pizero_mass         = tree.fqpi0mass[0]

    if not (fitqun_pizero_likelihood-fitqun_electron_likelihood) < (175 - 0.875*fitqun_pizero_mass): return False

    if rtype == 3 and tree.amomm[ipmax] > 200 : return True
    for inmue in range(tree.nmue):
        evect = [ (tree.pos[xx] - tree.epos[3*inmue+xx])**2 for xx in range(3) ]
        edist = sqrt(sum(evect))
        if edist > 100: return True
    return False











def select_mr_nuecc_highmom(tree) : # highest momentum ring
    #if tree.nmue>0 : return False # new cut
    rtype = ring_type(tree) 
    ipmax = get_pmax_ring(tree)
    if rtype == 2 and tree.amome[ipmax] > 200: return True
    return False

def select_2r_nuecc_highmom_fitqun(tree) : # highest momentum ring w/ fiTQun
    #if tree.nmue>0 : return False # new cut --> no mu->e events 
    if tree.nring>2 : return False

    # fiTQun cut for highest mom ring (first seeded ring is the more energetic)
    #if tree.nring==2 : 
    if tree.fq2rpcflg[0] < 0 : return False # whether fiTQun believes the particle is exiting hte ID --> #MINUIT convergence check
    
    indMom = 0
    if tree.fq2rmom2[0]>tree.fq2rmom1[0] : indMom = 1
    if not tree.fq2rmom1[indMom]>200. : return False

    #print "indMom", indMom

    # selected_fitQun_16_9_14
    fitqun_electron_likelihood = tree.fq2rnll[indMom] # first seeded ring
    fitqun_pizero_likelihood = tree.fq2rnll[indMom+1] # first seeded ring
    fitqun_pizero_mass         = tree.fqpi0mass[0]

    # selected_fitQun_16_9_14 II
    #fitqun_electron_likelihood = min(tree.fq2rnll[0],tree.fq2rnll[2]) # first seeded ring
    #fitqun_pizero_likelihood = min(tree.fq2rnll[1],tree.fq2rnll[3]) # first seeded ring
    #fitqun_pizero_mass         = tree.fqpi0mass[0]

    if not (fitqun_pizero_likelihood-fitqun_electron_likelihood) < (175 - 0.875*fitqun_pizero_mass): return False

    rtype = ring_type(tree) 
    ipmax = get_pmax_ring(tree)
    if rtype == 2 and tree.amome[ipmax] > 200: return True
    return False





def reject1ring(tree) :
    nr = tree.nring 
    if nr<2 : return False
    return True


######################## Scaling factors ####################


nominalPOT = 1e21
#scales = {
#    "nue_x_numuflx":2.34490339e+00, 
#    "nuebar"       :4.66653529e-03, 
#    "numu"         :2.23465271e+00, 
#    "nue"          :4.46572986e-02, 
#    "numubar"      :7.66312101e-02
#}
scales = {          # Nexp    N truefid  Nfiles     
    "numu"         :1418.68 / 634998.0 * 1000.,
    "numubar"      :48.6039 / 319640.0 * 504.,
    "nue"          :28.312  / 633715.0 * 1000.,
    "nuebar"       :2.98885 / 319689.0 * 504.,
    "nue_x_numuflx":1490.15 / 320252.0 * 504.,
    "ncel"         :1.0
}

#Normalizing: 1418.68/634998
#Normalizing: 28.312/633715
#Normalizing: 48.6039/319640
#Normalizing: 1490.15/320252
#Normalizing: 2.98885/319689

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
     


