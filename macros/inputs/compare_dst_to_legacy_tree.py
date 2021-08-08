#!/usr/bin/env python

#
# Convert the minituple, mtuple from the SK group into a VALOR DST format to be read by GSSFM::FillFromDST
#
# NOTE: Be aware that the sample, beam and reac mode for each event is stored as an int (in the output tree),
# so if any of these enumerations change, the MC file will have to be regenerated
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# April 2019

import ROOT
import glob
from array import array
import sys
import os
import collections
from collections import namedtuple

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))
from UtlPy import c_to_pyroot_array_hack

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()




#_______________________________________________________________________

def main(sdk_filename, t2k_filename, sample_str, beam_str, mode_str):

  sample = util.SampleFromString(sample_str)
  if sample == conv.ksUndefined:
    print "Did not recognise the input sample. Exiting."
    sys.exit()

  beam = util.BeamConfFromString(beam_str)
  if beam == conv.kbUndefined:
    print "Did not recognise the input beam. Exiting."
    sys.exit()

  mode = util.ReacModeFromString(mode_str)
  if mode == conv.kmUndefined:
    print "Did not recognise the input mode. Exiting."
    sys.exit()

  print "Read SDK and T2K events for the events with sample =", util.SampleAsString(sample), \
    ", beam =", util.BeamConfAsString(beam), ", mode = ", util.ModeAsString(mode)

  # Open the input files and get the trees
  sdk_file = ROOT.TFile(sdk_filename, "read")
  t2k_file = ROOT.TFile(t2k_filename, "read")
  tree_sdk = sdk_file.Get("selected_events_tree")
  tree_t2k = t2k_file.Get("selected_events_tree")

  # Set up lists and a named tuple to hold the info from each selected event
  events_sdk = []
  events_t2k = []
  event_info = namedtuple('event', ['sample', 'beam', 'reac_mode', 'enu_true', 'enu_reco', 'Theta_reco', 'wght', 'wght_pot', 'wght_tun'])

  # Loop over the SDK tree and extract the selected events
  print "Extracting SDK selected events..."
  n_duplicates_sdk = 0
  for event in tree_sdk:
    theta_reco = 0
    if util.IsELike(sample):
      theta_reco = event.kv_reco2

    this_event = event_info(
      sample     = event.sample,
      beam       = event.beam,
      reac_mode  = event.mode,
      enu_true   = event.kv_true1,
      enu_reco   = event.kv_reco1,
      Theta_reco = theta_reco,
      wght       = event.wght,
      wght_pot   = event.wght_pot,
      wght_tun   = event.wght_tun
    )

    if event.sample != sample:
      print "WARNING: SDK sample ["+ util.SampleAsString(event.sample)+ "] doesn't match the input sample. Something's gone wrong in the event selection."
    if event.beam != beam:
      print "WARNING: SDK beam ["+ util.BeamConfAsString(event.beam)+ "] doesn't match the input beam. Something's gone wrong in the event selection."
    if event.mode != mode:
      print "WARNING: SDK mode ["+ util.ModeAsString(event.mode)+ "] doesn't match the input mode. Something's gone wrong in the event selection."

    # Check for duplicate events
    for event_stored in events_sdk:
      if namedtuples_are_equal(event_stored, this_event):
        n_duplicates_sdk += 1
        break

    events_sdk.append(this_event)
  print "  ...finished extracting SDK selected events. Found", n_duplicates_sdk, "duplicate events"

  # Loop over the T2K legacy tree and extract the selected events
  print "Extracting T2K selected events..."
  n_duplicates_t2k = 0
  for event in tree_t2k:
    t2k_sample = t2k_sample_to_sdk_sample(event.sample_str)
    theta_reco = 0
    if util.IsELike(sample):
      theta_reco = event.kv_reco2
    this_event = event_info(
      sample     = t2k_sample_to_sdk_sample(event.sample_str),
      beam       = t2k_beam_to_sdk_beam(event.sample_str),
      reac_mode  = mode,            # TODO: Ideally this should be converted from the legacy T2K mode and checked against the input mode, but for now I'm just going to trust that I set the mode correctly when extracting these events to make the input trees
      enu_true   = event.enu_true,
      enu_reco   = event.enu_reco,
      Theta_reco = theta_reco,
      wght       = event.wght,
      wght_pot   = event.wght_pot,
      wght_tun   = event.wght_tun
    )
    if t2k_sample_to_sdk_sample(event.sample_str) != sample:
      print "WARNING: T2K sample ["+ util.SampleAsString(sample)+ "] doesn't match the input sample. Something's gone wrong in the event selection."
    if t2k_beam_to_sdk_beam(event.sample_str) != beam:
      print "WARNING: T2K beam ["+ util.BeamConfAsString(sample)+ "] doesn't match the input beam. Something's gone wrong in the event selection."

    # Check for duplicate events
    for event_stored in events_t2k:
      if namedtuples_are_equal(event_stored, this_event):
        n_duplicates_t2k += 1
        break

    events_t2k.append(this_event)
  print "  ...finished extracting T2K selected events. Found", n_duplicates_t2k, "duplicate events"

  print "Number of selected SDK events =", len(events_sdk)
  print "Number of selected T2K events =", len(events_t2k)
  if len(events_sdk) != len(events_t2k):
    print "WARNING: The number of SDK and T2K events are not equal!"
  n_tot = len(events_t2k)

  # Loop over the SDK event list and try to find matching events in the T2K event list, then eliminate any matches
  n_matches = 0
  cpy_events_sdk = list(events_sdk)
  for i_sdk_event in range(0, len(events_sdk)):
    sdk_event = events_sdk[i_sdk_event]

    # If all kinematic values match, then we consider this to be the same event
    sdk_event_matches_t2k_event = False
    cpy_events_t2k = list(events_t2k)
    for i_t2k_event in range(0, len(events_t2k)):
      t2k_event = events_t2k[i_t2k_event]
      if    ROOT.TMath.Abs(sdk_event.enu_true   - t2k_event.enu_true)   < conv.kASmallNum \
        and ROOT.TMath.Abs(sdk_event.enu_reco   - t2k_event.enu_reco)   < conv.kASmallNum \
        and ROOT.TMath.Abs(sdk_event.Theta_reco - t2k_event.Theta_reco) < conv.kASmallNum:

        # The SDK event kinematics matches the T2K event kinematics, so compare the weights
        mismatch = False
        event_info_str = "For event [Et ="+ str(sdk_event.enu_true)+ ", Er ="+ str(sdk_event.enu_reco)+ ", th ="+ str(sdk_event.Theta_reco)+ "]:"

        if ROOT.TMath.Abs(sdk_event.wght     - t2k_event.wght)     > conv.kASmallNum:
          mismatch = True
          event_info_str += "\n  wght mismatch: SDK = "+ str(sdk_event.wght)+ " T2K = "+ str(t2k_event.wght)
        if ROOT.TMath.Abs(sdk_event.wght_pot - t2k_event.wght_pot) > conv.kASmallNum:
          mismatch = True
          event_info_str += "\n  wght_pot mismatch: SDK = "+ str(sdk_event.wght_pot)+ " T2K = "+ str(t2k_event.wght_pot)
        if ROOT.TMath.Abs(sdk_event.wght_tun - t2k_event.wght_tun) > conv.kASmallNum:
          mismatch = True
          event_info_str += "\n  wght_tun mismatch: SDK = "+ str(sdk_event.wght_tun)+ " T2K = "+ str(t2k_event.wght_tun)

        if mismatch:
          print event_info_str, "\n"
        else:
          # The SDK and T2K events match perfectly, so remove them from their respective lists
          sdk_event_matches_t2k_event = True
          cpy_events_t2k.remove(t2k_event)
          n_matches += 1
          continue
    events_t2k = cpy_events_t2k
    if sdk_event_matches_t2k_event:
      # The SDK and T2K events match perfectly, so remove them from their respective lists
      cpy_events_sdk.remove(sdk_event)
  events_sdk = cpy_events_sdk

  print "Matching number of events =", n_matches, "out of", n_tot

  # Print the first 5 left-over events:
  if len(events_t2k) > 0:
    print "\n\n\nLeft over T2K events (prints up to 5 events):"
    for i_t2k_event in range(0, len(events_t2k)):
      if i_t2k_event > 4:
        break
      print events_t2k[i_t2k_event]

  if len(events_sdk) > 0:
    print "\n\n\nLeft over SDK events (prints up to 5 events):"
    for i_sdk_event in range(0, len(events_sdk)):
      if i_sdk_event > 4:
        break
      print events_sdk[i_sdk_event]

#_______________________________________________________________________
def t2k_sample_to_sdk_sample(sample_str):
  if sample_str == "1Rmu" or sample_str == "RHC1Rmu":
    return conv.ks1Rmu
  elif sample_str == "1Re" or sample_str == "RHC1Re":
    return conv.ks1Re
  elif sample_str == "MultiRe":
    return conv.ks1Re1dcy
  else:
    print "ERROR: could not a convert T2K legacy sample into an SDK sample. Exiting"
    sys.exit()

#_______________________________________________________________________
def t2k_beam_to_sdk_beam(sample_str):
  if sample_str == "1Rmu" or sample_str == "1Re" or sample_str == "MultiRe":
    return conv.kbJPARC_FHC
  elif sample_str == "RHC1Rmu" or sample_str == "RHC1Re":
    return conv.kbJPARC_RHC
  else:
    print "ERROR: could not a convert T2K legacy sample into an SDK beam. Exiting"
    sys.exit()

#_______________________________________________________________________
def namedtuples_are_equal(nt1, nt2):
  if ROOT.TMath.Abs(nt1.sample     - nt2.sample    ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.beam       - nt2.beam      ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.reac_mode  - nt2.reac_mode ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.enu_true   - nt2.enu_true  ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.enu_reco   - nt2.enu_reco  ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.Theta_reco - nt2.Theta_reco) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.wght       - nt2.wght      ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.wght_pot   - nt2.wght_pot  ) < conv.kASmallNum and \
     ROOT.TMath.Abs(nt1.wght_tun   - nt2.wght_tun  ) < conv.kASmallNum:
    return True
  return False

# def t2k_mode_to_sdk_mode(mode_str):
#   if mode_str == "numuccqe"       : return conv.kmNuMuCCQE
#   if mode_str == "numucc1pi"      : return conv.
#   if mode_str == "numucccoh"      : return conv.
#   if mode_str == "numumec"        : return conv.
#   if mode_str == "numuccoth"      : return conv.
#   if mode_str == "numunc1piC"     : return conv.
#   if mode_str == "numunc1pi0"     : return conv.
#   if mode_str == "numunccoh"      : return conv.
#   if mode_str == "numuncoth"      : return conv.
#   if mode_str == "oscnumuccqe"    : return conv.
#   if mode_str == "oscnumucc1pi"   : return conv.
#   if mode_str == "oscnumucccoh"   : return conv.
#   if mode_str == "oscnumumec"     : return conv.
#   if mode_str == "oscnumuccoth"   : return conv.
#   if mode_str == "numubarccqe"    : return conv.kmNuMuBarCCQE
#   if mode_str == "numubarcc1pi"   : return conv.
#   if mode_str == "numubarcccoh"   : return conv.
#   if mode_str == "numubarmec"     : return conv.
#   if mode_str == "numubarccoth"   : return conv.
#   if mode_str == "numubarnc1piC"  : return conv.
#   if mode_str == "numubarnc1pi0"  : return conv.
#   if mode_str == "numubarnccoh"   : return conv.
#   if mode_str == "numubarncoth"   : return conv.
#   if mode_str == "oscnumubarccqe" : return conv.
#   if mode_str == "oscnumubarcc1pi": return conv.
#   if mode_str == "oscnumubarcccoh": return conv.
#   if mode_str == "oscnumubarmec"  : return conv.
#   if mode_str == "oscnumubarccoth": return conv.
#   if mode_str == "nueccqe"        : return conv.kmNuECCQE
#   if mode_str == "nuecc1pi"       : return conv.
#   if mode_str == "nuecccoh"       : return conv.
#   if mode_str == "nuemec"         : return conv.
#   if mode_str == "nueccoth"       : return conv.
#   if mode_str == "nuenc1piC"      : return conv.
#   if mode_str == "nuenc1pi0"      : return conv.
#   if mode_str == "nuenccoh"       : return conv.
#   if mode_str == "nuencoth"       : return conv.
#   if mode_str == "oscnueccqe"     : return conv.
#   if mode_str == "oscnuecc1pi"    : return conv.
#   if mode_str == "oscnuecccoh"    : return conv.
#   if mode_str == "oscnuemec"      : return conv.
#   if mode_str == "oscnueccoth"    : return conv.
#   if mode_str == "nuebarccqe"     : return conv.kmNuEBarCCQE
#   if mode_str == "nuebarcc1pi"    : return conv.
#   if mode_str == "nuebarcccoh"    : return conv.
#   if mode_str ==  "nuebarmec"     : return conv.
#   if mode_str == "nuebarccoth"    : return conv.
#   if mode_str == "nuebarnc1piC"   : return conv.
#   if mode_str == "nuebarnc1pi0"   : return conv.
#   if mode_str == "nuebarnccoh"    : return conv.
#   if mode_str == "nuebarncoth"    : return conv.
#   if mode_str == "oscnuebarccqe"  : return conv.
#   if mode_str == "oscnuebarcc1pi" : return conv.
#   if mode_str == "oscnuebarcccoh" : return conv.
#   if mode_str ==  "oscnuebarmec"  : return conv.
#   if mode_str == "oscnuebarccoth" : return conv.
#   if mode_str == "numunc1gamma"   : return conv.
#   if mode_str == "numubarnc1gamma": return conv.
#   if mode_str == "nuenc1gamma"    : return conv.
#   if mode_str == "nuebarnc1gamma" : return conv.

#_______________________________________________________________________
# def get_legacy_mc_tree(sample, beam, mode):
#   # file name depends on sample, tree name depends on flavour
#   file_name = os.environ["VALOR_EXPT"] + "/data/"
#   if sample == conv.ks1Rmu and beam == conv.kbJPARC_FHC:
#     file_name += "official/2018/superk/mc/mtuples_v2/valor_format/t2ksk.14a.neut5.3.2.13a_tuned_v1r1.250ka.numuselec.root"
#   elif sample == conv.ks1Rmu and beam == conv.kbJPARC_RHC:
#     file_name += "official/2018/superk/mc/mtuples_v2/valor_format/t2ksk.14a.neut5.3.2.13a_tuned_v1r1.m250ka.numuselec.root"
#   elif sample == conv.ks1Re and beam == conv.kbJPARC_FHC:
#     file_name += "official/2018/superk/mc/mtuples_v2/valor_format/t2ksk.14a.neut5.3.2.13a_tuned_v1r1.250ka.nueselec.root"
#   elif sample == conv.ks1Re and beam == conv.kbJPARC_RHC:
#     file_name += "official/2018/superk/mc/mtuples_v2/valor_format/t2ksk.14a.neut5.3.2.13a_tuned_v1r1.m250ka.nueselec.root"
#   elif sample == conv.ks1Re1dcy and beam == conv.kbJPARC_FHC:
#     file_name += "official/2018/superk/mc/mtuples_v2/valor_format/t2ksk.14a.neut5.3.2.13a_tuned_v1r1.250ka.cc1piselec.root"
#   else:
#     print "Could not get legacy file name. Exiting."
#     sys.exit()
#
#   tree_name = ""
#   if util.IsNuMu(mode):
#     tree_name = "tnumu"
#   elif util.IsNuMuBar(mode):
#     tree_name = "tnumubar"
#   elif util.IsNuE(mode):
#     tree_name = "tnue"
#   elif util.IsNuEBar(mode):
#     tree_name = "tnuebar"
#   elif util.IsNuMuToNuE(mode):
#     tree_name = "tsignue"
#   elif util.IsNuMuBarToNuEBar(mode):
#     tree_name = "tsignuebar"
#   else:
#     print "Could not get legacy tree name. Exiting."
#     sys.exit()
#
#   file = ROOT.TFile(file_name, "read")
#   tree = file.Get(tree_name)
#   ROOT.SetOwnership(tree, False)  # This stops the object being deleted from memory when it goes out of scope
#   return tree

#_______________________________________________________________________
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Compare input SK MC trees between the SDK and legacy T2K for a selected NEUT reac mode")
  
  parser.add_option('-i', '--in-sdk' , dest='sdk_selec_filename' ,                      default=""     , help="Name of the selected events file from the SDK.")
  parser.add_option('-t', '--in-t2xk', dest='t2xk_selec_filename',                      default=""     , help="Name of the selected events file from legacy T2K.")
  parser.add_option('-m', '--mode '  , dest='mode'               ,                      default=""     , help="SDK mode.")
  parser.add_option('-s', '--sample' , dest='sample'             ,                      default=""     , help="SDK sample.")
  parser.add_option('-b', '--beam'   , dest='beam'               ,                      default=""     , help="SDK beam.")
  parser.add_option('-f', '--flavour', dest='flav'               ,                      default=""     , help="SDK flavour.")
  parser.add_option('-o', '--osc'    , dest='osc'                , action='store_true', default="False", help="Treat the reac mode as oscillated?")

  o,other = parser.parse_args()

  main(o.sdk_selec_filename, o.t2xk_selec_filename, o.sample, o.beam, o.mode)
  if other:
    print "Warning - ignored arguments: ", other
