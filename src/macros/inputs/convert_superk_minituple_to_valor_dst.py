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

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))
from UtlPy import c_to_pyroot_array_hack

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util  =  valor.util
conv  =  valor.conv
opt   =  valor.Opt.Instance()
env   =  valor.Env.Instance()

python_undefined = 2L*sys.maxint+1

# At this point, we don't have an MC, so can't load in the config files to read the beam and sample lists, so hard-code them here
beams = [conv.kbJPARC_FHC, conv.kbJPARC_RHC]
samples = collections.OrderedDict()
samples[conv.kbJPARC_FHC] = [conv.ks1Rmu, conv.ks1Re, conv.ksMultiRe]
samples[conv.kbJPARC_RHC] = [conv.ks1Rmu, conv.ks1Re]

types = ["nuebar_x_nuebar"  ,  # The naming conventions used by the input SK MC files
         "nue_x_nue"        ,
         "numubar_x_nuebar" ,
         "numubar_x_numubar",
         "numu_x_nue"       ,
         "numu_x_numu"      ,
        ]

true_int_flav = dict()    # Map of the input file name conventions to their corresponding true nu flavour (at interaction time)
true_int_flav["nuebar_x_nuebar"]   = conv.kNuEBar
true_int_flav["nue_x_nue"]         = conv.kNuE
true_int_flav["numubar_x_nuebar"]  = conv.kNuEBar
true_int_flav["numubar_x_numubar"] = conv.kNuMuBar
true_int_flav["numu_x_nue"]        = conv.kNuE
true_int_flav["numu_x_numu"]       = conv.kNuMu

is_osc = dict()         # Map of the input file name conventions to whether the file corresponds to oscillated neutrinos or not
is_osc["nuebar_x_nuebar"]   = False
is_osc["nue_x_nue"]         = False
is_osc["numubar_x_nuebar"]  = True
is_osc["numubar_x_numubar"] = False
is_osc["numu_x_nue"]        = True
is_osc["numu_x_numu"]       = False

orig_treename = "mtuple"
rhc_file_name_conv = ".rhc."   # All RHC file names contain this string, search for this before FHC
fhc_file_name_conv = ".fhc."    # All FHC file names contain this string
t2k_beam_dir = [0.669764, -0.742179, 0.024223]  # This is a direct copy of kBeamDir from VALOR-SDK/src/Conventions/Constants.h. The reason for this is that pyROOT often completely fails to read 1D C arrays.  # TODO: use a getter function instead

# NOTE: See https://www.t2k.org/t2ksk/code/ntuplevariables  for reading in any other variables from the input MC trees

# TODO: Make the sample/mode/beam enumerations more robust. Currently they're just stored in the trees as the integer
# TODO: of the enumeration, which can silently cause problems with filling MCTs when the enumerations change (e.g. by
# TODO: adding another sample or mode). So in the TFile with the tree, store a lookup table (+ additional info e.g. the
# TODO: generation date) that relates the integer enumeration to AsString(enumeration). Then in GSSFM::FillFromDST read
# TODO: the int from the tree, use the lookup table (map<int,string>?) to convert this to a string, and then call
# TODO: FromString(string) to get the current VALOR enumeration.

#_______________________________________________________________________

def main(fhc_dir, rhc_dir, output_filename, modelist_file):

  # Get the modes
  if not modelist_file:
    print "ERROR: could not open the input mode list file!"
    sys.exit()
  ml = valor.ModeList(util.CompleteConfigPath(modelist_file))
  modes = ml.Modes()

  files_fhc = glob.glob(fhc_dir+"/*fhc*root")
  files_rhc = glob.glob(rhc_dir+"/*rhc*root")
  print "Registered the following FHC files, will select the relevant files\n", [f.rsplit("/")[-1] for f in files_fhc], "\n\n"
  print "Registered the following RHC files, will select the relevant files\n", [f.rsplit("/")[-1] for f in files_rhc], "\n\n"
  files = files_fhc + files_rhc

  f = files[0].rsplit("/")[-1]  # Gets the filename of the first file (strips the path)

  appropriate_files = get_appropriate_files(files)
  print "Reading from the following input files:\n", [f.rsplit("/")[-1] for f in appropriate_files], "\n"
  if len(appropriate_files) != len(list(set(appropriate_files))):
    print "ERROR: Duplicate input files detected!"
    util.Die(-1)

  print "Creating the output file", output_filename, "\n"
  output_file = ROOT.TFile(output_filename, "RECREATE")

  # Set up the output tree
  tree_out = ROOT.TTree("valor_t2k_dst","valor_t2k_dst")
  sample     = array('i', [-1])   # For why an array is used see https://root.cern.ch/how/how-write-ttree-python
  detector   = array('i', [-1])
  beam       = array('i', [-1])
  reac_mode  = array('i', [-1])
  enu_true   = array('d', [-1])
  enu_reco   = array('d', [-1])
  Theta_reco = array('d', [-1])
  wght       = array('d', [-1])
  wght_pot   = array('d', [-1])
  tree_out.Branch("sample"       , sample        , "sample/I"    )
  tree_out.Branch("detector"     , detector      , "detector/I"  )
  tree_out.Branch("beam"         , beam          , "beam/I"      )
  tree_out.Branch("reac_mode"    , reac_mode     , "reac_mode/I" )
  tree_out.Branch("enu_true"     , enu_true      , "enu_true/D"  )
  tree_out.Branch("enu_reco"     , enu_reco      , "enu_reco/D"  )
  tree_out.Branch("Theta_reco"   , Theta_reco    , "Theta_reco/D")
  tree_out.Branch("wght"         , wght          , "wght/D"      )
  tree_out.Branch("wght_pot"     , wght_pot      , "wght_pot/D"  )

  # Loop over the input files and fill the output tree
  iev = 0
  nev_input = [0 for i in appropriate_files]
  entries_per_beam_sample_mode = [[[0 for i in range(0, conv.kmOther)] for j in range(0, conv.kMaxNSamples)] for k in range(0, conv.kNBeamsDefined)]
  files_found_selected_mode_in = dict()
  for ifile in range(0, len(appropriate_files)):
    file_in_name = appropriate_files[ifile]
    file_in = ROOT.TFile(file_in_name, "read")
    in_tree = file_in.Get(orig_treename)
    nev_input[ifile] += in_tree.GetEntries()

    sample[0]       = get_selection(file_in_name)
    detector[0]     = conv.kdSuperK
    beam[0]         = get_beam(file_in_name)
    true_nu_flavour = get_flavour(file_in_name)
    is_oscillated   = get_osc(file_in_name)
    this_scale      = get_scale(file_in)
    wght_pot[0]     = get_pot_weight(this_scale)

    for event in in_tree:
      neut_mode     = event.mode
      reac_mode[0]  = NeutModeToValorMode(neut_mode, true_nu_flavour, is_oscillated)
      enu_true[0]   = event.pnu[0]                                                  # Already in GeV/c
      Theta_reco[0] = calc_theta_degrees(event.fq1rdir, event.fqnse, sample[0])     # event.fq1rdir is a 3D list that contains info on the lepton direction and is dependant on the sample
      enu_reco[0]   = calc_enu_reco(event, sample[0], reac_mode[0], Theta_reco[0])  # Already in GeV
      wght[0]       = group_misc_weights(event, neut_mode)                          # Group all of the standard multiplicative NIWG weights into a single weight

      if enu_reco[0] < 0. \
         or (is_oscillated and not util.ModeIsCC(reac_mode[0])):
        # Skip oscillated NC events, because NC is always just the sum of events in numu, nue, nutau. So the oscillation
        # is already accounted for.
        nev_input[ifile] -= 1
        continue

      # Debug selected events
      # if sample[0] == conv.ks1Rmu and beam[0] == conv.kbJPARC_FHC:
      #   if reac_mode[0] == conv.kmNuENC1piC           :
      #     if file_in_name in files_found_selected_mode_in:
      #         files_found_selected_mode_in[file_in_name] += 1
      #     else:
      #         files_found_selected_mode_in[file_in_name] = 1
          # print "Inserting to MC tree --> MC event #" + str(iev) + \
          #       "\n sample       =", str(sample[0]      ) , "(" + util.SampleAsString(sample[0])  + ")" + \
          #       "\n beam         =", str(beam[0]        ) , "(" + util.BeamConfAsString(beam[0])  + ")" + \
          #       "\n imode        =", str(reac_mode[0]   ) , "(" + util.ModeAsString(reac_mode[0]) + ")" + \
          #       "\n wght         =", str(wght[0]        ) + \
          #       "\n is NC 1gamma?" , str(NeutIsNC1Gamma(neut_mode)) + \
          #       "\n is 2p2h?"      , str(NeutIs2p2h(neut_mode)) + \
          #       "\n  True:", str(enu_true[0]) + \
          #       "\n  Reco:", str(enu_reco[0]), ",", str(Theta_reco[0]), "\n\n"

      entries_per_beam_sample_mode[beam[0]][sample[0]][reac_mode[0]] += 1
      iev += 1
      tree_out.Fill()

  output_file.cd()
  nev_output_tot = tree_out.GetEntries()
  tree_out.Write("valor_t2k_dst")
  output_file.Close()
  print "Written the MC DST to the output file", output_filename

  # Print debug info
  print "The number of events in each mode are:"
  for beam in beams:
    print "\n\nbeam =", util.BeamConfAsString(beam)
    for sample in samples[beam]:
      print "\n  sample =", util.SampleAsString(sample)
      tot_entries_sample = 0
      for mode in modes:
        if mode == -1 or mode == python_undefined:
          # Temporary fix until SampleReacMode.h is removed (as it is a fixed length, with unused modes being set to undefined)
          continue
        nevents = entries_per_beam_sample_mode[beam][sample][mode]
        tot_entries_sample += nevents
        print "    " + util.ModeAsString(mode) + ": " + str(nevents)
      print "Total", util.BeamConfAsString(beam), util.SampleAsString(sample), \
        " number of events across all modes =", tot_entries_sample

  print "The number of events in each input file are (minus any events with Erec < 0):"
  nev_input_tot = 0
  for ifile in range(0, len(appropriate_files)):
    nev_input_tot += nev_input[ifile]
    print " ", appropriate_files[ifile], ":", nev_input[ifile]
  print "The total number of events across all input files is = ", nev_input_tot
  print "The total number of events in the output DST = ", nev_output_tot

  if len(files_found_selected_mode_in) > 0:
    print "The selected mode/sample/beam combination was found in events from the following files:\n", files_found_selected_mode_in

  print "Done!"


#_______________________________________________________________________
def get_appropriate_files(files):
  appropriate_files = []
  for f in files:
    for partial_name in types:
      if partial_name in f:
        appropriate_files.append(f)
        continue
  return appropriate_files

#_______________________________________________________________________
# Get the POT scale corresponding to the input file
def get_scale(file_in):
  scale_hist = file_in.Get("norm")
  if not scale_hist:
    print "ERROR: could not find the scaling hist named \'{}\'".format("norm"), "in the input file:", file_in.GetName()
    util.Die(-1)
  scale_factor = scale_hist.GetBinContent(1)
  print "Using scaling factor of \'{}\'".format(scale_factor), "for input file", file_in.GetName()
  return scale_factor

#_______________________________________________________________________
# Get the selection corresponding to the input file
def get_selection(file_in_name):
  selection = conv.ksUndefined
  if "cc1piselec" in file_in_name:
    selection = conv.ksMultiRe
  elif "nueselec" in file_in_name:
    selection = conv.ks1Re
  elif "numuselec" in file_in_name:
    selection = conv.ks1Rmu
  else:
    print "FATAL: Could not determine the selection for input file", file_in_name
    util.Die(-1)
  return selection    # In python the enumerations are just treated as ints, so there's no need to cast to an int

#_______________________________________________________________________
# Get the beam corresponding to the input file
def get_beam(file_in_name):
  beam = conv.kbUndefined
  if rhc_file_name_conv in file_in_name:
    beam = conv.kbJPARC_RHC
  elif fhc_file_name_conv in file_in_name:
    beam = conv.kbJPARC_FHC
  else:
    print "FATAL: Could not determine the horn current for input file ", file_in_name
    util.Die(-1)
  return beam

#_______________________________________________________________________
# Get the true nu flavour corresponding to the input file
def get_flavour(file_in_name):
  true_nu_flavour = conv.kNuFlavourUndefined
  for filename_partial in types:
    if filename_partial in file_in_name:
      true_nu_flavour = true_int_flav[filename_partial]
      break
  if true_nu_flavour == conv.kNuFlavourUndefined or true_nu_flavour == 2**64 -1:
    print "FATAL: Could not determine true nu flavour for file", file_in_name
    util.Die(-1)
  return true_nu_flavour

#_______________________________________________________________________
# Determine if the input file corresponds to oscillated or un-oscillated events
def get_osc(file_in_name):
  for partial_filename in types:
    if partial_filename in file_in_name:
      return is_osc[partial_filename]
  print "FATAL: Could not determine whether the neutrinos are oscillated for file", file_in_name
  util.Die(-1)

#_______________________________________________________________________
# This returns 1 / MC POT. When the tree is later read in GSSFM::FillFromDST, this will be multiplied by the data POT to get the true POT weight
def get_pot_weight(this_scale):
  return  this_scale / 1E21

#_______________________________________________________________________
# Takes the scalar product between the beam and lepton directions
def calc_theta_degrees(lep_dir_in, size, sample):
  # lep_dir contains info on the lepton direction and is dependant on the sample
  lep_dir = c_to_pyroot_array_hack(lep_dir_in, size, 7, 3)
  costheta = -999.
  if util.SampleIsELike(sample):
    costheta = t2k_beam_dir[0] * lep_dir[0][1][0] + t2k_beam_dir[1] * lep_dir[0][1][1] + t2k_beam_dir[2] * lep_dir[0][1][2]
  if util.SampleIsMuLike(sample):
    costheta = t2k_beam_dir[0] * lep_dir[0][2][0] + t2k_beam_dir[1] * lep_dir[0][2][1] + t2k_beam_dir[2] * lep_dir[0][2][2]

  if  costheta < -1.00001 or costheta > +1.00001:
    print "ERROR: costheta = " << costheta << " --> outside allowed range"
    util.Die(1)

  LepTheta        = ROOT.TMath.ACos(costheta)          # radians
  LepThetaDegrees = LepTheta / ROOT.TMath.Pi() * 180.  # degrees
  if LepTheta < 0.:
    print "ERROR: Theta lepton is negative: " << LepThetaDegrees << " degrees"
    util.Die(1)

  return LepThetaDegrees

#_______________________________________________________________________
# Calculate the reconstructed neutrino energy, taking into account the Coulomb shift for CC events.
# See https://www.t2k.org/asg/oagroup/meeting/2018/2018-02-22/coulombcorrections for details.
# The Erec calculation uses the QE approximation, see e.g. TN315v5 p31.
def calc_enu_reco(event, sample, reac_mode, Theta_reco):
    ereco_MeV = None

    Vnuc = 27.0                          # Nominal binding energy, MeV
    mn   = conv.kNeutronMassGeV  * 1000. # MeV
    mp   = conv.kProtonMassGeV   * 1000. # MeV
    me   = conv.kElectronMassGeV * 1000. # MeV
    mm   = conv.kMuonMassGeV     * 1000. # MeV
    md   = conv.kDeltaMassGeV    * 1000. # MeV

    mass = None
    mom  = None
    if util.SampleIsELike(sample):
        mass = me
        mom = event.fqmome
    else:
        mass = mm
        mom = event.fqmomm

    mom_shift_nu    = -4.3   # MeV
    mom_shift_nubar =  3.3   # MeV
    mom_shift = None
    if "bar" in util.FlavourAsString(util.ReacModeToFinalFlavor(reac_mode)).lower():
        mom_shift = mom_shift_nubar
    else:
        mom_shift = mom_shift_nu
    shifted_mom = mom + mom_shift

    elep = ROOT.TMath.Sqrt(mass*mass +shifted_mom*shifted_mom)
    theta_reco_rads = Theta_reco * (ROOT.TMath.Pi() / 180.0)   # Theta_reco is input in degrees, so convert to radians here

    if util.IsNC(reac_mode):
        # Don't apply the Coulomb shift to NC events
        ereco_MeV = event.erec
    else:
        if sample != conv.ksMultiRe:
            ereco_MeV  = (mn - Vnuc)*elep - mass*mass/2.
            ereco_MeV += mn*Vnuc - Vnuc*Vnuc/2.
            ereco_MeV += (mp*mp - mn*mn)/2.
            ereco_MeV /= (mn - Vnuc - elep + shifted_mom * ROOT.TMath.Cos(theta_reco_rads))
        else:
            ereco_MeV  = (mp*elep-mass*mass/2.)
            ereco_MeV += (md*md - mp*mp)/2.
            ereco_MeV /= (mp - elep + shifted_mom * ROOT.TMath.Cos(theta_reco_rads))

    # Convert from MeV to GeV and set any overflow events with Erec > 29.9GeV to Erec = 29.9GeV
    return ROOT.TMath.Min(ereco_MeV / 1000., 29.9)

#_______________________________________________________________________
# Group all of the standard multiplicative NIWG and flux weights etc. into a single weight
def group_misc_weights(event, neut_mode):
  weight_NC1Gamma = 1.
  weight_flux_tuning = event.wgtflx    # NOTE (2020-02-19): Before run 10, we would calculate & apply the flux tuning weights at runtime using tuning histograms; however, these weights are now in the mtuple

  # The NC1gamma production cross-section in NEUT is off by factor 2
  # NOTE (2020-02-19): This is still true as of NEUT 5.4.0, so will be used in the run 10 analysis
  if NeutIsNC1Gamma(neut_mode):
    weight_NC1Gamma = 2.

  # NOTE: If, in the future, you want to put in fake data weights, do so here. See L1597 of SuperkMC.cxx and util::UtlMakeFakeData (in VALOR-T2K-LBL) for examples

  return weight_NC1Gamma * weight_flux_tuning

#___________________________________________________________________________
def NeutIsCCQE(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 1

#___________________________________________________________________________
def NeutIsCC1piC(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 11 or neut_reaction_mode == 13

#___________________________________________________________________________
def NeutIsCC1pi0(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 12

#___________________________________________________________________________
def NeutIsCCcoh(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 16

#___________________________________________________________________________
def NeutIs2p2h(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 2

#___________________________________________________________________________
def NeutIsCCDIS(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 26

#___________________________________________________________________________
def NeutIsCCMpi(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 21
#___________________________________________________________________________
def NeutIsCCmisc(neut_reaction_mode):
  # CCmisc=CCdiffractive+CC1K+CC1eta+CC1gamma
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return (neut_reaction_mode == 15 or
          neut_reaction_mode == 17 or
          neut_reaction_mode == 22 or
          neut_reaction_mode == 23 )   

# this function is not used anymore, though can be put back in case of need
#___________________________________________________________________________
#def NeutIsCCoth(neut_reaction_mode):
  #this is a catch all for everything CC that is not QE, 1pi, or coh
  # needs to be updated if another subcategory is added
 # neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  #return ( neut_reaction_mode == 17 or
   #        neut_reaction_mode == 21 or
    #       neut_reaction_mode == 22 or
     #      neut_reaction_mode == 23 or
      #     neut_reaction_mode == 26 )

#___________________________________________________________________________
def NeutIsNC1pi0(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 31 or neut_reaction_mode == 32

#___________________________________________________________________________
def NeutIsNC1Gamma(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 38 or neut_reaction_mode == 39

#___________________________________________________________________________
def NeutIsNC1piC(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 33 or neut_reaction_mode == 34

#___________________________________________________________________________
def NeutIsNCcoh(neut_reaction_mode):
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return neut_reaction_mode == 36

#___________________________________________________________________________
def NeutIsNCoth(neut_reaction_mode):
  #this is a catch all for everything NC that is not 1pi0, 1piC or coh
  # needs to be updated if another subcategory is added
  neut_reaction_mode = ROOT.TMath.Abs(neut_reaction_mode)
  return ( neut_reaction_mode == 35 or
           neut_reaction_mode == 41 or
           neut_reaction_mode == 42 or
           neut_reaction_mode == 43 or
           neut_reaction_mode == 44 or
           neut_reaction_mode == 45 or
           neut_reaction_mode == 46 or
           neut_reaction_mode == 51 or
           neut_reaction_mode == 52 )

#___________________________________________________________________________
def NeutIsCC(neut_reaction_mode):
  return NeutIsCCQE(neut_reaction_mode) or NeutIsCC1piC(neut_reaction_mode) or NeutIsCC1pi0(neut_reaction_mode) or \
         NeutIs2p2h(neut_reaction_mode) or NeutIsCCcoh(neut_reaction_mode) or NeutIsCCoth(neut_reaction_mode)

#___________________________________________________________________________
def NeutIsNC(neut_reaction_mode):
  return NeutIsNC1piC(neut_reaction_mode) or NeutIsNC1pi0(neut_reaction_mode) or NeutIsNC1Gamma(neut_reaction_mode) or \
         NeutIsNCcoh(neut_reaction_mode) or NeutIsNCoth(neut_reaction_mode)

#___________________________________________________________________________
# Convert a Neut MC reaction mode into a VALOR SDK reaction mode
# See https://t2k.org/docs/technotes/351/v02 tables 5 and 6 (p105) for an explanation of what each neut reac mode number is (for the 2018 analysis)
def NeutModeToValorMode(neut_mode, flav,  is_oscillated):
  valor_mode = conv.kmUndefined
  abs_mode = ROOT.TMath.Abs(neut_mode)   # +ve neut modes are for neutrinos, -ve neut modes are the equivalent antineutrino modes

  if flav == conv.kNuFlavourUndefined:
    print "WARNING: Tried to convert a neut reac mode to a valor reac mode, but the input flavour was undefined!"
    return conv.kmUndefined
  elif flav == conv.kNuTau or flav == conv.kNuTauBar:
    print "WARNING: NuTau and NuTauBar not supported when converting from neut reac mode to valor reac mode"
    return conv.kmUndefined

  if NeutIsCCQE(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCCQE
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECCQE
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCCQE
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCCQE
  elif NeutIsCC1piC(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCC1piC
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECC1piC
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCC1piC
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCC1piC
  elif NeutIsCC1pi0(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCC1pi0
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECC1pi0
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCC1pi0
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCC1pi0
  elif NeutIs2p2h(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCCMEC
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECCMEC
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCCMEC
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCCMEC
  elif NeutIsCCcoh(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCCcoh
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECCcoh
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCCcoh
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCCcoh
  elif NeutIsCCDIS(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCCDIS
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECCDIS
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCCDIS
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCCDIS
  elif NeutIsCCMpi(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCCMpi
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECCMpi
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCCMpi
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCCMpi
  elif NeutIsCCmisc(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuCCmisc
    if flav == conv.kNuE:
      valor_mode = conv.kmNuECCmisc
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarCCmisc
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarCCmisc
  elif NeutIsNC1piC(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuNC1piC
    if flav == conv.kNuE:
      valor_mode = conv.kmNuENC1piC
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarNC1piC
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarNC1piC
  elif NeutIsNC1pi0(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuNC1pi0
    if flav == conv.kNuE:
      valor_mode = conv.kmNuENC1pi0
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarNC1pi0
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarNC1pi0
  elif NeutIsNCcoh(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuNCcoh
    if flav == conv.kNuE:
      valor_mode = conv.kmNuENCcoh
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarNCcoh
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarNCcoh
  elif NeutIsNC1Gamma(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuNC1Gamma
    if flav == conv.kNuE:
      valor_mode = conv.kmNuENC1Gamma
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarNC1Gamma
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarNC1Gamma
  elif NeutIsNCoth(abs_mode):
    if flav == conv.kNuMu:
      valor_mode = conv.kmNuMuNCoth
    if flav == conv.kNuE:
      valor_mode = conv.kmNuENCoth
    if flav == conv.kNuMuBar:
      valor_mode = conv.kmNuMuBarNCoth
    if flav == conv.kNuEBar:
      valor_mode = conv.kmNuEBarNCoth
  else:
    print "WARNING: Could not convert a neut reac mode with code " << neut_mode << ", to a valor reac mode"
    return conv.kmUndefined

  if is_oscillated and (flav == conv.kNuE or flav == conv.kNuEBar) and not NeutIsNC(abs_mode):
    # We don't keep track of NC oscillated modes, as e.g. the nue template for NC covers nue -> nue, nue -> numu and
    # nue -> nutau (they would all be detected as the same)
    valor_mode = util.NuEModeToOscNuEMode(valor_mode)

  return valor_mode

#_______________________________________________________________________
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Reformat and optionally reduce input SK MC trees")
  
  parser.add_option('-o', '--out-file' , dest='output_filename', default="", help="Name for output MC DST file.")
  parser.add_option('-f', '--fhc-dir'  , dest='fhc_dir'        , default="", help="Directory containing all FHC input files")
  parser.add_option('-r', '--rhc-dir'  , dest='rhc_dir'        , default="", help="Directory containing all RHC input files")
  parser.add_option('-m', '--mode-list', dest='modelist_file'  , default="", help="XML file containing a list of active modes")

  o,other = parser.parse_args()

  main(o.fhc_dir, o.rhc_dir, o.output_filename, o.modelist_file)

  if other:
    print "Warning - ignored arguments: ", other
