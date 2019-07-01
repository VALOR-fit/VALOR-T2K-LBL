#!/usr/bin/python

import ROOT
import os
import sys
from collections import namedtuple
import re

fmwk_var = "VALOR_FMWK"
sys.path.append(os.path.abspath(os.environ[fmwk_var] + "/src/valor/Utilities/"))
from UtlPy import c_to_pyroot_array_hack

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()
arridx= util.ArrIdx.Instance()


_first_xsec_index       = -999
_first_flux_index       = -999
_first_SKDetFSI_index   = -999

_name_flux_counter     = {}
_name_flux_counter_rhc = {}
_name_flux_counter["numu"]        = 0
_name_flux_counter["numubar"]     = 0
_name_flux_counter["nue"]         = 0
_name_flux_counter["nuebar"]      = 0
_name_flux_counter_rhc["numu"]    = 0
_name_flux_counter_rhc["numubar"] = 0
_name_flux_counter_rhc["nue"]     = 0
_name_flux_counter_rhc["nuebar"]  = 0

_possible_reac_type_names = [
    "CC QE",
    "CC 1pi",
    "CC coh",
    "CC 2p2h",
    "CC oth",
    "NC 1piC",
    "NC 1pi0",
    "NC coh",
    "NC oth",
    "NC 1gamma"
]


# Pyroot can't handle reading some 1D C++ arrays, so hardcode the array instead
# These are the PhysParamId_t integer values of each MultiRing SKDetFSI param
MultiRingSKDetFSIParams = [
    345,
    346,
    348,
    349,
    351,
    352,
    354,
    355,
]

#______________________________________________________________________________________________________________________________


def main(filename, t2k_legacy_config, overwrite):
    global _first_xsec_index
    global _first_flux_index
    global _first_SKDetFSI_index

    # Check if the output file already exists
    exists = os.path.isfile(filename)
    if not exists:
        print "Creating new file", filename
    else:
        if overwrite and filename != "":
            print "Overwriting file", filename
            os.remove(filename)
        else:
            print "A file named", filename, "already exists and the --overwrite option was not specified. Exiting."
            util.Die(-1)

    skexp = valor.SuperkExpectation()
    skexp.Read(os.environ[env_var] + "/config/" + t2k_legacy_config)

    outfile = open(filename, "w")
    outfile.write("# List of (\'SDK\', \'LEGACY\'), parameters")

    # Get the T2K param name and convert to an SDK param name
    params = valor.util.NamedParameterList("nuosc3f+nuisance")
    reac_types_all_okay = True
    N = params.GetNumberParameters()
    for ip in range(0, N):
        p = params.WhichParam(ip)
        # Skip non-systematic params and SCCA and SCCV (which aren't used in the analysis)
        if not params.IsNuisance(p) or p == conv.kNuisBANFFSCCA or p == conv.kNuisBANFFSCCV \
                or p == conv.kRHCNuisBANFFSCCA or p == conv.kRHCNuisBANFFSCCV:
            continue

        # Work out indexes of the first xsex, flux, and SKDetFSI params
        if util.IsBANFF(p):
            if util.IsFlux(p):
                if _first_flux_index < 0:
                    _first_flux_index = ip
            elif _first_xsec_index < 0:
                _first_xsec_index = ip
        elif _first_SKDetFSI_index < 0:
            _first_SKDetFSI_index = ip

        name_t2k = util.AsStringCompact(p)
        name_sdk = get_name(ip, p)
        outfile.write("    ( \'" + name_sdk + "\', \'" + name_t2k +"\' ),\n" )

    outfile.close()
    print "Successfully wrote the parameter file", filename
    print "Done!"
    return


def get_name(ip, p):
    global _name_flux_counter
    global _name_flux_counter_rhc
    name = ""
    if util.IsBANFF(p):
        if util.IsFlux(p):
            # Use a more descriptive name for flux params than e.g. banff00.
            flavour_str = get_flavour(ip, p, 0, 0, 0)
            counter = 0
            if util.IsRHCSystematic(p):
                counter = _name_flux_counter_rhc[flavour_str]
                _name_flux_counter_rhc[flavour_str] += 1
            else:
                counter = _name_flux_counter[flavour_str]
                _name_flux_counter[flavour_str] += 1
            replacement_str = "_" + flavour_str + "{0:0>2}".format(str(counter))
            name = util.AsStringCompact(p)
            name = name.replace("banff", "flux")
            name = re.sub(r'\d\d', replacement_str, name)
        else:
            name = util.AsStringNIWG(p)
    else:
        # Is SKDetFSI
        name = util.AsStringCompact(p)
    # Many param names start with 'f_', remove this
    name = name.replace("f_", "")
    if not name:
        print "ERROR: name not found for param", p, ". Exiting."
        util.Die(-1)
    return name



def get_flavour(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring):
    apply_to_modes = None
    if util.IsBANFF(p):
        if util.IsFlux(p):
            nu_PDG = None
            if util.IsRHCSystematic(p):
                index = ip - _first_flux_index - conv.kNumOfFluxParamsInBANFFfit
                nu_PDG = conv.kNuPdgCodeOfRHCFluxParamInBANFFfit[index]
            else:
                index = ip - _first_flux_index
                nu_PDG = conv.kNuPdgCodeOfFluxParamInBANFFfit[index]
            if nu_PDG == conv.kPdgNuMu:
                return "numu"
            elif nu_PDG == conv.kPdgNuMuBar:
                return "numubar"
            elif nu_PDG == conv.kPdgNuE:
                return "nue"
            elif nu_PDG == conv.kPdgNuEBar:
                return "nuebar"
            else:
                print "ERROR: unrecognised nu flavour PDG code", nu_PDG, "for ip =", ip, ", param name =", util.ParamAsString(p), ". Exiting."
                util.Die(-1)
        else:
            # Use the 3D array kXsApplyXSecSystToSKMode to decide whether each xsec param is applied to each sample
            apply_to_modes = get_apply_to_modes_xsec(ip, p, apply_table_xsec)
    else:
        apply_to_modes = get_apply_to_modes_skdetfsi(ip, p, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring)

    flavours = {}
    flavours["numu"] = False
    flavours["numubar"] = False # No osc numu/numubar templates
    flavours["nue"] = False
    flavours["nuebar"] = False
    for mode in apply_to_modes:
        if apply_to_modes[mode]:
            if util.IsNuMu(mode):
                flavours["numu"] = True
            if util.IsNuMuBar(mode):
                flavours["numubar"] = True
            if util.IsNuE(mode) or util.IsOscNuE(mode):
                flavours["nue"] = True
            if util.IsNuEBar(mode) or util.IsOscNuEBar(mode):
                flavours["nuebar"] = True

    has_mu = flavours["numu"] or flavours["numubar"]
    has_e = flavours["nue"] or flavours["nuebar"]
    all_mu = flavours["numu"] and flavours["numubar"]
    all_e = flavours["nue"] and flavours["nuebar"]

    if has_mu and has_e:
        if all_mu and all_e:
            return "all"
        elif flavours["numu"] and flavours["nue"] and not (flavours["numubar"] or flavours["nuebar"]):
            return "nu any"
        elif flavours["numubar"] and flavours["nuebar"]  and not (flavours["numu"] or flavours["nue"]):
            return "nubar any"
        else:
            print "Not expecting this flavor apply combination. Exiting.", flavours
            util.Die(-1)
    elif has_mu:
        if all_mu:
            return "numu any"
        else:
            len_true = 0
            for flav in flavours:
                if flavours[flav] == True:
                    len_true += 1
            if len_true == 1:
                if flavours["numu"]:
                    return "numu"
                elif flavours["numubar"]:
                    return "numubar"
                else:
                    print "ERROR: Shouldn't be here"
                    util.Die(-1)
            else:
                print "Not expecting this flavor apply combination. Exiting.", flavours
                util.Die(-1)
    elif has_e:
        if all_e:
            return "nue any"
        else:
            len_true = 0
            for flav in flavours:
                if flavours[flav] == True:
                    len_true += 1
            if len_true == 1:
                if flavours["nue"]:
                    return "nue"
                elif flavours["nuebar"]:
                    return "nuebar"
                else:
                    print "ERROR: Shouldn't be here"
                    util.Die(-1)
            else:
                print "Not expecting this flavor apply combination. Exiting.", flavours
                util.Die(-1)
    else:
        print "WARNING: Not applying param", util.ParamAsString(p), "to any true nu flavour"
        return "none"


# Turns the 3D list of apply_table[sample][mode][param] into a dictionary with keys = modes and values = whether the
# param p is applied to this mode
# Assumes apply_table has structure apply_table[sample][mode][param]
def get_apply_to_modes_xsec(ip, p, apply_table):
    param_index = None
    if util.IsRHCSystematic(p):
        param_index = 1 + ip - _first_xsec_index - conv.kNumOfXSecSyst    # 1 + because the first param in the table is RPA, which is depreciated and the BeRPA params are used instead. Iterating over the params does not give the RPA param.
    else:
        param_index = 1 + ip - _first_xsec_index

    # Initialise all modes to False
    apply_to_modes = {}
    for sample_index in range(0, conv.kns):
        s = util.GetSKSample(sample_index)
        if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:
            # Don't yet have multi-ring mu-like samples
            continue
        for mode_index in range(0, conv.knm):
            m = util.GetMode(sample_index, mode_index)
            if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                apply_to_modes[m] = False

    for sample_index in range(0, conv.kns):
        s = util.GetSKSample(sample_index)
        if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:
            # Don't yet have multi-ring mu-like samples
            continue
        for mode_index in range(0, conv.knm):
            m = util.GetMode(sample_index, mode_index)
            if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                apply_to_modes[m] = apply_to_modes[m] or apply_table[sample_index][mode_index][param_index]
    return apply_to_modes


# Turns the 3D list of apply_table[selection][mode][param] into a dictionary with keys = modes and values = whether the
# param p is applied to this mode
# Assumes apply_table has structure apply_table[selection][mode][param]
# This is specific to the SKDetFSI apply table, which groups FHC and RHC 1Rmu samples together (apply_table[0][:][:])
# and FHC and RHC 1Re samples together (apply_table[1][:][:]) and has a separate table for multi-ring samples.
# Whether the param is for FHC or RHC should be determined using util.IsRHCSystematic(p)
def get_apply_to_modes_skdetfsi(ip, p, apply_table_single_ring, apply_table_multi_ring):
    apply_to_modes = {}
    if p != conv.kNuisSuperkEnergyScale and p != conv.kRHCNuisSuperkEnergyScale:
        # Initialise all modes to False
        for sample_index in range(0, conv.kns):
            s = util.GetSKSample(sample_index)
            if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1 \
                    or (util.IsFHCSample(s) and util.IsRHCSystematic(p)) \
                    or (util.IsRHCSample(s) and util.IsFHCSystematic(p)):
                # Don't yet have multi-ring mu-like samples; there are separate SKDetFSI params for each beam current
                continue
            for mode_index in range(0, conv.knm):
                m = util.GetMode(sample_index, mode_index)
                if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                    apply_to_modes[m] = False

        for sample_index in range(0, conv.kns):
            s = util.GetSKSample(sample_index)
            if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1 \
                    or (util.IsFHCSample(s) and util.IsRHCSystematic(p)) \
                    or (util.IsRHCSample(s) and util.IsFHCSystematic(p)):
                # Don't yet have multi-ring mu-like samples
                continue
            selection_index = 0     # The SKDetFSI apply table is grouped not by sample, but by selection (0 = mu-like, 1 = e-like)
            if util.Is1ReSample(s) or util.IsMultiReSample(s):
                selection_index = 1

            param_index = None
            table = None
            if p not in MultiRingSKDetFSIParams:
                if util.IsRHCSystematic(p):
                    param_index = ip - _first_SKDetFSI_index - conv.kNumOfSKDetFSIParams
                else:
                    param_index = ip - _first_SKDetFSI_index
                table = apply_table_single_ring
            else:
                if s == conv.kSKsFCFV1RingMuLike or s == conv.kSKsFCFVMultiRingMuLike or util.IsRHCSample(s):
                    continue
                param_index = None
                if util.IsRHCSystematic(p):
                    print "ERROR: Something's gone wrong, there shouldn't be any RHC MultiRing params! Exiting"
                    util.Die(-1)
                else:
                    param_index = ip - _first_SKDetFSI_index - 2*conv.kNumOfSKDetFSIParams
                table = apply_table_multi_ring
                if selection_index == 0:
                    print "ERROR: SKDetFSI param is trying to access the MultiRmu apply table, this should not happen! Exiting. Sample =", util.SampleAsString(s)
                    util.Die(-1)

            for mode_index in range(0, conv.knm):
                m = util.GetMode(sample_index, mode_index)
                if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                    apply_to_modes[m] = apply_to_modes[m] or table[selection_index][mode_index][param_index]
    else:
        # The SK energy scale isn't in the apply tables, so put in by hand
        for sample_index in range(0, conv.kns):
            s = util.GetSKSample(sample_index)
            if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1 \
                    or (util.IsFHCSample(s) and util.IsRHCSystematic(p)) \
                    or (util.IsRHCSample(s) and util.IsFHCSystematic(p)):
                # Don't yet have multi-ring mu-like samples
                continue
            for mode_index in range(0, conv.knm):
                m = util.GetMode(sample_index, mode_index)
                if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                    apply_to_modes[m] = True
    return apply_to_modes



if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Create an SDK parameter XML file with hard-coded T2K-specific parameter information.")
  parser.add_option('-f', '--file'  , dest='filename', default="systematic_params.xml", help="Name of the output parameter config file.")
  parser.add_option('-c', '--config', dest='config'  , default="analysis_configuration.xml", help="T2K legacy config file.")
  parser.add_option("-w", "--overwrite", dest="overwrite", action="store_true", help="If this option is specified, any existing file with the same output name will be deleted and overwritten.")

  o,other = parser.parse_args()
  main(o.filename, o.config, o.overwrite)
  
  if other:
    print "Warning - ignored arguments: ",other
