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


_xmleng = 0
_xmldoc = 0
_default_root_node_name = "param_list"
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

_param_fields = namedtuple('param_fields',
   [
       'name_simple',
       'name_long',
       'name_rootex',
       'name_tex',
       'kine_ranges',
       'type',
       'category',
       'apply_detector',
       'apply_beam',
       'apply_selection',
       'apply_flav',
       'apply_osc',
       'apply_reac',
       'apply_target',
       'limits'
   ])


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
    global _xmleng
    global _xmldoc
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

    # Set up the XML file
    _xmleng = ROOT.TXMLEngine()
    _xmldoc = _xmleng.NewDoc()
    if not _xmldoc:
        print "FATAL: Could not create the XML document in file", filename, ". Exiting."
        util.Die(-1)

    # Create a root node if one doesn't exist already
    root_node = try_make_root_node("T2K_Run1-9")

    skexp = valor.SuperkExpectation()
    skexp.Read(os.environ[env_var] + "/config/" + t2k_legacy_config)

    # pyroot can't handle c++ arrays, see c_to_pyroot_array_hack definition for details
    apply_table_xsec               = c_to_pyroot_array_hack(conv.kXsApplyXSecSystToSKMode       , conv.kns, conv.knm, conv.kNumOfXSecSyst)                 # apply_table[samples][modes][xsec params] = bool
    apply_table_SKDetFSI_FHC       = c_to_pyroot_array_hack(conv.kSKDetFSIParamToSKMode         , 2       , conv.knm, conv.kNumOfSKDetFSIParams)           # apply_table[mu-like, or e-like][modes][skdetfsi params] = bool
    apply_table_SKDetFSI_RHC       = c_to_pyroot_array_hack(conv.kRHCSKDetFSIParamToSKMode      , 2       , conv.knm, conv.kNumOfRHCSKDetFSIParams)        # apply_table[mu-like, or e-like][modes][skdetfsi params] = bool
    apply_table_SKDetFSI_multiring = c_to_pyroot_array_hack(conv.kMultiRingSKDetFSIParamToSKMode, 2       , conv.knm, conv.kNumOfSKDetFSIParamsMultiRing)  # apply_table[mu-like, or e-like][modes][skdetfsi params] = bool

    # lists of apply mask types (used for printing)
    reac_types_list = []
    selec_list      = []
    flav_list       = []

    # Read the attributes of each T2K systematic parameter and write to the XML file
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
        apply_table_SKDetFSI = apply_table_SKDetFSI_FHC
        if util.IsRHCSystematic(p):
            apply_table_SKDetFSI = apply_table_SKDetFSI_RHC

        # Set all of the possible XML fields
        fields = _param_fields(
            name_simple     = get_name(ip, p),
            name_long       = util.ParamAsString(p),
            name_rootex     = util.ParamAsStringPaperRooTex(p),
            name_tex        = util.AsStringTex(p),
            kine_ranges     = {},
            type            = get_type(ip, p),
            category        = [],
            apply_detector  = "SuperK",
            apply_beam      = get_beam_mode(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring),
            apply_selection = get_selection(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring),
            apply_flav      = get_flavour(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring),
            apply_osc       = get_osc(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring),
            apply_reac      = get_reac_type(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring),
            apply_target    = get_target(p),
            limits          = get_limits(p)
        )
        set_kine_range(ip, p, fields)
        set_category(p, fields)
        print "For param", util.ParamAsString(p), "the following fields have been set:\n", fields

        # Collect additional info
        reac_types_list.append(fields.apply_reac)
        selec_list.append(fields.apply_selection)
        flav_list.append(fields.apply_flav)
        reac_types_all_okay = reac_types_all_okay and check_reac_types(ip, p, apply_table_xsec, apply_table_SKDetFSI)

        # Write all parameter fields to the output XML file
        add_param(root_node, fields)

        del fields
        print "\n\n"

    _xmleng.SaveDoc(_xmldoc, filename)
    _xmleng.FreeDoc(_xmldoc)

    print "The following reaction type apply masks are used. Make sure they exist in FitParamApply\n", list(set(reac_types_list)), "\n"
    print "The following selection apply masks are used. Make sure they exist in FitParamApply\n", list(set(selec_list)), "\n"
    print "The following nu flavour apply masks are used. Make sure they exist in FitParamApply\n", list(set(flav_list)), "\n"

    if not reac_types_all_okay:
        print "WARNING: At least one of the parameters have different applied reaction types for different samples. This may be a problem for specifying parameters in this way."

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
        # TODO: use a more descriptive name for SKDetFSI params than e.g. skdetfsi000?
        name = util.AsStringCompact(p)
    # Many param names start with 'f_', remove this
    name = name.replace("f_", "")
    if not name:
        print "ERROR: name not found for param", p, ". Exiting."
        util.Die(-1)
    return name


def get_type(ip, p):
    global _first_xsec_index
    if p == conv.kNuisSuperkEnergyScale or p == conv.kRHCNuisSuperkEnergyScale:
        return "lateral"
    elif util.IsBANFF(p) and not util.IsFlux(p):
        # Xsec params could have splines
        if util.IsRHCSystematic(p):
            index = 1 + ip - _first_xsec_index - conv.kNumOfXSecSyst    # 1 + because the first param in the table is RPA, which is depreciated and the BeRPA params are used instead. Iterating over the params does not give the RPA param.
            if conv.kRHCXSecSystHasWghtSplines[index]:
                return "response"
            else:
                return "normalisation"
        else:
            index = 1 + ip - _first_xsec_index  # 1 + because the first param in the table is RPA, which is depreciated and the BeRPA params are used instead. Iterating over the params does not give the RPA param.
            if conv.kXSecSystHasWghtSplines[index]:
                return "response"
            else:
                return "normalisation"
    else:
        return "normalisation"


def get_beam_mode(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring):
    if util.IsBANFF(p):
        if util.IsFlux(p):
            if util.IsRHCSystematic(p):
                return "RHC"
            else:
                return "FHC"
        else:
            # Use the 3D array kXsApplyXSecSystToSKMode to decide whether each xsec param is applied to each sample and therefore beam current
            apply_to_samples = get_apply_to_samples_xsec(ip, p, apply_table_xsec)
            apply_FHC     = False
            apply_RHC     = False
            for s in apply_to_samples:
                if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:
                    continue
                if apply_to_samples[s]:
                    if util.IsFHCSample(s):
                        apply_FHC = True
                    else:
                        apply_RHC = True
            # Get beam current
            if apply_FHC and apply_RHC:
                return "all"
            elif apply_FHC:
                return "FHC"
            elif apply_RHC:
                return "RHC"
            else:
                return "none"
    else:
        apply_to_selections = get_apply_to_selections_skdetfsi(ip, p, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring)
        apply_to_any_selection = apply_to_selections["1Rmu"] or apply_to_selections["1Re"] or apply_to_selections["MultiRe"]
        if apply_to_any_selection:
            if util.IsEnergyScale(p):
                # SKDetFSI params are grouped into FHC params, RHC params and Escale, so the only param to be applied
                # to both FHC and RHC is Escale
                return "all"
            elif util.IsRHCSystematic(p):
                return "RHC"
            else:
                return "FHC"
        else:
            return "none"


def get_selection(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring):
    if util.IsBANFF(p):
        if util.IsFlux(p):
            return "all"
        else:
            # Use the 3D array kXsApplyXSecSystToSKMode to decide whether each xsec param is applied to each sample
            apply_to_samples = get_apply_to_samples_xsec(ip, p, apply_table_xsec)
            apply_1Rmu    = False
            apply_1Re     = False
            apply_MultiRe = False
            for s in apply_to_samples:
                if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:
                    continue
                if apply_to_samples[s]:
                    if util.Is1RMuSample(s):
                        apply_1Rmu = True
                    elif util.Is1ReSample(s):
                        apply_1Re = True
                    else:
                        apply_MultiRe = True
            # Get selection
            if apply_1Rmu and apply_1Re and apply_MultiRe:
                return "all rings"
            elif apply_1Rmu and not (apply_1Re or apply_MultiRe):
                return "1Rmu"
            elif apply_1Re and apply_MultiRe and not apply_1Rmu:
                return "e-like rings"
            elif apply_1Re and not apply_MultiRe:
                return "1Re"
            elif apply_MultiRe and not apply_1Re:
                return "MultiRe"
            elif apply_1Rmu and apply_1Re and not apply_MultiRe:
                return "1R any"
            elif not apply_1Rmu and not apply_1Re and not apply_MultiRe:
                return "none"
            else:
                print "Shouldn't be here!"
                util.Die(-1)
    else:
        apply_to_selections = get_apply_to_selections_skdetfsi(ip, p, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring)
        if apply_to_selections["1Rmu"] and apply_to_selections["1Re"] and apply_to_selections["MultiRe"]:
            return "all"
        elif apply_to_selections["MultiRe"]:
            return "MultiRe"
        elif apply_to_selections["1Rmu"] and apply_to_selections["1Re"]:
            return "1R_all"
        elif apply_to_selections["1Rmu"]:
            return "1Rmu"
        elif apply_to_selections["1Re"]:
            return "1Re"
        else:
            return "none"
    print "Shouldn't be here!"
    util.Die(-1)


# Turns the 3D list of apply_table[sample][mode][param] into a dictionary with keys = samples and values = whether the
# param p is applied to this sample ( = is applied to any of the modes)
# Assumes apply_table has structure apply_table[sample][mode][param]
def get_apply_to_samples_xsec(ip, p, apply_table):
    apply_to_samples = {}
    param_index = None
    if util.IsRHCSystematic(p):
        param_index = 1 + ip - _first_xsec_index - conv.kNumOfXSecSyst    # 1 + because the first param in the table is RPA, which is depreciated and the BeRPA params are used instead. Iterating over the params does not give the RPA param.
    else:
        param_index = 1 + ip - _first_xsec_index
    for sample_index in range(0, conv.kns):
        s = util.GetSKSample(sample_index)
        if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:
            # Don't yet have multi-ring mu-like samples
            continue
        applies_to_this_sample = False
        for param_apply_list_for_mode in apply_table[sample_index]:
            applies_to_this_sample = applies_to_this_sample or param_apply_list_for_mode[param_index]
        apply_to_samples[s] = applies_to_this_sample
    return apply_to_samples


# Turns the 3D list of apply_table[selection][mode][param] into a dictionary with keys = samples and values = whether the
# param p is applied to this selection ( = is applied to any of the modes)
# Assumes apply_table has structure apply_table[selection][mode][param]
# This is specific to the SKDetFSI apply table, which groups FHC and RHC 1Rmu samples together (apply_table[0][:][:])
# and FHC and RHC 1Re samples together (apply_table[1][:][:]) and has a separate table for multi-ring samples.
# Whether the param is for FHC or RHC should be determined using util.IsRHCSystematic(p)
def get_apply_to_selections_skdetfsi(ip, p, apply_table_single_ring, apply_table_multi_ring):
    param_index = None
    apply_1Rmu_like = False
    apply_1Re_like = False
    apply_MultiRe_like = False
    if p != conv.kNuisSuperkEnergyScale and p != conv.kRHCNuisSuperkEnergyScale:
        if p not in MultiRingSKDetFSIParams:
            if util.IsRHCSystematic(p):
                param_index = ip - _first_SKDetFSI_index - conv.kNumOfSKDetFSIParams
            else:
                param_index = ip - _first_SKDetFSI_index
            for param_apply_list_for_mode in apply_table_single_ring[0]:
                apply_1Rmu_like = apply_1Rmu_like or param_apply_list_for_mode[param_index]
            for param_apply_list_for_mode in apply_table_single_ring[1]:
                apply_1Re_like = apply_1Re_like or param_apply_list_for_mode[param_index]
        else:
            if util.IsRHCSystematic(p):
                print "ERROR: Something's gone wrong, there shouldn't be any RHC MultiRing params! Exiting"
                util.Die(-1)
            else:
                param_index = ip - _first_SKDetFSI_index - 2*conv.kNumOfSKDetFSIParams
            for param_apply_list_for_mode in apply_table_multi_ring[1]: # skip multi-ring mu-like sample
                apply_MultiRe_like = apply_MultiRe_like or param_apply_list_for_mode[param_index]
    else:
        # The SK energy scale isn't in the apply tables, so put in by hand
        apply_1Rmu_like = True
        apply_1Re_like = True
        apply_MultiRe_like = True
    apply_to_selections = {}
    apply_to_selections["1Rmu"]    = apply_1Rmu_like
    apply_to_selections["1Re"]     = apply_1Re_like
    apply_to_selections["MultiRe"] = apply_MultiRe_like
    return apply_to_selections


def get_reac_type(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring):
    global _possible_reac_type_names
    apply_to_modes = None
    if util.IsBANFF(p):
        if util.IsFlux(p):
            return "all"
        else:
            # Use the 3D array kXsApplyXSecSystToSKMode to decide whether each xsec param is applied to each sample
            apply_to_modes = get_apply_to_modes_xsec(ip, p, apply_table_xsec)
    else:
        apply_to_modes = get_apply_to_modes_skdetfsi(ip, p, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring)

    reac_types = []
    is_applied_to_all = True
    for mode in apply_to_modes:
        is_applied_to_all = is_applied_to_all and apply_to_modes[mode]
    if is_applied_to_all:
        reac_types.append("all")
    else:
        for mode in apply_to_modes:
            if apply_to_modes[mode]:
                full_mode_name = util.ModeAsString(mode)
                for reac_type_name in _possible_reac_type_names:
                    if reac_type_name in full_mode_name:
                        reac_types.append(reac_type_name)
    reac_types = list(set(reac_types))  # remove duplicate entries

    if len(reac_types) == 0:
        return "none"
    elif len(reac_types) == 1:
        return reac_types[0]
    elif "CC QE" in reac_types and "CC coh" in reac_types and "CC oth" in reac_types and "CC 1pi" in reac_types \
            and "CC 2p2h" in reac_types and len(reac_types) == 5:
        return "all CC"
    elif "NC oth" in reac_types and "NC 1gamma" in reac_types and "NC 1pi0" in reac_types and "NC coh" in reac_types \
            and "NC 1piC" in reac_types and len(reac_types) == 5:
        return "all NC"
    elif "NC 1pi0" in reac_types and "CC 1pi" in reac_types and "NC 1piC" in reac_types and len(reac_types) == 3:
        return "all 1pi"    # TODO: does CC 1pi include both charged and neutral pions? ktApplyToAll1pi assumes both charged and neutral
    elif "CC QE" in reac_types and "CC 2p2h" in reac_types and len(reac_types) == 2:
        return "CC QE like"
    elif "CC coh" in reac_types and "CC 1pi" in reac_types and "CC oth" in reac_types and len(reac_types) == 3:
        return "CC nQE"
    else:
        print "ERROR: unrecognised combination of reaction types:", reac_types
        util.Die(-1)


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


def get_osc(ip, p, apply_table_xsec, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring):
    apply_to_modes = None
    if util.IsBANFF(p):
        if util.IsFlux(p):
            return "all"
        else:
            # Use the 3D array kXsApplyXSecSystToSKMode to decide whether each xsec param is applied to each sample
            apply_to_modes = get_apply_to_modes_xsec(ip, p, apply_table_xsec)
    else:
        apply_to_modes = get_apply_to_modes_skdetfsi(ip, p, apply_table_SKDetFSI, apply_table_SKDetFSI_multiring)

    masks = {}
    masks["osc"] = False
    masks["unosc"] = False # No osc numu/numubar templates
    for mode in apply_to_modes:
        if apply_to_modes[mode]:
            if util.IsNuMu(mode) or util.IsNuMuBar(mode) or util.IsNuE(mode) or util.IsNuEBar(mode):
                masks["unosc"] = True
            if util.IsOscNuEBar(mode) or util.IsOscNuE(mode):
                masks["osc"] = True

    if masks["osc"] and masks["unosc"]:
        return "all"
    elif masks["osc"]:
        return "osc"
    elif masks["unosc"]:
        return "unosc"
    else:
        print "WARNING: Not applying param", util.ParamAsString(p), "to any mode"
        return "none"



def get_target(p):
    # TODO: For now the SDK doesn't use the nuclear target, so just set to apply to all
    return "all"
    # if p == conv.kNuisBANFFNorm2p2h or p == conv.kNuisBANFFFermiMom or p == conv.kNuisBANFFNormCCcoh \
    #         or p == conv.kNuisBANFFNorm2p2hBar or p == conv.kNuisBANFFShape2p2h or p == conv.kNuisBANFFNorm2p2hCtoO \
    #         or p == conv.kRHCNuisBANFFFermiMom or p == conv.kRHCNuisBANFFNormCCcoh \
    #         or p == conv.kRHCNuisBANFFShape2p2h or p == conv.kRHCNuisBANFFNorm2p2hCtoO:
    #     return "Oxygen"
    # else:
    #     return "Oxygen, Hydrogen"   # TODO: make FitParamApply work with this


def get_limits(p):
    lims = valor.T2KPhysicalParameterLimits.Instance()
    lims.SetDefaults()
    min = lims.MinValue(p)
    max = lims.MaxValue(p)
    #if min == -9999 or max == 9999:
    #    print "ERROR: Should not get the limits [", min, ",", max, "] for param", get_name(p)
    limits_str = ""
    if min == -9999:
        limits_str += "-inf"
    else:
        limits_str += str(min)
    limits_str += ", "
    if max == 9999:
        limits_str += "inf"
    else:
        limits_str += str(max)
    return limits_str


    # Get kinematic ranges of parameter
def set_kine_range(ip, p, fields):
    global _first_flux_index
    global _first_SKDetFSI_index
    range_low  = -999
    range_high = -999
    if util.IsBANFF(p):
        if util.IsFlux(p):
            if util.IsRHCSystematic(p):
                index = ip - _first_flux_index - conv.kNumOfFluxParamsInBANFFfit # - num FHC params because FHC and RHC ranges are in the same indexes of 2 different arrays
                ranges = c_to_pyroot_array_hack(conv.kEtrueRangeOfRHCFluxParamInBANFFfit, conv.kNumOfRHCParamsInBANFFfit, 2)   # pyroot can't handle c++ arrays, see c_to_pyroot_array_hack definition for details
                range_low  = ranges[index][0]
                range_high = ranges[index][1]
            else:
                index = ip - _first_flux_index
                ranges = c_to_pyroot_array_hack(conv.kEtrueRangeOfFluxParamInBANFFfit, conv.kNumOfParamsInBANFFfit, 2)
                range_low  = ranges[index][0]
                range_high = ranges[index][1]
            if range_low != -999 and range_high != -999:
                fields.kine_ranges["Etrue"] = str(range_low) + ", " + str(range_high)
        else:
            # interaction params
            fields.kine_ranges["Etrue"] = "0, inf"
    elif util.IsSKDetFSI(p):
        if not (p == conv.kNuisSuperkEnergyScale or p == conv.kRHCNuisSuperkEnergyScale):
            ranges = None
            index  = None
            if p not in MultiRingSKDetFSIParams:
                # Param is not MultiReSKDetFSI
                if util.IsRHCSystematic(p):
                    index = ip - _first_SKDetFSI_index - conv.kNumOfSKDetFSIParams
                    ranges = c_to_pyroot_array_hack(conv.kRHCSKDetFSIParamsErecoCut, conv.kNumOfRHCSKDetFSIParams, 2)
                else:
                    index = ip - _first_SKDetFSI_index
                    ranges = c_to_pyroot_array_hack(conv.kSKDetFSIParamsErecoCut, conv.kNumOfSKDetFSIParams, 2)
            else:
                # Param is MultiReSKDetFSI
                index = ip - _first_SKDetFSI_index - 2*conv.kNumOfSKDetFSIParams
                if util.IsRHCSystematic(p):
                    print "ERROR: Something's gone wrong, there shouldn't be any RHC MultiRing params! Exiting"
                    util.Die(-1)
                ranges = c_to_pyroot_array_hack(conv.kMultiRingSKDetFSIParamsErecoCut, conv.kNumOfSKDetFSIParamsMultiRing, 2)
            range_low  = ranges[index][0]
            range_high = ranges[index][1]
            if range_low != -999 and range_high != -999:
                fields.kine_ranges["Ereco"] = str(range_low) + ", " + str(range_high)
        else:
            # SK energy scale
            fields.kine_ranges["Etrue"] = "0, inf"


def set_category(p, fields):
    if util.IsBANFF(p):
        fields.category.append("BANFF")
        if util.IsFlux(p):
            fields.category.append("flux")
        else:
            fields.category.append("int")
    else:
        fields.category.append("SKDetFSI")


# Check if a root node already exists. If not, create one.
def try_make_root_node(list_name):
    root_node = _xmleng.DocGetRootElement(_xmldoc)
    try:
        _xmleng.IsEmptyNode(root_node)  # If there's no ROOT node, TMXLEngine will always throw an error
        print "ERROR: A root node already existed. This should not happen!"
    except:
        root_node = _xmleng.NewChild(0, 0, _default_root_node_name)
        _xmleng.NewAttr(root_node, 0, "name", list_name)
        _xmleng.DocSetRootElement(_xmldoc, root_node)
        print "No root node found. Creating a root node with attribute \'name\' =", list_name
    return root_node


def add_param(root_node, fields):
    # Check whether a parameter with this simple name already exists in the file
    check_param_exists(root_node, fields.name_simple)

    # This structure for a parameter is according to Parameter.h
    param_node        = add_node(root_node , "param"         , ""                     , False, True )
    name_node         = add_node(param_node, "name"          , ""                     , True , True )
    simple_name_node  = add_node(name_node , "simple"        , fields.name_simple     , False, False)
    long_name_node    = add_node(name_node , "long"          , fields.name_long       , False, False)
    rootex_name_node  = add_node(name_node , "rootex"        , fields.name_rootex     , False, False)
    tex_name_node     = add_node(name_node , "tex"           , fields.name_tex        , False, False)
    range_node        = add_node(param_node, "apply_to_range", ""                     , True , True )
    for kv in fields.kine_ranges:
        kine_node     = add_node(range_node, kv              , fields.kine_ranges[kv] , False, False)
    type_node         = add_node(param_node, "type"          , fields.type            , True , False)
    cats_node         = add_node(param_node, "categories"    , ""                     , True , True )
    if len(fields.category) > 0:
        for category in fields.category:
            cat_node  = add_node(cats_node , "category"      , category               , False, False)
    else:
        cat_node      = add_node(cats_node , "category"      , ""                     , False, False)
    apply_node        = add_node(param_node, "apply_to"      , ""                     , True , True )
    apply_det_node    = add_node(apply_node, "detector"      , fields.apply_detector  , False, False)
    apply_beam_node   = add_node(apply_node, "beam"          , fields.apply_beam      , False, False)
    apply_selec_node  = add_node(apply_node, "selection"     , fields.apply_selection , False, False)
    apply_reac_node   = add_node(apply_node, "reac_type"     , fields.apply_reac      , False, False)
    apply_flav_node   = add_node(apply_node, "flav"          , fields.apply_flav      , False, False)
    apply_flav_node   = add_node(apply_node, "osc"           , fields.apply_osc       , False, False)
    apply_target_node = add_node(apply_node, "target"        , fields.apply_target    , False, False)
    limits_node       = add_node(param_node, "limits"        , fields.limits          , False, False)
    _xmleng.AddRawLine(root_node, "\n\n")     # Additional blank lines after each parameter, for readability


def check_param_exists(root_node, simple_name):
    all_param_nodes = util.FindAllNodes(_xmleng, root_node)
    for param_node in all_param_nodes:
        if util.NodeExists(_xmleng, param_node, "name/simple"):
            param_name = util.XML2String(_xmleng, param_node, "name/simple")
            param_name_check = util.TrimSpaces(param_name)
            if param_name_check == simple_name:
                print "FATAL: A parameter with simple name", simple_name, "already exists. Simple names must be unique. Exiting"
                util.Die(-1)


# Improve readability and make it easier to manaully insert node content (changes e.g. <node\> to <node>  <\node>)
def add_node(parent_node, node_name, node_content_in, blank_line_after = False, should_always_be_empty = False):
    node_content_print = ""
    if not should_always_be_empty:
        node_content_print = "  " + node_content_in + "  "
    new_node = _xmleng.NewChild(parent_node, 0, node_name, node_content_print)
    if blank_line_after:
        _xmleng.AddRawLine(parent_node, "")  # spacing between groups of nodes
    return new_node


# For a given parameter, get the applied reac type(s) for each sample and check whether this is consistent across
# all samples. If different samples have different applied reac types, this will be a problem for this
# parameter specification, as the applied sample(s) is specified independently of the applied reac type(s).
# Assumes apply_table has structure apply_table[sample][mode][param]
def check_reac_types(ip, p, apply_table_xsec, apply_table_skdetfsi_single_ring):
    #print "Checking whether we can represent the reac type apply mask independent of the sample (for param", util.ParamAsString(p), ")"
    global _possible_reac_type_names
    param_index = None
    apply_table = None

    if util.IsBANFF(p):
        if util.IsFlux(p):
            return True     # Flux params are applied to all reac types, so are okay
        else:
            if util.IsRHCSystematic(p):
                param_index = 1 + ip - _first_xsec_index - conv.kNumOfXSecSyst    # 1 + because the first param in the table is RPA, which is depreciated and the BeRPA params are used instead. Iterating over the params does not give the RPA param.
            else:
                param_index = 1 + ip - _first_xsec_index
            apply_table = apply_table_xsec
    elif util.IsSKDetFSI(p):
        if p == conv.kNuisSuperkEnergyScale or p == conv.kRHCNuisSuperkEnergyScale:
            return True     # the SK E-scale is just applied to everything, so is okay
        if p not in MultiRingSKDetFSIParams:
            if util.IsRHCSystematic(p):
                param_index = ip - _first_SKDetFSI_index - conv.kNumOfSKDetFSIParams
            else:
                param_index = ip - _first_SKDetFSI_index
            apply_table = apply_table_skdetfsi_single_ring
        else:
            return True     # because only applies to one sample

    all_okay = True
    for mode_index in range(0, conv.knm):
        apply_reac_types_to_sample = {}
        for sample_index in range(0, conv.kns):
            s = util.GetSKSample(sample_index)
            m = util.GetMode(sample_index, mode_index)
            if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1 \
                    or m == conv.kSKmUndefined or m == 2**64 -1:  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                # Don't yet have multi-ring mu-like samples
                continue

            # The xsec and SKDetFSI tables use different first indexes. For xsec use 1st index = sample, for SKDetFSI use 1st index = 0 for mu-like, 1 or e-like
            apply_table_first_index = sample_index
            if util.IsSKDetFSI(p):
                if p not in MultiRingSKDetFSIParams:
                    if util.Is1RMuSample(s):
                        apply_table_first_index = 0
                    else:
                        apply_table_first_index = 1
                else:
                    if util.Is1RMuSample(s) or util.Is1ReSample(s):
                        continue
                    else:
                        apply_table_first_index = 1

            reac_types = []
            full_mode_name = util.ModeAsString(m)
            for reac_type_name in _possible_reac_type_names:
                if reac_type_name in full_mode_name and apply_table[apply_table_first_index][mode_index][param_index]:
                    reac_types.append(reac_type_name)
            apply_reac_types_to_sample[s] = reac_types.sort()

        first = True
        first_sample_apply = None
        for sample_index in range(0, conv.kns):
            s = util.GetSKSample(sample_index)
            m = util.GetMode(sample_index, mode_index)
            if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1 \
                    or m == conv.kSKmUndefined or m == 2**64 -1:  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                # Don't yet have multi-ring mu-like samples
                continue
            if first:
                first_sample_apply = apply_reac_types_to_sample[s]
                first = False
            elif apply_reac_types_to_sample[s] != first_sample_apply:
                print "WARNING: For param", util.ParamAsString(p), " and mode", util.ModeAsString(m), ", the applied reac types are not the same for each sample:\n", apply_reac_types_to_sample
                all_okay = False
                break
    #if all_okay:
        #print "All applied reac types are okay for all samples used by this parameter"
    return all_okay


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
