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


_default_root_node_name = "systematic_list"
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

# Stupid pyroot can't read many 1D C arrays, so I got this array in cint using
# BANFF* banff = skexp->GetBANFF()
# const conv::PhysParamId_t * banff_pars = banff->GetParameters()
_banff_pars = [
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    201,
    202,
    203,
    204,
    205,
    206,
    207,
    208,
    209,
    210,
    211,
    212,
    213,
    214,
    215,
    216,
    217,
    218,
    219,
    220,
    221,
    222,
    223,
    224,
    225,
    84,
    85,
    86,
    76,
    77,
    87,
    88,
    89,
    90,
    91,
    92,
    93,
    94,
    95,
    96,
    97,
    98,
    99,
    100,
    101,
    102,
    103,
    104,
    105
]

# Stupid pyroot can't read many 1D C arrays, so I got this array in cint using
# SuperkDetectorFSI* skdetfsi = skexp->GetSuperkDetectorFSI()
# const conv::PhysParamId_t * sk_pars = skdetfsi->GetParameters()
_skdetfsi_pars = [
    132,
    133,
    134,
    135,
    136,
    137,
    138,
    139,
    140,
    141,
    142,
    143,
    144,
    145,
    146,
    147,
    148,
    149,
    302,
    303,
    304,
    305,
    306,
    307,
    308,
    309,
    310,
    311,
    312,
    313,
    314,
    315,
    316,
    317,
    318,
    319,
    345,
    346,
    348,
    349,
    351,
    352,
    354,
    355,
    31#,
    # 0,
    # 369,
    # 0,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
    # 0,
    # -1,
]

#______________________________________________________________________________________________________________________________


def main(tag, t2k_legacy_config, overwrite):
    skexp = valor.SuperkExpectation()
    skexp.Read(os.environ[env_var] + "/config/" + t2k_legacy_config)
    set_firsts()
    out_filename_banff = tag.replace("[NAME]", "banff")
    out_filename_sk    = tag.replace("[NAME]", "sk")
    write_input_param_list(skexp, out_filename_banff, overwrite, True)
    write_input_param_list(skexp, out_filename_sk   , overwrite, False)



def write_input_param_list(skexp, outfilename, overwrite, is_banff):
    global _banff_pars
    global _skdetfsi_pars

    # Check if the output file already exists
    exists = os.path.isfile(outfilename)
    if not exists:
        print "Creating new file", outfilename
    else:
        if overwrite and outfilename != "":
            print "Overwriting file", outfilename
            os.remove(outfilename)
        else:
            print "A file named", outfilename, "already exists and the --overwrite option was not specified. Exiting."
            util.Die(-1)

    # Set up the XML file
    xmleng = ROOT.TXMLEngine()
    xmldoc = xmleng.NewDoc()
    if not xmldoc:
        print "FATAL: Could not create the XML document in file", outfilename, ". Exiting."
        util.Die(-1)

    # Create a root node if one doesn't exist already
    root_node_attr_name = "banff"
    if not is_banff:
        root_node_attr_name = "sk"
    root_node = try_make_root_node(xmleng, xmldoc, root_node_attr_name+"_Run1-9")

    params = []
    param_names = []
    if is_banff:
        #params = c_to_pyroot_array_hack(skexp.GetBANFF().GetParameters(), skexp.GetBANFF().Size())   # pyroot can't handle c++ arrays, see c_to_pyroot_array_hack definition for details
        params = _banff_pars
    else:
        #params = c_to_pyroot_array_hack(skexp.GetSuperkDetectorFSI().GetParameters(), skexp.GetSuperkDetectorFSI().Size)
        params = _skdetfsi_pars

    # Read the attributes of each T2K systematic parameter and write to the XML file
    for p in params:
        simple_name = "\t" + get_name(p) + "\t"
        param_names.append(simple_name)
        xmleng.NewChild(root_node, 0, "parameter", simple_name)

    print "For", root_node_attr_name, ", param list =\n", param_names

    xmleng.SaveDoc(xmldoc, outfilename)
    xmleng.FreeDoc(xmldoc)
    del xmleng
    del xmldoc

    print "Successfully wrote the parameter list file", outfilename
    print "Done!"
    return



def get_name(p):
    global _name_flux_counter
    global _name_flux_counter_rhc
    name = ""
    if util.IsBANFF(p):
        if util.IsFlux(p):
            # Use a more descriptive name for flux params than e.g. banff00.
            flavour_str = get_flavour_flux(p)
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


def get_flavour_flux(p):
    ip = get_param_index(p)
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
    print "ERROR: Tried to get the flavour of a flux param, but input a non-flux param", util.ParamAsString(p)
    util.Die(-1)
    return "none"



# Check if a root node already exists. If not, create one.
def try_make_root_node(xmleng, xmldoc, list_name):
    root_node = xmleng.DocGetRootElement(xmldoc)
    try:
        xmleng.IsEmptyNode(root_node)  # If there's no ROOT node, TMXLEngine will always throw an error
        print "ERROR: A root node already existed. This should not happen!"
    except:
        root_node = xmleng.NewChild(0, 0, _default_root_node_name)
        xmleng.NewAttr(root_node, 0, "name", list_name)
        xmleng.DocSetRootElement(xmldoc, root_node)
        print "No root node found. Creating a root node with attribute \'name\' =", list_name
    return root_node


def get_param_index(p):
    params = valor.util.NamedParameterList("nuosc3f+nuisance")
    for ip in range(0, params.GetNumberParameters()):
        p_try = params.WhichParam(ip)
        if p_try == p:
            return ip
    return -1


# Work out indexes of the first xsex, flux, and SKDetFSI params
def set_firsts():
    global _first_xsec_index
    global _first_flux_index
    global _first_SKDetFSI_index
    params = valor.util.NamedParameterList("nuosc3f+nuisance")
    for ip in range(0, params.GetNumberParameters()):
        p = params.WhichParam(ip)
        if util.IsBANFF(p):
            if util.IsFlux(p):
                if _first_flux_index < 0:
                    _first_flux_index = ip
            elif _first_xsec_index < 0:
                _first_xsec_index = ip
        elif _first_SKDetFSI_index < 0:
            _first_SKDetFSI_index = ip


if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Create an SDK parameter XML file with hard-coded T2K-specific parameter information.")
  parser.add_option('-t', '--tag'      , dest='tag'      , default="t2k_[NAME]_params_run1-9d.xml", help="Name of the output parameter config file.")
  parser.add_option('-c', '--config'   , dest='config'   , default="analysis_configuration.xml"   , help="T2K legacy config file.")
  parser.add_option("-w", "--overwrite", dest="overwrite", action="store_true"                    , help="If this option is specified, any existing file with the same output name will be deleted and overwritten.")

  o,other = parser.parse_args()
  main(o.tag, o.config, o.overwrite)
  
  if other:
    print "Warning - ignored arguments: ",other
