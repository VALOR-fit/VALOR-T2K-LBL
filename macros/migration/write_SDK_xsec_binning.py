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

modelist_names = [
    "mulike_samples__rightsign_numu+oscnue_modes",
    "mulike_samples__wrongsign_numu+oscnue_modes",
    "mulike_samples__rightsign_nue_modes",
    "mulike_samples__wrongsign_nue_modes",
    "elike_samples__rightsign_numu+oscnue_modes",
    "elike_samples__wrongsign_numu+oscnue_modes",
    "elike_samples__rightsign_nue_modes",
    "elike_samples__wrongsign_nue_modes"
]


#______________________________________________________________________________________________________________________________

# TODO: Do we have more than one combination of kinevars for xsec binning?

def main(tag, overwrite):
    #nbins_etrue_all = []
    #nbins_theta_all = []
    edges_etrue_all = []
    edges_ereco_all = []
    edges_theta_all = []
    modes_with_same_binning = []

    # Store the xsec binning schemes
    for sample_index in range(0, conv.kns):
        s = util.GetSKSample(sample_index)
        if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:  # kSKsUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
            # Don't yet have multi-ring mu-like samples
            continue
        nbins_theta = util.NumXsecQ2Bins(s)
        edges_theta = c_to_pyroot_array_hack(util.XsecQ2BinEdges(s), nbins_theta+1) # Need to set the size of the read-write buffer, otherwise pyroot thinks it has 2**31 -1 entries
        nbins_ereco = conv.kNXsErecoBins[sample_index]
        edges_ereco = c_to_pyroot_array_hack(util.XsecErecoBinEdges(s), nbins_ereco+1)
        for mode_index in range(0, conv.knm):
            m = util.GetMode(sample_index, mode_index)
            if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                nbins_etrue = util.NumVariableXsEtrueBins(s,m)
                edges_etrue = c_to_pyroot_array_hack(util.VariableXsEtrueBinEdges(s,m), nbins_etrue+1)

                print "For sample", util.SampleAsString(s), "and mode", util.ModeAsString(m), \
                    ": \n  nbins theta =", nbins_theta, ", nbins etrue =", nbins_etrue, ", nbins ereco =", nbins_ereco, \
                    "\n  edges theta =", edges_theta, \
                    "\n  edges etrue =", edges_etrue, \
                    "\n  edges ereo =" , edges_ereco

                binning_index = get_common_binning_index(edges_etrue_all, edges_etrue, edges_theta_all, edges_theta, edges_ereco_all, edges_ereco)
                if binning_index == -1:
                    # We haven't saved this binning scheme yet, so save it and say what mode it corresponds to
                    print "Recording this binning scheme"
                    edges_theta_all.append(edges_theta) # Does this insert the a whole list into the new element, or does it append each element of the list idividually? (want the former)
                    edges_etrue_all.append(edges_etrue)
                    edges_ereco_all.append(edges_ereco)
                    modes_with_same_binning.append([[s, m]])
                else:
                    # We have already recorded this binning scheme, so just record that this mode uses this binning.
                    # To save the mode in the right position, need both the current etrue and xsec binning to match a
                    # previously added binning scheme
                    print "Binning scheme already exists, adding mode to list"
                    modes_with_same_binning[binning_index].append([s,m])
                print "\n\n"

    if not ( len(edges_etrue_all) == len(edges_theta_all)
            and len(modes_with_same_binning) == len(edges_theta_all)
            and len(edges_ereco_all) == len(edges_theta_all) ):
        print "ERROR: Number of binning schemes recorded for each kinevar / number of mode lists. Exiting."
        util.Die(-1)

    # Write the Xsec binning schemes to XML files
    for ischeme in range (0, len(edges_theta_all)):
        mode_list_str = get_mode_sample_list_str(modes_with_same_binning[ischeme])
        fname = tag + "__" + modelist_names[ischeme] + ".xml"
        print "Writing binning to XML file", fname, ", which applied to the following modes:\n", mode_list_str
        exists = os.path.isfile(fname)
        if exists and not overwrite:
            print "ERROR: Output file", fname, " already exists and the overwrite option was not specified. Exiting."
            util.Die(-1)
        f = open(fname, "w")
        write_binning_xml_content(edges_etrue_all[ischeme], edges_theta_all[ischeme], edges_ereco_all[ischeme], mode_list_str, f)

    print "Done!"
    return



def get_mode_sample_list_str(sample_mode_pair_list):
    mode_list_str = ""
    for pair in sample_mode_pair_list:
        mode_list_str += util.SampleAsString(pair[0]) + "\t\t" + util.ModeAsString(pair[1]) + "\n"
    return mode_list_str



def get_common_binning_index(edges_etrue_all, edges_etrue, edges_theta_all, edges_theta, edges_ereco_all, edges_ereco):
    etrue_indicies = get_indicies_of_matching_list(edges_etrue_all, edges_etrue)
    ereco_indicies = get_indicies_of_matching_list(edges_ereco_all, edges_ereco)
    theta_indicies = get_indicies_of_matching_list(edges_theta_all, edges_theta)
    if len(etrue_indicies) == 0 or len(theta_indicies) == 0 or len(ereco_indicies) == 0:
        return -1   # This binning scheme has not yet been recorded
    common_index = -1
    found = False
    for etrue_index in etrue_indicies:
        for theta_index in theta_indicies:
            for ereco_index in ereco_indicies:
                if etrue_index == theta_index and etrue_index == ereco_index:
                    if not found:
                        common_index = etrue_index
                        found = True
                    else:
                        print "ERROR: Should not be more than one common index. Exiting."
                        util.Die(-1)
    return common_index



def get_indicies_of_matching_list(edges_all, edges_search):
    matching_indicies = []
    for i in range(0, len(edges_all)):
        if edges_all[i] == edges_search:
            matching_indicies.append(i)
    return matching_indicies



def write_binning_xml_content(edges_etrue, edges_theta, edges_ereco, mode_list_str, f):
    f.write( "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n\n")
    f.write( "<!--\nTODO: Make sure this xsec binning file is only applied to the SDK equivalent of the following legacy T2K modes:\n")
    f.write( mode_list_str + "-->\n\n")
    f.write( "<binning>\n\n")
    f.write( "\t<space>\n\t\tEtrue,Ereco,Thetareco\n\t</space>\n\n")
    f.write( "\t<type>\n\t\tedge-to-edge\n\t</type>\n\n")
    write_single_kinevar_e2e_binning("Etrue", edges_etrue, f)
    write_single_kinevar_e2e_binning("Ereco", edges_ereco, f)
    write_single_kinevar_e2e_binning("thetareco", edges_theta, f)
    f.write( "</binning>\n")



def write_single_kinevar_e2e_binning(name, edges, f):
    print_line_length = 5
    f.write( "\t<"+name+"_binning>\n")
    f.write( "\t\t<bin_size> variable </bin_size>\n")
    f.write( "\t\t<nbins> " + str(len(edges)-1) + " </nbins>\n")
    f.write( "\t\t<bin_edges>\n")
    edges_str = "\t\t\t"
    for i in range(0, len(edges)):
        edges_str += str(edges[i])
        if i != len(edges) -1:
            edges_str += ",\t"
            if i != 0 and i % print_line_length == 0:
                edges_str += "\n\t\t\t"
    f.write( edges_str + "\n")
    f.write( "\t\t</bin_edges>\n")
    f.write( "\t</"+name+"_binning>\n\n")



if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Create an SDK binning XML file for T2K xsec binning schemes.")
  parser.add_option('-t', '--tag', dest='tag', default="t2k_xsec_binning", help="Tag used in naming the output XML files")
  parser.add_option('-o', '--overwrite', dest='overwrite', action="store_true", help="If this option is specified, any existing file with the same output name will be deleted and overwritten.")

  o,other = parser.parse_args()
  main(o.tag, o.overwrite)
  
  if other:
    print "Warning - ignored arguments: ",other
