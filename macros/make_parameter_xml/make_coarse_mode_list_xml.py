#!/usr/bin/python

import ROOT
import os
from collections import namedtuple
import sys

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
es    = valor.ExperimentSetup.Instance()

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

python_undefined = 2L*sys.maxint+1
python_undefined_short = 2**32 -1

_xmleng = 0
_xmldoc = 0

_param_fields = namedtuple('param_fields',
                           [
                               'name',
                               'line_col',
                               'fill_col',
                               'fill_sty',
                           ])

#______________________________________________________________________________________________________________________________


def main(tag, overwrite):
    global _xmleng
    global _xmldoc

    beams = es.BeamConfs()
    for b in beams:
        samples = es.SamplesFromBeam(b)
        for s in samples:
            filename = tag + "_" + es.ExperimentName() + "_" + util.BeamConfAsString(b) + "_" + util.SampleAsString(s) + ".xml"

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
            root_node_name = util.BeamConfAsString(b) + "_" + util.SampleAsString(s)
            root_node = make_root_node(root_node_name)

            # Read the attributes of each coarse node
            coarse_modes = es.ReactionModesCoarseFromBeamAndSample(b, s)
            for mc in coarse_modes:
                # Skip any undefined coarse modes
                if mc == python_undefined_short or mc == python_undefined or mc == -1:
                    continue

                formatting = es.PlotStyleFromBeamAndSampleAndCoarseMode(b, s, mc)

                # Set all of the possible XML fields
                fields = _param_fields(
                    name = util.CoarseModeAsString(mc),
                    line_col = formatting[0],
                    fill_col = formatting[1],
                    fill_sty = formatting[2]
                )
                print "The following fields have been set:\n", fields

                # Write all parameter fields to the output XML file
                add_param(root_node, fields)

                del fields
                print "\n\n"

            _xmleng.SaveDoc(_xmldoc, filename)
            _xmleng.FreeDoc(_xmldoc)
            print "Successfully wrote the parameter file", filename


    print "Done!"
    return


    # Check if a root node already exists. If not, create one.
def make_root_node(list_name):
    root_node = _xmleng.NewChild(0, 0, "coarse_mode_list")
    _xmleng.NewAttr(root_node, 0, "name", list_name)
    _xmleng.DocSetRootElement(_xmldoc, root_node)
    print "Creating a root node with attribute \'name\' =", list_name
    return root_node


def add_param(root_node, fields):
    maxlen_node_names = len("coarse_mode")
    maxlen_node_content = max(len(str(fields.name)), len(str(fields.line_col)), len(str(fields.fill_col)), len(str(fields.fill_sty)),)

    mc_node       = add_node(root_node, "coarse_mode", ""                  , False, True,  0                , 0)
    name_node     = add_node(mc_node  , "name"       , fields.name         , False, False, maxlen_node_names, maxlen_node_content)
    line_col_node = add_node(mc_node  , "line_col"   , str(fields.line_col), False, False, maxlen_node_names, maxlen_node_content)
    fill_col_node = add_node(mc_node  , "fill_col"   , str(fields.fill_col), False, False, maxlen_node_names, maxlen_node_content)
    fill_sty_node = add_node(mc_node  , "fill_sty"   , str(fields.fill_sty), False, False, maxlen_node_names, maxlen_node_content)
    _xmleng.AddRawLine(root_node, "\n\n")     # Additional blank lines after each coarse mode, for readability


# Improve readability and make it easier to manually insert node content (changes e.g. <node\> to <node>  <\node>)
def add_node(parent_node, node_name, node_content_in, blank_line_after = False, should_always_be_empty = False, maxlen_node = 0, maxlen_content = 0):
    node_content_print = ""
    if not should_always_be_empty:
        if maxlen_node == 0 or maxlen_content == 0:
            node_content_print = "  " + node_content_in + "  "
        else:
            node_content_print = " "*(maxlen_node-len(node_name)) + node_content_in + " "*(2+maxlen_content-len(node_content_in))
    new_node = _xmleng.NewChild(parent_node, 0, node_name, node_content_print)
    if blank_line_after:
        _xmleng.AddRawLine(parent_node, "")  # spacing between groups of nodes
    return new_node


if __name__=="__main__":

    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..",description="Create an SDK parameter XML file with hard-coded T2K-specific parameter information.")
    parser.add_option('-t', '--tag'  , dest='tag', default="coarse_mode_list", help="Tag used to name the output config file.")
    parser.add_option("-w", "--overwrite", dest="overwrite", action="store_true", help="If this option is specified, any existing file with the same output name will be deleted and overwritten.")

    o,other = parser.parse_args()
    main(o.tag, o.overwrite)

    if other:
        print "Warning - ignored arguments: ",other
