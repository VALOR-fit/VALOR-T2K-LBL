#!/usr/bin/env python

import argparse
import ROOT
import sys
import re
from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

_nbins = 50
_text_size = 0.03

def compare(in_filenames, labels, colours, styles, tree_names, out_filename, param_map = None, apply_param_map = None, toy_limit = -1, param_names = None):

    if len(in_filenames) != len(tree_names):
        print "ERROR: input filename and tree name lists have different lengths."
        sys.exit()

    if not apply_param_map:
        apply_param_map = [False]*len(in_filenames)

    # Use the branch names from the first tree to determine what branches we'll look at in ALL trees
    first_file = ROOT.TFile(in_filenames[0], "READ")
    first_tree = first_file.Get(tree_names[0])
    branches = first_tree.GetListOfBranches()
    if not param_names:
        param_names = [branch.GetName() for branch in branches]

    print "Considering params:", param_names

    for param_name in param_names:
        print "Considering param: ", param_name
        if param_map:
            if convert_param_name(param_name, param_map) == "" or convert_param_name(param_name, param_map) == "-":
                continue

        # First get all the values for this param across all trees, so that we can get the hist ranges
        param_vals_all = []

        for itree, (filename, treename) in enumerate(zip(in_filenames, tree_names)):
            file = ROOT.TFile(filename, "READ")
            if file.GetListOfKeys().Contains(treename):
                tree = file.Get(treename)
                ntoys = toy_limit
                if toy_limit == -1:
                    ntoys = tree.GetEntries()
                if toy_limit > tree.GetEntries():
                    ntoys = tree.GetEntries()
                param_vals = [-999]*ntoys
                for itoy, toy in enumerate(tree):
                    name = param_name
                    if param_name and param_map and apply_param_map[itree]:
                        name = convert_param_name(name, param_map)
                    param_vals[itoy] = getattr(toy, name)

                    # Put special cases here, where the param value must be converted
                    if param_name == "sk_e_scale" and "legacy" in filename:
                        param_vals[itoy] += 1

                    if itoy == ntoys -1:
                        break
                param_vals_all.append(param_vals)
                del param_vals
                del tree
            else:
                print "ERROR: Could not find requested tree:", treename
                sys.exit()
            file.Close()
            del file

        # Check that the number of read param values is the same in all trees
        if len(tree_names) > 1:
            for itree in range(1, len(tree_names)):
                n_this_tree = len(param_vals_all[itree])
                n_prev_tree = len(param_vals_all[itree-1])
                if n_this_tree != n_prev_tree:
                    print "WARNING: The tree named", tree_names[itree], "has", n_this_tree, \
                        "entries, but the tree named", tree_names[itree-1], "has", n_prev_tree

        for itree in range(0, len(tree_names)):
            print "Tree #", itree, "contains", n_this_tree, "events"

        min_val = min([min(vals) for vals in param_vals_all])
        max_val = max([max(vals) for vals in param_vals_all])
        if min_val < 0:
            min_val *= 1.05
        else:
            min_val *= 0.95
        if max_val < 0:
            max_val *= 0.95
        else:
            max_val *= 1.05

        # 2 sample K-S test:
        if len(tree_names) == 2:
            gof = ROOT.Math.GoFTest(len(param_vals_all[0]), array('d', param_vals_all[0]),
                                    len(param_vals_all[1]), array('d', param_vals_all[1]))
            pval = gof.KolmogorovSmirnov2SamplesTest()
            print "The toy throw distributions for param", param_name, \
                " have a K-S test p-value of", pval
            del gof

        # Fill the hists for this param
        hists_to_draw = []
        for i in range(0, len(tree_names)):
            hname = param_name + "_" + str(i)
            hist = ROOT.TH1D(hname, hname, _nbins, min_val, max_val)
            hist.SetDirectory(0)
            for iv, val in enumerate(param_vals_all[i]):
                hist.Fill(val)
            hist.Scale(1./hist.Integral())
            hists_to_draw.append(hist)

        hist_max = max([h.GetMaximum() for h in hists_to_draw])

        # Draw the distribution of throws for this param
        canv_name = "canvas_" + param_name
        canvas = ROOT.TCanvas(canv_name, canv_name, 700, 500)
        first = True
        for i, hist in enumerate(hists_to_draw):
            draw_opt = "HIST"
            if first:
                first = False
                hist.GetXaxis().SetTitle(param_name)
                hist.GetYaxis().SetTitle("Entries")
                hist.SetMaximum(1.1 * hist_max)
            else:
                draw_opt += ",SAME"
            hist.SetLineColor(colours[i])
            hist.SetLineStyle(styles[i])
            hist.Draw(draw_opt)

        # Draw the legend
        lymax = 0.89
        lymin = lymax - (1 + len(hists_to_draw)) * _text_size
        lxmin = 0.15
        lxmax = 0.3
        legend = ROOT.TLegend(lxmin, lymin, lxmax, lymax)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(132)
        legend.SetTextSize(_text_size)
        if len(labels) > 0:
            for i, hist in enumerate(hists_to_draw):
                legend.AddEntry(hist, labels[i], "l")
        legend.Draw()

        # Save to the output file
        save_name = out_filename
        if len(param_names) > 1:
            if param_name == param_names[0]:
                save_name += "("
            elif param_name == param_names[len(param_names) -1]:
                save_name += ")"
        canvas.SaveAs(save_name, "RECREATE")

        # Free up memory
        del hists_to_draw[:]
        del param_vals_all[:]
        del legend
        del canvas


def convert_param_name(param_name, arg_file_name):
    if arg_file_name == "":
        return param_name
    parameters = list_of_tuples(arg_file_name)
    comp_param_name = ""
    for iparam,(name1,name2) in enumerate(parameters):
        if name1 == param_name:
            comp_param_name = name2
            break
    return comp_param_name


def list_of_tuples(filename):
    tuplelist = []
    with open(filename) as file:
        for line in file:
            filtered = re.sub('[\"\'\(]|(\)\s*\,)|[\)]|(#.*$)|^\s*$|\s', '', line)
            if re.match(r'^\s*$', filtered):
                continue
            tuplelist.append(tuple(filtered.split(',')))

    return tuplelist


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename"       , dest="in_filenames" , nargs="+"            , help="Input filename(s) for contour(s).")
    parser.add_argument("-l", "--label"          , dest="labels"       , nargs="+"            , help="Legend label(s) for contours(s). This should have either 0 length or the same number of arguments as provided to --filename.")
    parser.add_argument("-c", "--colour"         , dest="colours"      , nargs="+"            , help="Custom line colour(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename.")
    parser.add_argument("-s", "--style"          , dest="styles"       , nargs="+"            , help="Custom line style(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename.")
    parser.add_argument("-t", "--tree-name"      , dest="tree_names"   , nargs="+"            , help="Name(s) of the trees in the input file(s). If used, this should have the same number of arguments as provided to --filename.")
    parser.add_argument("-o", "--output"         , dest="out_filename"                        , help="Output PDF filename")
    parser.add_argument("-m", "--param-map"      , dest="param_map"                           , help="File name of the parameter name map, used to convert the param names read from the first file into parm names for the other files.")
    parser.add_argument("-a", "--apply-param-map", dest="apply_map_str", nargs="+"            , help="1 to apply the parameter map to this file, 0 to not.")
    parser.add_argument(      "--toy-limit"      , dest="toy_limit"    , type=int, default=-1 , help="Only consider the first X toys in each tree. set to -1 (default) to look at all toys.")
    parser.add_argument(      "--params"         , dest="params"       , nargs="+"            , help="A list of parameters to plot. If this is not specified, all params are plot.")
    args = parser.parse_args()
    labels = []
    if args.labels:
        labels = [i for i in args.labels]
        if not len(args.in_filenames) == len(labels):
            print "Error: The number of filenames and labels should match"
            sys.exit()

    colours = []
    if args.colours:
        colours = [int(i) for i in args.colours]
    styles = []
    if args.styles:
        styles = [int(i) for i in args.styles]
    apply_map = [False]*len(args.in_filenames)
    if args.apply_map_str:
        apply_map = [bool(int(i)) for i in args.apply_map_str]
        if not len(args.in_filenames) == len(apply_map):
            print "Error: The number of --filename and --apply-param-map entries should match"
            sys.exit()

    compare(args.in_filenames, labels, colours, styles, args.tree_names, args.out_filename, args.param_map, apply_map, args.toy_limit, args.params)


