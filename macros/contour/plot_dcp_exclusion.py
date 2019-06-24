#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys
from glob import glob

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

def plot_contours_dcp(dchisq_objects, labels, pdf_filename, is_t2k_only, is_inverted=False, trues=None, bestfits=None, input_params=None, tag=None):
    canvas = ROOT.TCanvas("canvas1d", "canvas1d", 700, 500)
    canvas.SetRightMargin(0.01)
    canvas.SetBottomMargin(0.2)
    canvas.SetLeftMargin(0.12)

    #setup style options
    if len(dchisq_objects) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1
    line_colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kGray]
    line_styles = [ROOT.kSolid, ROOT.kDashed, ROOT.kDotted]

    #create & fill the histograms
    exclusions = []
    bfdiffs = []
    expercent = []
    max_val = 0.0 #get maximum of y-axis
    max_val_bf = [0.0, 0.0, 0.0] #get maximum of y-axis
    for i, dchisq_objs in enumerate(dchisq_objects):
        exclusions.append(ROOT.TH1D("hexclusions"+str(i), ";True #delta_{CP} [rad];#sigma=#sqrt{#chi^{2}_{true} - Min(#chi^{2}_{#delta_{CP}=0}, #chi^{2}_{#delta_{CP}=#pi})}", 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi()))
        bfdiffs.append([ROOT.TH1D("hbdiff"+str(i), ";True #delta_{CP} [rad];#delta_{CP}^{Best fit} - #delta_{CP}^{True} [rad]", 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi()),
                       ROOT.TH1D("hbdiffbin"+str(i), ";True #delta_{CP} [rad];#delta_{CP}^{Best fit} - #delta_{CP}^{True} [bins]", 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi()),
                       ROOT.TH1D("hbddiffchisq"+str(i), ";True #delta_{CP} [rad];#Delta#chi^{2}_{Best fit} - #Delta#chi^{2}_{True}", 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
                        ])
        assert len(bfdiffs[i]) == len(max_val_bf)
        percent = 0
        for j, d in enumerate(dchisq_objs):
            ibin = j + 1
            assert(ibin == d.FindBin(trues[i][j][0]))
            true = d.GetBinContent(ibin)
            bf   = d.GetBinContent(d.FindBin(bestfits[i][j][0]))
            pipl = d.GetBinContent(101)
            pimi = d.GetBinContent(1)
            zero = d.GetBinContent(50)
            mini = min(pipl, min(pimi, zero))
            sigma = ROOT.TMath.Sqrt(mini - true)
            print ibin, "bf %.2f true %.2f pipl %.2f pimi %.2f zero %.2f mini %.2f sigma %.2f" % (bf, true, pipl, pimi, zero, mini, sigma)
            exclusions[i].SetBinContent(ibin, 0 if ROOT.TMath.IsNaN(sigma) else sigma)
            diff = bestfits[i][j][0] - trues[i][j][0]
            bfbin = d.FindBin(bestfits[i][j][0])
            diffbin = bfbin - ibin
            print diff, "=", bestfits[i][j][0], "-", trues[i][j][0], diffbin, "bins"
            bfdiffs[i][0].SetBinContent(ibin, diff)
            bfdiffs[i][1].SetBinContent(ibin, diffbin)
            bfdiffs[i][2].SetBinContent(ibin, bf - true)
            if(sigma >= 5):
              percent += 1
        percent = 100. * percent / float(len(dchisq_objs))
        expercent.append(percent)
        print percent, "percent"

        #exclusions styles
        exclusions[i].SetLineColor(line_colors[i])
        exclusions[i].SetLineStyle(line_styles[i % len(line_styles)])
        if not i:
            exclusions[i].GetYaxis().SetTitleOffset(exclusions[i].GetYaxis().GetTitleOffset() * 0.94)
        exclusions[i].Draw("l" + ("SAME" if i else ""))

        #bfdiffs styles
        for ib, b in enumerate(bfdiffs[i]):
          b.SetLineColor(line_colors[i])
          b.SetLineStyle(line_styles[i % len(line_styles)])
          if not i:
            b.GetYaxis().SetTitleOffset(exclusions[i].GetYaxis().GetTitleOffset() * 0.94)
          thismax = max(abs(b.GetMaximum()), abs(b.GetMinimum()))
          if (ib == 1 and abs(thismax) > 10):
            thismax = 10
          elif(ib == 0 and abs(thismax) > 1):
            thismax = 1
          max_val_bf[ib] = max(thismax, max_val_bf[ib])
          
        line_color += 1
        max_val = max(max_val, exclusions[i].GetMaximum())

    #set y-axis maximum
    max_val = np.ceil(max_val) * 1.2
    exclusions[0].GetYaxis().SetRangeUser(0, max_val)
    for ib, bmax in enumerate(max_val_bf):
       bmax = np.ceil(bmax) * 1.2
       print bmax
       bfdiffs[0][ib].GetYaxis().SetRangeUser(-bmax, +bmax)

    #info box
    l = ROOT.TLatex()
    l.SetTextSize(0.05)
    #print input_params
    l.DrawLatexNDC(.01, .1, "#Deltam^{2}_{32} = %.2f #times 10^{-3}" % (input_params[0][0][conv.kPhypDm23sq] * 1E3))
    l.DrawLatexNDC(.01, .03, "#Deltam^{2}_{21} = %.2f #times 10^{-5}" % (input_params[0][0][conv.kPhypDm12sq] * 1E5))
    l.DrawLatexNDC(.29, .1, "sin^{2}#theta_{23} = %.2f" % (input_params[0][0][conv.kPhypMixing23]))
    l.DrawLatexNDC(.29, .03, "sin^{2}#theta_{12} = %.2f" % (input_params[0][0][conv.kPhypMixing12]))
    l.DrawLatexNDC(.51, .1, "sin^{2}#theta_{13} = %.7f" % (input_params[0][0][conv.kPhypMixing13]))
    l.DrawLatexNDC(.51, .03, "Hierarchy: " + input_params[0][0]["MH"])
    l.DrawLatexNDC(.15, .96, tag)

    #legend
    if len(dchisq_objects) > 0:
        lymax = 0.93
        lymin = lymax - (1 + len(dchisq_objects) + 1) * 0.03
        #in the canvas x goes from 0.12 to 0.99
        #half way is 0.435 + 0.12
        #
        xmid = 0.435 + 0.12
        xw = 0.1
        legend = ROOT.TLegend(xmid - xw, lymin, xmid + xw, lymax);
        #legend.SetFillStyle(0);
        legend.SetFillStyle(1001);
        legend.SetFillColor(10);
        legend.SetBorderSize(0);
        legend.SetTextFont(132);
        legend.SetTextSize(0.03);
        for i, exclusion in enumerate(exclusions):
            legend.AddEntry(exclusion, labels[i] + " (%.0f%%)" % (expercent[i]));

    critline = ROOT.TLine()
    critline.SetLineColor(ROOT.kBlack)
    critline.SetLineStyle(ROOT.kDashed)
    for nsigma in [3,5]:
        critline.DrawLine(-np.pi, nsigma, np.pi, nsigma)

    if len(dchisq_objects) > 1:
        legend.Draw()
    canvas.Print(pdf_filename)

    #and draw the bfdiffs
    for ib in xrange(len(max_val_bf)):
      for i, bfdiff in enumerate(bfdiffs):
        bfdiff[ib].Draw("l" + ("SAME" if i else ""))
      legend.Draw()
      critline.SetLineStyle(ROOT.kSolid)
      critline.DrawLine(-np.pi, 0, np.pi, 0)
      #print input_params
      l.DrawLatexNDC(.01, .1, "#Deltam^{2}_{32} = %.2f #times 10^{-3}" % (input_params[0][0][conv.kPhypDm23sq] * 1E3))
      l.DrawLatexNDC(.01, .03, "#Deltam^{2}_{21} = %.2f #times 10^{-5}" % (input_params[0][0][conv.kPhypDm12sq] * 1E5))
      l.DrawLatexNDC(.29, .1, "sin^{2}#theta_{23} = %.2f" % (input_params[0][0][conv.kPhypMixing23]))
      l.DrawLatexNDC(.29, .03, "sin^{2}#theta_{12} = %.2f" % (input_params[0][0][conv.kPhypMixing12]))
      l.DrawLatexNDC(.51, .1, "sin^{2}#theta_{13} = %.7f" % (input_params[0][0][conv.kPhypMixing13]))
      l.DrawLatexNDC(.51, .03, "Hierarchy: " + input_params[0][0]["MH"])
      l.DrawLatexNDC(.15, .96, tag)
      canvas.Print(pdf_filename)

    #and draw all the raw dchisq...
    debug_pdf_filename = pdf_filename.replace('.pdf', '_debug.pdf')
    canvas.Print(debug_pdf_filename + '[')
    line = ROOT.TLine()
    line.SetLineStyle(ROOT.kDashed)
    line3 = ROOT.TLine(-ROOT.TMath.Pi(), 9, +ROOT.TMath.Pi(), 9)
    line3.SetLineStyle(ROOT.kDashed)
    ymax = 25
    for j in xrange(101):
      legend.Clear()
      for i in xrange(len(dchisq_objects)):
        dchisq_objects[i][j].SetLineColor(line_colors[i])
        dchisq_objects[i][j].SetLineStyle(line_styles[i % len(line_styles)])
        if not i:
          dchisq_objects[i][j].GetYaxis().SetRangeUser(0, ymax)
        dchisq_objects[i][j].Draw("" if not i else "SAME")
        true = trues[i][j][0]
        legend.AddEntry(dchisq_objects[i][j], labels[i] + " (%.0f%%)" % (expercent[i]), "L")
      line.DrawLine(true, 0, true, ymax)
      line3.Draw()
      legend.Draw()
      l.DrawLatexNDC(.01, .1, "#Deltam^{2}_{32} = %.2f #times 10^{-3}" % (input_params[0][0][conv.kPhypDm23sq] * 1E3))
      l.DrawLatexNDC(.01, .03, "#Deltam^{2}_{21} = %.2f #times 10^{-5}" % (input_params[0][0][conv.kPhypDm12sq] * 1E5))
      l.DrawLatexNDC(.29, .1, "sin^{2}#theta_{23} = %.2f" % (input_params[0][0][conv.kPhypMixing23]))
      l.DrawLatexNDC(.29, .03, "sin^{2}#theta_{12} = %.2f" % (input_params[0][0][conv.kPhypMixing12]))
      l.DrawLatexNDC(.51, .1, "sin^{2}#theta_{13} = %.7f" % (input_params[0][0][conv.kPhypMixing13]))
      l.DrawLatexNDC(.51, .03, "Hierarchy: " + input_params[0][0]["MH"])
      l.DrawLatexNDC(.15, .96, tag)
      canvas.Print(debug_pdf_filename)
    canvas.Print(debug_pdf_filename + ']')

def get_local_histo(filename):
    in_file_tree = ROOT.TFile(filename.replace("_contour",""))
    tree = in_file_tree.Get('fitnt')
    entry = tree.GetEntry(0)
    trues = (tree.truepar0, tree.truepar1)
    bfs = (tree.bestfit0, tree.bestfit0)
    inputs = tree.input_params
    input_params = {}
    input_params[conv.kPhypMixing23] = inputs.Value(conv.kPhypMixing23)
    input_params[conv.kPhypMixing13] = inputs.Value(conv.kPhypMixing13)
    input_params[conv.kPhypMixing12] = inputs.Value(conv.kPhypMixing12)
    input_params[conv.kPhypDm23sq]   = inputs.Value(conv.kPhypDm23sq)
    input_params[conv.kPhypDm12sq]   = inputs.Value(conv.kPhypDm12sq)
    input_params["MH"] = str(tree.trueMH.c_str())
    print input_params
    in_file = ROOT.TFile(filename, "READ")
    if in_file.GetListOfKeys().Contains("h_deltachisq"):
        obj = in_file.Get("h_deltachisq")
        obj.SetDirectory(0)
        in_file.Close()
        return [obj, 1, trues, bfs, input_params]
    elif in_file.GetListOfKeys().Contains("dchisq_map_2D"):
        obj = in_file.Get("dchisq_map_2D")
        obj.SetDirectory(0)
        in_file.Close()
        return [obj, 2, trues, bfs]
    else:
        print "Error: No dchisq surface found"
        sys.exit()

def get_global_histos(filename):
    if "normal" in filename:
        filename_nh = filename
        filename_ih = filename.replace("normal", "inverted")
    else:
        filename_ih = filename
        filename_nh = filename.replace("inverted", "normal")
    nh_file = ROOT.TFile(filename_nh, "READ")
    ih_file = ROOT.TFile(filename_ih, "READ")
    if nh_file.GetListOfKeys().Contains("h_deltachisq"):
        obj_nh = nh_file.Get("h_chisq")
        obj_nh.SetDirectory(0)
        obj_ih = ih_file.Get("h_chisq")
        obj_ih.SetDirectory(0)
        nh_file.Close()
        ih_file.Close()
        minval = min(obj_nh.GetMinimum(), obj_ih.GetMinimum())
        for i in range(1, obj_nh.GetNbinsX() + 1):
            obj_nh.SetBinContent(i, obj_nh.GetBinContent(i) - minval)
            obj_ih.SetBinContent(i, obj_ih.GetBinContent(i) - minval)
        return (obj_nh, obj_ih, 1)
    elif nh_file.GetListOfKeys().Contains("dchisq_map_2D"):
        obj_nh = nh_file.Get("chisq_map_2D")
        obj_nh.SetDirectory(0)
        obj_ih = ih_file.Get("chisq_map_2D")
        obj_ih.SetDirectory(0)
        nh_file.Close()
        ih_file.Close()
        minval = min(obj_nh.GetMinimum(), obj_ih.GetMinimum())
        for i in range(1, obj_nh.GetNbinsX() + 1):
            for j in range(1, obj_nh.GetNbinsY() + 1):
                obj_nh.SetBinContent(i, j, obj_nh.GetBinContent(i, j) - minval)
                obj_ih.SetBinContent(i, j, obj_ih.GetBinContent(i, j) - minval)
        return obj_nh, obj_ih, 2
    else:
        print "Error: No dchisq surface found"
        sys.exit()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-d", "--directories", dest="in_directories", nargs="+", help="Input director(ies) for contour(s).")
    parser.add_argument("-l", "--label", dest="labels", nargs="+", help="Legend label(s) for contours(s). This should have the same number of arguments as provided to --directories.")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Legend label(s) for contours(s). This should have the same number of arguments as provided to --directories.")
    parser.add_argument("-o", "--output", dest="out_filename", help="Output PDF filename")
    parser.add_argument("--global", dest="do_global", action="store_true")
    parser.add_argument("--t2k", dest="is_t2k_only", action="store_true")
    args = parser.parse_args()
    if not (len(args.in_directories) == 1 or len(args.in_directories) == len(args.labels)):
        print "Error: The number of directories and labels should match (or can be 0 if exactly 1 directories)"
        sys.exit()

    dchisq_objects = []
    dchisq_nh_objects = []
    dchisq_ih_objects = []
    trues = []
    bfs = []
    input_params = []
    for i, directory in enumerate(args.in_directories):
        dchisq_objects.append([])
        trues.append([])
        bfs.append([])
        input_params.append([])
        for ibin in range(101):
            for f in sorted(glob(directory + '/fit_*job0_%d-_0_0_contour.root' % ibin)):
                if not args.do_global:
                    obj, dim, true, bf, inputs = get_local_histo(f)
                    dchisq_objects[i].append(obj)
                    trues[i].append(true)
                    bfs[i].append(bf)
                    input_params[i].append(inputs)
                else:
                    obj_nh, obj_ih, dim = get_global_histos(filename)
                    dchisq_nh_objects.append(obj_nh)
                    dchisq_ih_objects.append(obj_ih)
                assert(dim == 1)
                print f, "\t", trues[i][ibin][0]
        print len(dchisq_objects[i]), "contours found in directory", directory

    if not args.do_global:
        c = ROOT.TCanvas()
        c.SaveAs(args.out_filename + '[')
        plot_contours_dcp(dchisq_objects, args.labels, args.out_filename, args.is_t2k_only, trues=trues, bestfits=bfs, input_params=input_params, tag=args.tag)
        c.SaveAs(args.out_filename + ']')
    else:
        nh_filename = args.out_filename
        ih_filename = args.out_filename.replace("normal", "inverted")
        plot_contours_1d(dchisq_nh_objects, args.labels, nh_filename, args.is_t2k_only)


