#!/usr/bin/env python
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

to3dp = "{0:.3f}"
orig_plot_opt = ""

text_entry_height = 0.08

npads_x = 2
npads_y = 5

def main(file_patterns,labels=[None],do_systematics=False,print_stats=False):
  
  objects_to_keep = []
  
  ch = ROOT.TChain("fitnt")
  ch.Add(file_patterns[0])
  assert(ch.GetEntries())
  ch.GetEntry(0)
  
  fparm0 = ch.best_fit_params
  
  while len(labels)<len(file_patterns):
    labels.append("")
    
  draw_legend = any(labels)
  ymax = 0.9
  ymin = ymax
  for l in labels:
    if l: ymin -= text_entry_height
  legend = ROOT.TLegend(0.7,ymin,0.9,0.9)
  
  stat_boxes = []
  ymin = ymax
  for ipattern,pattern in enumerate(file_patterns):
    this_ymax = ymin
    ymin -= text_entry_height
    box = ROOT.TPaveText(0.1,ymin,0.4,this_ymax,"ndc")
    box.SetTextColor(colors[ipattern])
    box.SetShadowColor(0)
    stat_boxes.append(box)
  
  params_to_plot = []
  for p in xrange(conv.knppmax):
    if (not fparm0.Exists(p)): continue
    if (fparm0.IsFixed(p)): continue
    
    params_to_plot.append(p)
    assert(ch.input_params.Contains(p))
  
  pull_hists = dict()
  
  for ipattern,file_pattern in enumerate(file_patterns):
    pull_hists[ipattern] = dict()
    for p in params_to_plot:
      hname = "pulls_f"+`ipattern`+"_"+util.ParamAsStringSimple(p)
      h = ROOT.TH1D(hname,hname,50,-5,5)
      h.GetXaxis().SetTitle("Pull")
      h.GetYaxis().SetTitle("Fraction of fits")
      h.SetTitle(util.ParamAsStringPaperRooTex(p).replace(";"," - "))
      pull_hists[ipattern][p] = h
      
    ch = ROOT.TChain("fitnt")
    ch.Add(file_patterns[ipattern])
    
    entries = ch.GetEntries()
      
    for entry in ch:
      fparm = entry.best_fit_params
      in_params = entry.input_params
      for p in params_to_plot:
        bf_val = fparm.BestFitValue(p)
        in_val = in_params.Value(p)
        hesse_error = fparm.HesseError(p)
        pull = (bf_val-in_val)/hesse_error
        pull_hists[ipattern][p].Fill(pull)
  
    for h in pull_hists[ipattern].itervalues():
      h.Scale(1./entries)
  
  big_canvs = [ROOT.TCanvas("big_canv_0","big_canv_0",2100,2970)]
  
  curr_big_canv = big_canvs[0]
  curr_big_canv.Divide(npads_x,npads_y)
  ipad = 1
  npads = npads_x*npads_y
  
  for p in params_to_plot:
    legend.Clear()
    for box in stat_boxes: box.Clear()
    
    canv = curr_big_canv.cd(ipad)
    
    plot_opt = orig_plot_opt
    histmax = 1.1*max([h[p].GetMaximum() for h in pull_hists.itervalues()])
    for ipattern,pattern in enumerate(file_patterns):
      h = pull_hists[ipattern][p]
      h.SetMaximum(histmax)
      h.SetLineColor(colors[ipattern])
      h.Draw(plot_opt)
      
      legend.AddEntry(h,labels[ipattern],"l")
      stat_boxes[ipattern].AddText("mean="+to3dp.format(h.GetMean())+", RMS="+to3dp.format(h.GetRMS()))
      
      plot_opt = orig_plot_opt+"same"
    
    if draw_legend: legend.Draw()
    if print_stats:
      for box in stat_boxes:
        clone_box = box.Clone()
        clone_box.Draw()
        objects_to_keep.append(clone_box)
    
    if (do_systematics or not util.IsSystematic(p)):
      canv.Print("pull_"+util.ParamAsStringSimple(p)+".eps")
    
    title = ROOT.TPaveText(0.,0.93,1.,1.,"ndc")
    title.SetShadowColor(0)
    title.SetFillColor(0)
    title.SetFillStyle(0)
    title.SetLineColor(0)
    title.AddText(util.ParamAsStringPaperRooTex(p))
    objects_to_keep.append(title)
    title.Draw()
    
    ipad+=1
    
    if ipad>npads:
      bcname = "big_canv_"+`len(big_canvs)`
      curr_big_canv = ROOT.TCanvas(bcname,bcname,2100,2970)
      curr_big_canv.Divide(npads_x,npads_y)
      ipad = 1
      big_canvs.append(curr_big_canv)
    
  for icanv,canv in enumerate(big_canvs):
    canv.Print("many_pulls_"+`icanv`+".eps")

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] pattern1 [pattern2] ...",description="Make plots of physics parameter pulls")
  parser.add_option('-s','--do-systematics', dest='do_systematics'    , action="store_true", default=False, help="Draw the systematic parameters individually as well")
  parser.add_option('-l','--labels', dest='labels', default="", help="Legend entries, comma seperated")
  parser.add_option('-p','--print_stats', dest='do_stats', action="store_true", default=False, help="Draw the Gaussian fit parameters as well")
  
  options,patterns = parser.parse_args()
  labels = options.labels.split(",")
  
  main(patterns,labels,options.do_systematics,options.do_stats)
