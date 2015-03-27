#!/usr/bin/python
import ROOT

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


def main(file_patterns,do_systematics=False):
  ch = ROOT.TChain("fitnt")
  ch.Add(file_patterns[0])
  assert(ch.GetEntries())
  ch.GetEntry(0)
  
  fparm0 = ch.best_fit_params
  
  params_to_plot = []
  for p in xrange(conv.knppmax):
    if (not fparm0.Exists(p)): continue
    if (fparm0.IsFixed(p)): continue
    
    if util.IsSystematic(p) and not do_systematics: continue
    params_to_plot.append(p)
    assert(ch.input_params.Contains(p))
  
  pull_hists = dict()
  
  for ipattern,file_pattern in enumerate(file_patterns):
    pull_hists[ipattern] = dict()
    for p in params_to_plot:
      hname = "pulls_f"+`ipattern`+"_"+util.ParamAsStringSimple(p)
      h = ROOT.TH1D(hname,hname,50,-5,5)
      h.GetXaxis().SetTitle("Pull")
      pull_hists[ipattern][p] = h
      
    ch = ROOT.TChain("fitnt")
    ch.Add(file_patterns[ipattern])
      
    for entry in ch:
      fparm = entry.best_fit_params
      in_params = entry.input_params
      for p in params_to_plot:
        bf_val = fparm.BestFitValue(p)
        in_val = in_params.Value(p)
        hesse_error = fparm.HesseError(p)
        pull = (bf_val-in_val)/hesse_error
        pull_hists[ipattern][p].Fill(pull)
  
  
  for p in params_to_plot:
    cname = "canv_"+util.ParamAsStringSimple(p)
    canv = ROOT.TCanvas(cname,cname,1280,720)
    canv.cd()
    plot_opt = ""
    histmax = 1.1*max([h[p].GetMaximum() for h in pull_hists.itervalues()])
    for ipattern,pattern in enumerate(file_patterns):
      h = pull_hists[ipattern][p]
      h.SetLineColor(colors[ipattern])
      h.Draw(plot_opt)
      plot_opt = "same"
      
    canv.Print("pull_"+util.ParamAsStringSimple(p)+".eps")

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] pattern1 [pattern2] ...",description="Make plots of physics parameter pulls")
  parser.add_option('-s','--do-systematics', dest='do_systematics'    , action="store_true", default=False, help="Draw the systematic parameters as well")
  
  options,patterns = parser.parse_args()
  
  main(patterns,options.do_systematics)
