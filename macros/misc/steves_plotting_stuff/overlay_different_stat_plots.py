#!/usr/bin/python
import ROOT

R = ROOT

R.gROOT.SetBatch(1)

old_POT = "2.315 #times 10^{20} POT"
new_POT = "4.011 #times 10^{20} POT"

f_old = ROOT.TFile("valor_fit_contours_old.root")
f_new = ROOT.TFile("valor_fit_contours.root")

confidences = [("90",R.kSolid),
               ("68",R.kDashed)]

hierarchies = ["nh","ih"]

title_size = 0.08
label_size = 0.07

def main():
  for hierarchy in hierarchies:
    bf_old = f_old.Get("valor_"+hierarchy+"_bestfit")
    bf_new = f_new.Get("valor_"+hierarchy+"_bestfit")
    
    contours_old = []
    contours_new = []
    for conf,style in confidences:
      cont_old = f_old.Get("valor_"+hierarchy+"_"+conf+"_contour")
      cont_new = f_new.Get("valor_"+hierarchy+"_"+conf+"_contour")
      
      cont_old.SetLineStyle(style)
      cont_new.SetLineStyle(style)
      
      contours_old.append(cont_old)
      contours_new.append(cont_new)
    
    cname = "c_"+hierarchy
    canv = ROOT.TCanvas(cname,cname,1280,720)
    canv.SetRightMargin(0.04)
    canv.SetTopMargin(0.04)
    canv.SetBottomMargin(0.16)
    canv.SetLeftMargin(0.17)
    canv.cd()
    frame = canv.DrawFrame(0.2,1.3E-3,0.8,4.1E-3)
    frame.GetXaxis().SetTitle("sin^{2} #bar{#theta_{23}}")
    if hierarchy == "nh":
      frame.GetYaxis().SetTitle("|#bar{#Delta m^{2}_{32}}|")
    elif hierarchy=="ih":
      frame.GetYaxis().SetTitle("|#bar{#Delta m^{2}_{31}}|")
    else:
      frame.GetYaxis().SetTitle("#bar{#Delta m^{2}}")
    
    frame.GetXaxis().SetTitleOffset(0.8)
    frame.GetYaxis().SetTitleOffset(1.07)
    frame.GetXaxis().SetTitleSize(title_size)
    frame.GetYaxis().SetTitleSize(title_size)
    
    #~ frame.GetXaxis().SetLabelOffset(0.05)
    #~ frame.GetYaxis().SetLabelOffset(label_size)
    frame.GetXaxis().SetLabelSize(label_size)
    frame.GetYaxis().SetLabelSize(label_size)
    
    for cont in contours_old:
      cont.SetLineColor(R.kRed)
      cont.Draw("l same")
      
    for cont in contours_new:
      cont.SetLineColor(R.kBlack)
      cont.Draw("l same")
      
    bf_old.SetMarkerColor(contours_old[0].GetLineColor())
    bf_new.SetMarkerColor(contours_new[0].GetLineColor())
    
    for bf in [bf_old,bf_new]:
      bf.SetMarkerSize(3)
      bf.SetMarkerStyle(R.kFullStar)
      bf.Draw("p same")
      
    CopyMarkerToHist(bf_old,contours_old[0])
    CopyMarkerToHist(bf_new,contours_new[0])
      
    legend = ROOT.TLegend(0.6,0.7,0.96,0.96)
    legend.AddEntry(contours_old[0],old_POT)
    legend.AddEntry(contours_new[0],new_POT)
    
    #~ legend.SetBorderSize(0)
    legend.SetFillStyle (0)
    legend.Draw()
    canv.Update()
    
    canv.Print("comparison_of_different_POT_"+hierarchy+".eps")
  
  
def CopyMarkerToHist(marker,histogram):
  histogram.SetMarkerSize (marker.GetMarkerSize ())
  histogram.SetMarkerStyle(marker.GetMarkerStyle())
  histogram.SetMarkerColor(marker.GetMarkerColor())

if __name__=="__main__":
  main()
