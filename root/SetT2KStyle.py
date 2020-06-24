# --- T2K style ---

import ROOT

def SetT2KStyle(WhichStyle = 1, styleName = "T2K"):
  t2kStyle= ROOT.TStyle(styleName, "T2K approved plots style")

  # -- WhichStyle --
  # 1 = presentation large fonts
  # 2 = presentation small fonts
  # 3 = publication/paper
  FontStyle = 22
  FontSizeLabel = 0.035
  FontSizeTitle = 0.05
  YOffsetTitle = 1.3

  if WhichStyle == 1:
    FontStyle = 42
    FontSizeLabel = 0.05
    FontSizeTitle = 0.065
    YOffsetTitle = 1.19
  elif WhichStyle == 2:
    FontStyle = 42
    FontSizeLabel = 0.035
    FontSizeTitle = 0.05
    YOffsetTitle = 1.6
  else :
    FontStyle = 132
    FontSizeLabel = 0.035
    FontSizeTitle = 0.05
    YOffsetTitle = 1.6

  # use plain black on white colors
  t2kStyle.SetFrameBorderMode(0)
  t2kStyle.SetCanvasBorderMode(0)
  t2kStyle.SetPadBorderMode(0)
  t2kStyle.SetPadColor(0)
  t2kStyle.SetCanvasColor(0)
  t2kStyle.SetStatColor(0)
  t2kStyle.SetFillColor(0)
  t2kStyle.SetEndErrorSize(4)
  t2kStyle.SetStripDecimals(ROOT.kFALSE)
  t2kStyle.SetLegendBorderSize(0)
  t2kStyle.SetLegendFont(FontStyle)

  # set the paper & margin sizes
  t2kStyle.SetPaperSize(20, 26)
  t2kStyle.SetPadTopMargin(0.1)
  t2kStyle.SetPadBottomMargin(0.15)
  t2kStyle.SetPadRightMargin(0.13)

  # 0.075 . 0.13 for colz option
  t2kStyle.SetPadLeftMargin(0.16)

  #to include both large/small font options
  # Fonts, sizes, offsets
  t2kStyle.SetTextFont(FontStyle)
  t2kStyle.SetTextSize(0.08)
  t2kStyle.SetLabelFont(FontStyle, "x")
  t2kStyle.SetLabelFont(FontStyle, "y")
  t2kStyle.SetLabelFont(FontStyle, "z")
  t2kStyle.SetLabelFont(FontStyle, "t")
  t2kStyle.SetLabelSize(FontSizeLabel, "x")
  t2kStyle.SetLabelSize(FontSizeLabel, "y")
  t2kStyle.SetLabelSize(FontSizeLabel, "z")
  t2kStyle.SetLabelOffset(0.015, "x")
  t2kStyle.SetLabelOffset(0.015, "y")
  t2kStyle.SetLabelOffset(0.015, "z")
  t2kStyle.SetTitleFont(FontStyle, "x")
  t2kStyle.SetTitleFont(FontStyle, "y")
  t2kStyle.SetTitleFont(FontStyle, "z")
  t2kStyle.SetTitleFont(FontStyle, "t")
  t2kStyle.SetTitleSize(FontSizeTitle, "y")
  t2kStyle.SetTitleSize(FontSizeTitle, "x")
  t2kStyle.SetTitleSize(FontSizeTitle, "z")
  t2kStyle.SetTitleOffset(1.14, "x")
  t2kStyle.SetTitleOffset(YOffsetTitle, "y")
  t2kStyle.SetTitleOffset(1.2, "z")
  t2kStyle.SetTitleStyle(0)
  t2kStyle.SetTitleFontSize(0.06)	#0.08
  t2kStyle.SetTitleFont(FontStyle, "pad")
  t2kStyle.SetTitleBorderSize(0)
  t2kStyle.SetTitleX(0.1)
  t2kStyle.SetTitleW(0.8)

  # use bold lines and markers
  t2kStyle.SetMarkerStyle(20)
  t2kStyle.SetHistLineWidth(2.5)
  t2kStyle.SetLineStyleString(2, "[12 12]")

  # postscript dashes
  # get rid of X error bars and y error bar caps
  t2kStyle.SetErrorX(0.001)

  # do not display any of the standard histogram decorations
  t2kStyle.SetOptTitle(0)
  t2kStyle.SetOptStat(0)
  t2kStyle.SetOptFit(0)

  # put tick marks on top and RHS of plots
  t2kStyle.SetPadTickX(1)
  t2kStyle.SetPadTickY(1)

  # -- color --
  t2kStyle.SetFillColor(1)

  # make color fillings (not white)
  # - color setup for 2D -
  # - "cold"/ blue-ish -
  red   = [ 0.00, 0.00, 0.00 ]
  green = [ 1.00, 0.00, 0.00 ]
  blue  = [ 1.00, 1.00, 0.25 ]

  # - "warm" red-ish colors -
  # red   = [ 1.00, 1.00, 0.25 ]
  # green = [ 1.00, 0.00, 0.00 ]
  # blue  = [ 0.00, 0.00, 0.00 ]

  stops = [ 0.25, 0.75, 1.00 ]
  NRGBs = 3
  NCont = 500
  ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
  t2kStyle.SetNumberContours(NCont)

  # - Rainbow - #
  #t2kStyle.SetPalette(1) # use the rainbow color set

  return(t2kStyle)

#____________________________________________________________________________
def CenterHistoTitles(thisHisto):
  thisHisto.GetXaxis().CenterTitle()
  thisHisto.GetYaxis().CenterTitle()
  thisHisto.GetZaxis().CenterTitle()

#____________________________________________________________________________
def AddGridLinesToPad(thisPad):
  thisPad.SetGridx()
  thisPad.SetGridy()
  return(0)
#____________________________________________________________________________
