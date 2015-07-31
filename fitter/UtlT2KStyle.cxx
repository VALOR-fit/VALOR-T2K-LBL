//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________

#include <fstream>

#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>

#include "UtlT2KStyle.h"
#include "Utilities.h"
#include "Utilities/Rnd.h"
#include "Enumerations.h"

using std::cout;
using std::endl;
using std::ios;
using std::ifstream;

using namespace valor;

/**************** T2K SPECIFIC *********************************************/

//____________________________________________________________________________
void valor::util::SetDefaultSKSpectrumStyle (
   valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m, TH1D * hist)
{
  if(!hist) return;

  if(s == conv::kSKsUndefined) return;
  if(m == conv::kSKmUndefined) return;

  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);

  if(sidx == -1 || midx == -1) return;

  int fill_style  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttFillStyle];
  int fill_colour = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttFillColour];
  int line_style  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineStyle];
  int line_width  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineWidth];
  int line_colour = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineColour];

  if(fill_style != -1) {
    hist->SetFillStyle(fill_style);
  }
  if(fill_colour != -1) {
    hist->SetFillColor(fill_colour);
  }
  if(line_style != -1) {
    hist->SetLineStyle(line_style);
  }
  if(line_width != -1) {
    hist->SetLineWidth(line_width);
  }
  if(line_colour != -1) {
    hist->SetLineColor(line_colour);
  }
}
//____________________________________________________________________________


// NEW 2D TEMPLATES
void valor::util::SetDefaultSKSpectrumStyle (
   valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m, TH2D * hist)
{
  if(!hist) return;

  if(s == conv::kSKsUndefined) return;
  if(m == conv::kSKmUndefined) return;

  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);

  if(sidx == -1 || midx == -1) return;

  int fill_style  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttFillStyle];
  int fill_colour = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttFillColour];
  int line_style  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineStyle];
  int line_width  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineWidth];
  int line_colour = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineColour];

  if(fill_style != -1) {
    hist->SetFillStyle(fill_style);
  }
  if(fill_colour != -1) {
    hist->SetFillColor(fill_colour);
  }
  if(line_style != -1) {
    hist->SetLineStyle(line_style);
  }
  if(line_width != -1) {
    hist->SetLineWidth(line_width);
  }
  if(line_colour != -1) {
    hist->SetLineColor(line_colour);
  }
}
//






//________________________________________________________________________________________
void valor::util::FormatErecoHist(TH1 * hist)
{
  hist->SetDirectory(0);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);
  hist->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hist->GetYaxis()->SetTitle("Events per bin");
}
//________________________________________________________________________________________

// NEW 2D TEMPLATES
void valor::util::FormatErecoVsQ2Hist(TH2 * hist)
{
  hist->SetDirectory(0);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);
  hist->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hist->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hist->GetZaxis()->SetTitle("Events per bin");
}
//________________________________________________________________________________________
