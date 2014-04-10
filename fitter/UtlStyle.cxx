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

#include "UtlStyle.h"
#include "Utilities.h"
#include "Opt.h"
#include "Rnd.h"
#include "Enumerations.h"

using std::cout;
using std::endl;
using std::ios;
using std::ifstream;

using namespace t2kana;

//___________________________________________________________________________
void util::SetDefaultStyle(bool black_n_white)
{
  gROOT->SetStyle("Plain");

  gStyle -> SetPadTickX (1);
  gStyle -> SetPadTickY (1);

  //
  // Turn off all borders
  //
  gStyle -> SetCanvasBorderMode (0);
  gStyle -> SetFrameBorderMode  (0);
  gStyle -> SetPadBorderMode    (0);
  gStyle -> SetDrawBorder       (0);
  gStyle -> SetCanvasBorderSize (0);
  gStyle -> SetFrameBorderSize  (0);
  gStyle -> SetPadBorderSize    (0);
  gStyle -> SetTitleBorderSize  (0);

  //
  // Set the size of the default canvas
  //
  gStyle -> SetCanvasDefH (600);
  gStyle -> SetCanvasDefW (730);
  gStyle -> SetCanvasDefX  (10);
  gStyle -> SetCanvasDefY  (10);

  //
  // Set marker style
  //
  gStyle -> SetMarkerStyle (20);
  gStyle -> SetMarkerSize   (1);

  //
  // Set line widths
  //
  gStyle -> SetFrameLineWidth (1);
  gStyle -> SetFuncWidth      (2);
  gStyle -> SetHistLineWidth  (3);
  gStyle -> SetFuncColor      (2);
  gStyle -> SetFuncWidth      (3);

  //
  // Set margins
  //
  gStyle -> SetPadTopMargin    (0.10);
  gStyle -> SetPadBottomMargin (0.20);
  gStyle -> SetPadLeftMargin   (0.15);
  gStyle -> SetPadRightMargin  (0.03);

  // 
  // Set tick marks and turn off grids
  //
  gStyle -> SetNdivisions (505,"xyz");

  //
  // Adjust size and placement of axis labels
  //
  gStyle -> SetLabelSize   (0.050,  "xyz");
  gStyle -> SetLabelOffset (0.005,  "x"  );
  gStyle -> SetLabelOffset (0.005,  "y"  );
  gStyle -> SetLabelOffset (0.005,  "z"  );
  gStyle -> SetTitleSize   (0.060,  "xyz");
  gStyle -> SetTitleOffset (1.200,  "xz" );
  gStyle -> SetTitleOffset (1.000,  "y"  );

  //
  // Adjust fonts
  //
  gStyle -> SetLabelFont(42, "xyz");
  gStyle -> SetTitleFont(42, "xyz");
  gStyle -> SetStatFont (42);

  //
  // Set Data/Stat/... and other options
  //
  gStyle -> SetOptDate          (0);
  gStyle -> SetOptFile          (0);
  gStyle -> SetOptStat          (0);
  gStyle -> SetStatFormat       ("6.2f");
  gStyle -> SetFitFormat        ("8.4f");
  gStyle -> SetOptFit           (1);
  gStyle -> SetStatH            (0.20);
  gStyle -> SetStatStyle        (0);
  gStyle -> SetStatW            (0.30);
  gStyle -> SetStatX            (0.845);
  gStyle -> SetStatY            (0.845);
  gStyle -> SetOptTitle         (0);
  gStyle -> SetTitleX           (0.15);
  gStyle -> SetTitleW           (0.75);
  gStyle -> SetTitleY           (0.90);
  gStyle -> SetPalette          (1);
  gStyle -> SetLegendBorderSize (0);

  // 
  // Set paper size for life in the US or EU
  //
  gStyle -> SetPaperSize (TStyle::kA4);       //<-- tartes aux fraises
//gStyle -> SetPaperSize (TStyle::kUSLetter); //<-- donuts

  //
  // In B&W (papers)
  //
  if(black_n_white){
    const int ncol = 7;

    double red   [ncol];
    double green [ncol];
    double blue  [ncol];
    double stops [ncol];

    double dcol = -1/double(ncol);
    double gray = 1;
    for (int j = 0; j < ncol; j++) {
      // Define color with RGB equal to : gray, gray, gray
      stops[j] = double(j)/double(ncol-1);
      red  [j] = gray;
      blue [j] = gray;
      green[j] = gray;

      gray += dcol;
    }
    UInt_t totcol=50;
    TColor::CreateGradientColorTable(ncol,stops,red,green,blue,totcol);

    gStyle -> SetFuncWidth     (1);
    gStyle -> SetHistLineWidth (1);
    gStyle -> SetFuncColor     (1);
    gStyle -> SetFuncWidth     (1);
  }//bw

  gROOT->ForceStyle();
}
//___________________________________________________________________________
void util::Format(
   TGraph* gr, int lcol, int lsty, int lwid, int mcol, int msty, double msiz)
{
  if(!gr) return;

  if (lcol >= 0) gr -> SetLineColor   (lcol);
  if (lsty >= 0) gr -> SetLineStyle   (lsty);
  if (lwid >= 0) gr -> SetLineWidth   (lwid);


  if (mcol >= 0) gr -> SetMarkerColor (mcol);
  if (msty >= 0) gr -> SetMarkerStyle (msty);
  if (msiz >= 0) gr -> SetMarkerSize  (msiz);
}
//___________________________________________________________________________
void util::Format(
  TH1* hst, int lcol, int lsty, int lwid, int mcol, int msty, double msiz)
{
  if(!hst) return;

  if (lcol >= 0) hst -> SetLineColor   (lcol);
  if (lsty >= 0) hst -> SetLineStyle   (lsty);
  if (lwid >= 0) hst -> SetLineWidth   (lwid);

  if (mcol >= 0) hst -> SetMarkerColor (mcol);
  if (msty >= 0) hst -> SetMarkerStyle (msty);
  if (msiz >= 0) hst -> SetMarkerSize  (msiz);
}
//___________________________________________________________________________
void util::PrintBanner(void)
{
  string base_dir = string(gSystem->Getenv("VALOR_LBNX"));
  string banner = base_dir + string("/data/banner/banner.dat");
  util::PrintBanner(banner, 500);
}
//___________________________________________________________________________
void util::PrintBanner(string filename, UInt_t wait_msec)
{
  ifstream banner(filename.c_str(), ios::in);

  if( banner.is_open() ) {
      banner.seekg(0, ios::end);

      int    length = banner.tellg();
      char * buffer = new char[length];

      banner.seekg(0, ios::beg);
      banner.read(buffer, length);

      //cout << "\n\n" << buffer << "\n" << endl;
      cout << "\n\n";
      cout.write(buffer,length);
      cout << "\n" << endl;

      delete [] buffer;

      gSystem->Sleep(wait_msec); // watch the banner for a while
  }
}
//___________________________________________________________________________


/**************** T2K SPECIFIC *********************************************/

//____________________________________________________________________________
void t2kana::util::SetDefaultSKSpectrumStyle (
   t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, TH1D * hist)
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
//________________________________________________________________________________________
void t2kana::util::FormatErecoHist(TH1 * hist)
{
  hist->SetDirectory(0);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);
  hist->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hist->GetYaxis()->SetTitle("Events per bin");
}
//________________________________________________________________________________________
