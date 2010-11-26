//____________________________________________________________________________
/*

*/
//____________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TXMLEngine.h>
#include <TFile.h>

#include "year1_numu_osc_analysis/Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setfill;
using std::setprecision;
using std::ofstream;

using namespace t2kana;
using namespace t2kana::yr1numu::conv;

TRandom3 gR(89289281);

//____________________________________________________________________________
double t2kana::yr1numu::util::DeltaChisquare(
  ConfLevel_t cl, int nfitparams) 
{
  switch(cl) {
      case (kCLUndefined) : 
      {
	return -99999;  
        break;
      }
      case (kCL68) :
      {
        if      (nfitparams == 1)  { return  1.00;  }
        else if (nfitparams == 2)  { return  2.30;  }
        else if (nfitparams == 3)  { return  3.53;  }
        else                       { return -99999; }
        break;
      }
      case (kCL90) :
      {
        if      (nfitparams == 1)  { return  2.71;  }
        else if (nfitparams == 2)  { return  4.61;  }
        else if (nfitparams == 3)  { return  6.25;  }
        else                       { return -99999; }
        break;
      }
      case (kCL95) :
      {
        if      (nfitparams == 1)  { return  3.84;  }
        else if (nfitparams == 2)  { return  5.99;  }
        else if (nfitparams == 3)  { return  7.82;  }
        else                       { return -99999; }
        break;
      }
      case (kCL99) :
      {
        if      (nfitparams == 1)  { return  6.63;  }
        else if (nfitparams == 2)  { return  9.21;  }
        else if (nfitparams == 3)  { return 11.34;  }
        else                       { return -99999; }
        break;
      }
      default:
        return -99999;     
    }
  return -99999;
}
//____________________________________________________________________________
double t2kana::yr1numu::util::ChisqGauss(
  double ndata, double nmc)
{
  double chisq = 0.;
  if(ndata>0) {
    chisq = TMath::Power(ndata-nmc,2.)/ndata;
  }
  return chisq;
}
//____________________________________________________________________________
double t2kana::yr1numu::util::ChisqPoisson(
  double ndata, double nmc)
{
  const double epsilon = kASmallNum;

  if(TMath::Abs(nmc-ndata) < epsilon) {
    return 0.;
  }

  if(nmc   <= 0.) { nmc   = epsilon; }
  if(ndata <= 0.) { ndata = epsilon; }

  double chisq = 2. * TMath::Abs(nmc-ndata) + 
                 2. * ndata * TMath::Log(ndata/nmc);

  return chisq;
}
//____________________________________________________________________________
double t2kana::yr1numu::util::NuMuSurvProb2D(
  double sinsq2theta, double dmsq, double Et)
{
  if(Et <= 0.) return 0.;

  double osc   = TMath::Sin(1.267*dmsq*kBaseline/Et);
  double posc  = sinsq2theta * TMath::Power(osc,2);
  double psurv = 1-posc;

  psurv = TMath::Min(1.,psurv);
  psurv = TMath::Max(0.,psurv);

  return psurv;
}
//____________________________________________________________________________
double t2kana::yr1numu::util::NuMuSurvProb2D(
  double sinsq2theta, double dmsq, double Et, 
  t2kana::yr1numu::conv::SuperkMode_t m)
{
  if(m == kSKmNC) {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE || m == kSKmNuMuNonCCQE) {
    return NuMuSurvProb2D(sinsq2theta, dmsq, Et);
  }
  else {
    cerr << "NuMuSurvProb2D> Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
TMatrixD t2kana::yr1numu::util::CholeskyDecomposition(TMatrixD& cov_matrix)
{
// Perform a Cholesky decomposition of the input covariance matrix and
// return the lower triangular matrix
//
   int ncols = cov_matrix.GetNcols();
   int nrows = cov_matrix.GetNrows();

   assert(ncols==nrows);

   int n = nrows;

   TMatrixD L(n, n);

   for (int i = 0; i < n; ++i) {

     // calculate the diagonal term first
     L(i,i) = cov_matrix(i,i);
     for (int k = 0; k < i; ++k) {
        double tmp = L(k,i);
        L(i,i) -= tmp*tmp;
     }//k
     if(L(i,i) < 0) {
       cerr << "Decomposed covariance matrix not positive-definite" << endl;
       exit(1);
     }
     L(i,i) = TMath::Sqrt(L(i,i));
     // then the off-diagonal terms
     for (int j = i+1; j < n; ++j) {
        L(i,j) = cov_matrix(i,j);
        for (int k = 0; k < i; ++k) {
           L(i,j) -= L(k,i)*L(k,j);
        }
        L(i,j) /= L(i,i);
     }//j
  }//i

  // create the transpose of L
  TMatrixD LT(TMatrixD::kTransposed,L);

#ifdef __debug_level_1___
  cout << "Cholesky matrix L^{T}" << endl;
  for(int ip=0; ip<n; ip++) {
    for(int jp=0; jp<n; jp++) {
      cout << setw(15) << setfill(' ') << setprecision(6) << LT[ip][jp] << " | ";
    }//jp
    cout << endl;
  }//ip
#endif

  return LT;
}
//____________________________________________________________________________
TVectorD t2kana::yr1numu::util::CholeskyGenerateCorrelatedParams (
    TMatrixD& cholesky_triangular, TVectorD& mean_params) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();

  assert(ncols==nrows);
  assert(npars==nrows);

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = gR.Gaus();
  }
  g *= cholesky_triangular;

#ifdef __debug_level_3__
  ofstream foutp1("cholesky_uncorrelated_params.txt",fstream::app);
  for (int i = 0; i < n; ++i) {
    foutp1 << g[i] << "\t";
  }
  foutp1 << endl;
  foutp1.close();
#endif

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g(i);
     correlated_params[i] = v;
  }

#ifdef __debug_level_3__
  ofstream foutp2("cholesky_correlated_params.txt",fstream::app);
  for (int i = 0; i < n; ++i) {
    foutp2 << correlated_params[i] << "\t";
  }
  foutp2 << endl;
  foutp2.close();
#endif

  return correlated_params;
}
//____________________________________________________________________________
TVectorD t2kana::yr1numu::util::CholeskyGenerateCorrelatedParamVariations (
    TMatrixD& cholesky_triangular) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();

  assert(ncols==nrows);

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = gR.Gaus();
  }
  g *= cholesky_triangular;

  return g;
}
//____________________________________________________________________________
string t2kana::yr1numu::util::AsString(PhysParamId_t p)
{
  switch(p) {
    case ( kPhypSinsq2Theta23 ) : return "sinsq2theta23"; break;
    case ( kPhypDm23sq        ) : return "dm23sq";        break;
    default                     : return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::yr1numu::util::AsString(
    t2kana::yr1numu::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "FCFV1RingMuLike"; break;
     case ( kSKsFCFVmRingMuLike ) : return "FCFVmRingMuLike"; break;
     default                      : return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::yr1numu::util::AsString(
     t2kana::yr1numu::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE    ) : return "CCQE";    break;
     case ( kSKmNuMuNonCCQE ) : return "CCnonQE"; break;
     case ( kSKmNC          ) : return "NC";      break;
     default                  : return "--";
  }
  return "--";
}
//____________________________________________________________________________
void t2kana::yr1numu::util::SetDefaultStyle(bool black_n_white)
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
void t2kana::yr1numu::util::Format(
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
void t2kana::yr1numu::util::Format(
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
string t2kana::yr1numu::util::TrimSpaces(string input)
{
  if( input.find_first_not_of(" \n") != 0)
    input.erase( 0,  input.find_first_not_of(" \n")  );

  if( input.find_last_not_of(" \n") != input.length() )
    input.erase( input.find_last_not_of(" \n")+1, input.length() );

  return input;
}
//____________________________________________________________________________
vector<string> t2kana::yr1numu::util::Tokens (string input, string delimiter)
{
// split a string of 'delimiter' separated values and return each string
// part as a vector<string> element.

  vector<string> tokens;

  while(input.find_first_of(delimiter) < input.length()) {

    tokens.push_back( input.substr(0, input.find_first_of(delimiter)) );
    input = input.substr(input.find_first_of(delimiter)+1, input.length());
  }
  tokens.push_back(input);
  return tokens;
}
//___________________________________________________________________________
int t2kana::yr1numu::util::Sign(double tweaking_dial)
{
  if(TMath::Abs(tweaking_dial) < kASmallNum) return 0;

  if (tweaking_dial < 0) return -1;
  return +1;
}
//___________________________________________________________________________
XMLNodePointer_t t2kana::yr1numu::util::FindNode(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Find & return the XMLNodePointer_t at the specified path
//
  vector<string> node_path = Tokens(path,"/");

  unsigned int ndepth = node_path.size();
  unsigned int idepth = 0;

  XMLNodePointer_t curr_node = xmleng->GetChild(topnode);
  while(curr_node) {
    string curr_node_name = xmleng->GetNodeName(curr_node);
    if(node_path[idepth] == curr_node_name) {
        idepth++;
        if(idepth == ndepth) {
          return curr_node;
        } else {
          curr_node = xmleng->GetChild(curr_node);
        }
    } else {
       curr_node = xmleng->GetNext(curr_node);
    }
  }
  return 0;
}
//___________________________________________________________________________
TH1F * t2kana::yr1numu::util::XML2TH1F(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH1F instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TH1F >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH1F >> No node was found at input XML node path: "
         << path << endl;
    return 0;
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TH1F >> node = " << child_node_name << endl;
    cout << "utils/XML2TH1F >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TH1F >> filename : " << filename << endl;
  cout << "utils/XML2TH1F >> objname  : " << objname  << endl;
#endif

  TFile file(filename.c_str(), "read");
  TH1F * h = (TH1F*)file.Get(objname.c_str());
  if(!h) {
    cerr << "utils/XML2TH1F >> No " << objname << " in " << filename << endl;
    file.Close();
    return 0;
  }

  TH1F * hcopy = new TH1F(*h);
  hcopy->SetDirectory(0);

  file.Close();

  return hcopy;
}
//___________________________________________________________________________
TH1D * t2kana::yr1numu::util::XML2TH1D(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH1D instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TH1D >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH1D >> No node was found at input XML node path: "
         << path << endl;
    return 0;
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TH1D >> node = " << child_node_name << endl;
    cout << "utils/XML2TH1D >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TH1D >> filename : " << filename << endl;
  cout << "utils/XML2TH1D >> objname  : " << objname  << endl;
#endif

  TFile file(filename.c_str(), "read");
  TH1D * h = (TH1D*)file.Get(objname.c_str());
  if(!h) {
    cerr << "utils/XML2TH1D >> No " << objname << " in " << filename << endl;
    file.Close();
    return 0;
  }

  TH1D * hcopy = new TH1D(*h);
  hcopy->SetDirectory(0);

  file.Close();

  return hcopy;
}
//___________________________________________________________________________
TH2D * t2kana::yr1numu::util::XML2TH2D(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH2D instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TH2D >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH2D >> No node was found at input XML node path: "
         << path << endl;
    return 0;
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TH2D >> node = " << child_node_name << endl;
    cout << "utils/XML2TH2D >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TH1D >> filename : " << filename << endl;
  cout << "utils/XML2TH1D >> objname  : " << objname  << endl;
#endif

  TFile file(filename.c_str(), "read");
  TH2D * h2 = (TH2D*)file.Get(objname.c_str());
  if(!h2) {
    cerr << "utils/XML2TH2D >> No " << objname << " in " << filename << endl;
    file.Close();
    return 0;
  }

  TH2D * h2copy = new TH2D(*h2);
  h2copy->SetDirectory(0);

  file.Close();

  return h2copy;
}
//___________________________________________________________________________
TAxis *  t2kana::yr1numu::util::XML2BinInfo(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
//
//

#ifdef __debug_level_2__
  cout << "utils/XML2BinInfo >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH2D >> No node was found at input XML node path: "
         << path << endl;
    return 0;
  }

  string binsz = xmleng->GetAttr(found_node,"bin_size");
  if(binsz == "fixed") {
    int    nbins = 0;
    double min   = 0;
    double max   = 0;
    XMLNodePointer_t child_node = xmleng->GetChild(found_node);
    while(child_node) {
      string child_node_name = xmleng->GetNodeName(child_node);
      string content = xmleng->GetNodeContent(child_node);
      if(child_node_name == "nbins") { nbins = atoi(TrimSpaces(content).c_str()); }
      if(child_node_name == "min"  ) { min   = atof(TrimSpaces(content).c_str()); }
      if(child_node_name == "max"  ) { max   = atof(TrimSpaces(content).c_str()); }
      child_node = xmleng->GetNext(child_node);    
    }
#ifdef __debug_level_2__
    cout << "utils/XML2BinInfo >> nbins = " << nbins 
         << ", range = (" << min << ", " << max << ")" << endl;
#endif
    assert(nbins>0 && min<max);

    TAxis * axis = new TAxis(nbins,min,max);
    return axis;
  } 
  else   
  if(binsz == "variable") {
    cerr << "utils/XML2BinInfo >> Not handling variable size bins yet" << endl;
    exit(1);
  }
  else {
    cerr << "utils/XML2BinInfo >> Invalid bin_size attribute" << endl;
    exit(1);
  }

  return 0;
}
//___________________________________________________________________________
double t2kana::yr1numu::util::XML2Dbl(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get the specified node content as a floating-point number
//

#ifdef __debug_level_2__
  cout << "utils/XML2Dbl >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2Dbl >> No node was found at input XML node path: "
         << path << endl;
    exit(1);
  }
  string content = xmleng->GetNodeContent(found_node);
  double value = atof(content.c_str());
  return value;
}
//___________________________________________________________________________
t2kana::yr1numu::util::TH1DErrEnv *    
  t2kana::yr1numu::util::XML2TH1DErrEnv(
     TXMLEngine * eng, XMLNodePointer_t topnode, string path)
{
  TH1D * low  = XML2TH1D(eng,topnode,path + "/low" );
  TH1D * high = XML2TH1D(eng,topnode,path + "/high");

  assert(low);
  assert(high);

  TH1DErrEnv * error_envelope = new TH1DErrEnv(low,high);

  return error_envelope;
}
//___________________________________________________________________________
TMatrixDSym * t2kana::yr1numu::util::XML2TMatrixDSym(
     TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TMatrixDSym instance from the info in the input XML config file
//

#ifdef __debug_level_2__
  cout << "utils/XML2TMatrixDSym >> "
       << "Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TMatrixDSym >> "
         << "No node was found at input XML node path: " << path << endl;
    return 0;
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TMatrixDSym >> filename : " << filename << endl;
  cout << "utils/XML2TMatrixDSym >> objname  : " << objname  << endl;
#endif

  TFile file(filename.c_str(), "read");
  TMatrixDSym * matrix = (TMatrixDSym*)file.Get(objname.c_str());
  if(!matrix) {
    cerr << "utils/XML2TMatrixDSym >> No " 
         << objname << " in " << filename << endl;
    file.Close();
    return 0;
  }

  TMatrixDSym * matrix_copy = new TMatrixDSym(*matrix);

  file.Close();

  return matrix_copy;
}
//___________________________________________________________________________

