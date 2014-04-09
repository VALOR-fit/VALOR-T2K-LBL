//____________________________________________________________________________
/*!

\class     t2kana::utils

\brief     XML utilities

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_xml_utils_h_
#define _t2kana_xml_utils_h_

#include <string>
#include <vector>
#include <complex>

#include <TXMLEngine.h>
#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TSpline.h>

class TGraph;
class TH1;
class TH1F;
class TH1D;
class TH2D;
class TAxis;

using std::string;
using std::vector;

namespace t2kana {
namespace util {

  XMLNodePointer_t FindNode          (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  bool             NodeExists        (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);

  double           XML2Dbl           (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  vector<double>   XML2DblArray      (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  string           XML2String        (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);

  string           XML2RootFileName  (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  string           XML2RootObjName   (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  string           XML2RootObjType   (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);

  TFile *          XML2TFile         (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TH1F *           XML2TH1F          (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TH1D *           XML2TH1D          (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TH2D *           XML2TH2D          (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TTree *          XML2TTree         (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TVectorD *       XML2TVectorD      (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TMatrixDSym *    XML2TMatrixDSym   (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TMatrixD *       XML2TMatrixD      (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TMatrixD *       XML2TH2D2TMatrixD (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TSpline3 *       XML2TSpline3      (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TAxis *          XML2BinInfo       (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  
} // util
} // t2kana
  
#endif


