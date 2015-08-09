//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/

#include <TMath.h>

#include "Enumerations.h"
#include "Utilities.h"
#include "Binning.h"

using namespace valor;
using namespace valor::conv;
//________________________________________________________________________________________
int valor::util::NumErecoBins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kNErecoBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCNErecoBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCNErecoBins1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kNErecoBinsMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kNErecoBinsMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::ErecoBinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kErecoBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kErecoBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCErecoBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCErecoBinEdges1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kErecoBinEdgesMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kErecoBinEdgesMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double valor::util::ErecoBinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kErecoBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike      ) : return kErecoBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCErecoBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCErecoBinEdges1Re[i];  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kErecoBinEdgesMultiRmu[i]; break;
     case ( kSKsFCFVMultiRingELike  ) : return kErecoBinEdgesMultiRe[i];  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int valor::util::ErecoBinContains(valor::conv::SuperkSample_t s, double Ereco)
{
  int nedges = util::NumErecoBins(s);
  const double * edges = util::ErecoBinEdges(s);
  if (Ereco > edges[nedges] || Ereco < edges[0]) return -1;
  int bin =  TMath::BinarySearch(nedges,edges,Ereco);
  if (bin >= util::NumErecoBins(s)) return -1;
  return bin; //return in range 0->nbins-1
}
//___________________________________________________________________________
int valor::util::NumXsecErecoBins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kXsNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kXsNErecoBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCXsNErecoBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCXsNErecoBins1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kXsNErecoBinsMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kXsNErecoBinsMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::XsecErecoBinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kXsErecoBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kXsErecoBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCXsErecoBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCXsErecoBinEdges1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kXsErecoBinEdgesMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kXsErecoBinEdgesMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double valor::util::XsecErecoBinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kXsErecoBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike      ) : return kXsErecoBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCXsErecoBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCXsErecoBinEdges1Re[i];  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kXsErecoBinEdgesMultiRmu[i]; break;
     case ( kSKsFCFVMultiRingELike  ) : return kXsErecoBinEdgesMultiRe[i];  break;

     default                          : return 0;
  }
  return 0;
}















// NEW 2D TEMPLATES
int valor::util::NumQ2Bins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kNQ2Bins1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kNQ2Bins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCNQ2Bins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCNQ2Bins1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kNQ2BinsMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kNQ2BinsMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::Q2BinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kQ2BinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kQ2BinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCQ2BinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCQ2BinEdges1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kQ2BinEdgesMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kQ2BinEdgesMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double valor::util::Q2BinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kQ2BinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike      ) : return kQ2BinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCQ2BinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCQ2BinEdges1Re[i];  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kQ2BinEdgesMultiRmu[i]; break;
     case ( kSKsFCFVMultiRingELike  ) : return kQ2BinEdgesMultiRe[i];  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int valor::util::Q2BinContains(valor::conv::SuperkSample_t s, double Q2)
{
  int nedges = util::NumQ2Bins(s);
  const double * edges = util::Q2BinEdges(s);
  if (Q2 > edges[nedges] || Q2 < edges[0]) return -1;
  int bin =  TMath::BinarySearch(nedges,edges,Q2);
  if (bin >= util::NumQ2Bins(s)) return -1;
  return bin; //return in range 0->nbins-1
}
//___________________________________________________________________________
int valor::util::NumXsecQ2Bins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kXsNQ2Bins1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kXsNQ2Bins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCXsNQ2Bins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCXsNQ2Bins1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kXsNQ2BinsMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kXsNQ2BinsMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::XsecQ2BinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kXsQ2BinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kXsQ2BinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCXsQ2BinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCXsQ2BinEdges1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kXsQ2BinEdgesMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kXsQ2BinEdgesMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double valor::util::XsecQ2BinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kXsQ2BinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike      ) : return kXsQ2BinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCXsQ2BinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCXsQ2BinEdges1Re[i];  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kXsQ2BinEdgesMultiRmu[i]; break;
     case ( kSKsFCFVMultiRingELike  ) : return kXsQ2BinEdgesMultiRe[i];  break;

     default                          : return 0;
  }
  return 0;
}
//


//___________________________________________________________________________
string valor::util::GetXMLDatasetParamsName()
{
#ifndef __USE_2D_TEMPLATES__ //1D templates
  return "ereco_values";
#else //2D templates
#ifdef __ERECO_VS_Q2_TEMPLATES__
  return "ereco_q2_values";
#else
  return "ereco_theta_values";
#endif //__ERECO_VS_Q2_TEMPLATES__
#endif //__USE_2D_TEMPLATES__
}
//___________________________________________________________________________
string valor::util::GetAxisNameCompact(int iaxis)
{
  if(iaxis == 0)
    return "Ereco";
#ifdef __USE_2D_TEMPLATES__ //2D templates
  if(iaxis == 1)
#ifdef __ERECO_VS_Q2_TEMPLATES__
    return "Q2";
#else
  return "Theta";
#endif //__ERECO_VS_Q2_TEMPLATES__
#endif //__USE_2D_TEMPLATES__
  else
    return "???";  
}
//___________________________________________________________________________
string valor::util::GetAxisName(int iaxis)
{
  if(iaxis == 0)
    return "Reconstructed energy (GeV)";
#ifdef __USE_2D_TEMPLATES__ //2D templates
  if(iaxis == 1)
#ifdef __ERECO_VS_Q2_TEMPLATES__
    return "Q2 (GeV)";
#else
  return "Theta";
#endif //__ERECO_VS_Q2_TEMPLATES__
#endif //__USE_2D_TEMPLATES__
  else
    return "???";  
}
//___________________________________________________________________________



















//___________________________________________________________________________
int valor::util::NumEtrueBins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kNEtrueBins1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kNEtrueBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCNEtrueBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCNEtrueBins1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kNEtrueBinsMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kNEtrueBinsMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::EtrueBinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kEtrueBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike      ) : return kEtrueBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCEtrueBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCEtrueBinEdges1Re;  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kEtrueBinEdgesMultiRmu; break;
     case ( kSKsFCFVMultiRingELike  ) : return kEtrueBinEdgesMultiRe;  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double  valor::util::EtrueBinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return kEtrueBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike      ) : return kEtrueBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return kRHCEtrueBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike   ) : return kRHCEtrueBinEdges1Re[i];  break;
     case ( kSKsFCFVMultiRingMuLike ) : return kEtrueBinEdgesMultiRmu[i]; break;
     case ( kSKsFCFVMultiRingELike  ) : return kEtrueBinEdgesMultiRe[i];  break;

     default                          : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int valor::util::NumVariableXsEtrueBins(valor::conv::SuperkSample_t s,
                                         valor::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  //need to start counter at zero. Want nbins. Looping over bin edges.
  int nbins = -1;

  for(int it = 0; it <= valor::conv::kXsNEtrueBins; it++) {
    if(valor::conv::kXsEtrueBinEdges[is][im][it] < 0)
      break;
    nbins++;
    }

  return nbins;
}
//___________________________________________________________________________
const double * valor::util::VariableXsEtrueBinEdges(valor::conv::SuperkSample_t s,
                                                     valor::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  return valor::conv::kXsEtrueBinEdges[is][im];
}
//___________________________________________________________________________
int valor::util::NumVariableMCEtrueBins(valor::conv::SuperkSample_t s,
                                         valor::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  //need to start counter at zero. Want nbins. Looping over bin edges.
  int nbins = -1;

  for(int it = 0; it <= valor::conv::kMCNEtrueBins; it++) {
    if(valor::conv::kMCEtrueBinEdges[is][im][it] < 0)
      break;
    nbins++;
    }

  return nbins;
}
//___________________________________________________________________________
const double * valor::util::VariableMCEtrueBinEdges(valor::conv::SuperkSample_t s,
                                                     valor::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  return valor::conv::kMCEtrueBinEdges[is][im];
}
//___________________________________________________________________________
