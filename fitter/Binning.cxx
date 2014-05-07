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
     case ( kSKsFCFV1RingMuLike    ) : return kNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike     ) : return kNErecoBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCNErecoBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCNErecoBins1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::ErecoBinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kErecoBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike     ) : return kErecoBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCErecoBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCErecoBinEdges1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double valor::util::ErecoBinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kErecoBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike     ) : return kErecoBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCErecoBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCErecoBinEdges1Re[i];  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int valor::util::ErecoBinContains(valor::conv::SuperkSample_t s, double Ereco)
{
  int nedges = util::NumErecoBins(s)+1;
  const double * edges = util::ErecoBinEdges(s);
  if (Ereco > edges[nedges] || Ereco < edges[0]) return -1;
  int bin =  TMath::BinarySearch(nedges,edges,Ereco);
  if (bin == util::NumErecoBins(s)+1) return -1;
  return bin;
}
//___________________________________________________________________________
int valor::util::NumXsecErecoBins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kXsNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike     ) : return kXsNErecoBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCXsNErecoBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCXsNErecoBins1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::XsecErecoBinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kXsErecoBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike     ) : return kXsErecoBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCXsErecoBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCXsErecoBinEdges1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double valor::util::XsecErecoBinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kXsErecoBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike     ) : return kXsErecoBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCXsErecoBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCXsErecoBinEdges1Re[i];  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int valor::util::NumEtrueBins(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kNEtrueBins1Rmu; break;
     case ( kSKsFCFV1RingELike     ) : return kNEtrueBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCNEtrueBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCNEtrueBins1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * valor::util::EtrueBinEdges(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kEtrueBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike     ) : return kEtrueBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCEtrueBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCEtrueBinEdges1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double  valor::util::EtrueBinEdges(valor::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return kEtrueBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike     ) : return kEtrueBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCEtrueBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCEtrueBinEdges1Re[i];  break;

     default                      : return 0;
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
