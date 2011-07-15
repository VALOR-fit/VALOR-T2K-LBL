//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 05/07/2011 - NG
   First added in in v2.0.1

*/
//________________________________________________________________________________________

#include <TMath.h>

#include "NuOsc3Flavour.h"
#include "Conventions.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
NuOsc3Flavour::NuOsc3Flavour()
{
  this->Init();
}
//________________________________________________________________________________________
NuOsc3Flavour::~NuOsc3Flavour( )
{

}
//________________________________________________________________________________________
void NuOsc3Flavour::Init()
{
  // set default baseline
  this->SetBaseline(kBaseline);

  // set a default mass hierarchy
  this->SetMassHierarchy(kNmhNormal);

  // set default density 
  this->SetDensity(kCrustDensity);
}
//________________________________________________________________________________________
void NuOsc3Flavour::SetBaseline(double length)
{
  fBaseline = length;
}
//________________________________________________________________________________________
void NuOsc3Flavour::SetDensity(double density)
{
// Convert input Earth crust density (g cm^-3) to number density of electrons Ne (cm^-3)
// Ne = density * N_A * Z / A
// where N_A = Avogadro number (mol^-1) , Z = atomic number, A = atomic weight (g/mol)
// take Z/A = 1/2 since most elements in Earth's crust are light 

  fNe = density * kNA * 0.5;
}
//________________________________________________________________________________________
void NuOsc3Flavour::SetMassHierarchy(NuMassHierarchy_t h)
{
  fMassHierarchy = h;
}
//________________________________________________________________________________________
void NuOsc3Flavour::SetOscParameters(
  double sinsq2theta12, double sinsq2theta13, double sinsq2theta23,
  double dmsq12_in, double dmsq23_in, 
  double deltaCP)
{
  dmsq12 = dmsq12_in;
  dmsq23 = dmsq23_in;

  // Note: The input value of dmsq23 is the absolute one. The sign of dmsq23 depends on 
  // the mass hierarchy. In the normal mass hierarchy, dmsq23 is positive. In the inverted
  // mass hierarchy, dmsq23 is NEGATIVE. Set the appropriate sign.
  if(fMassHierarchy == kNmhInverted) {
   dmsq23 = -dmsq23;
  }

  // The input values of dmsq23 and dmsq12  uniquely define the value of dmsq13
  dmsq13 = dmsq12 + dmsq23;

  // Code returns NaN if sin^2(mixing angle) is unphysical.
  // Don't let sin^2(mixing angle) be below 0 or above 1.
  sinsq2theta12 = TMath::Min(1.,sinsq2theta12);
  sinsq2theta12 = TMath::Max(0.,sinsq2theta12);
  sinsq2theta13 = TMath::Min(1.,sinsq2theta13);
  sinsq2theta13 = TMath::Max(0.,sinsq2theta13);
  sinsq2theta23 = TMath::Min(1.,sinsq2theta23);
  sinsq2theta23 = TMath::Max(0.,sinsq2theta23);
 
  // Calculate sin(mixingAngle) from sin^2(2*mixingAngle) 
  double s12 = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta12)));
  double s13 = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta13)));
  double s23 = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta23)));
  // Calculate cos(mixingAngle) from sin(mixingAngle)     
  double c12 = sqrt(1 - s12 * s12);
  double c13 = sqrt(1 - s13 * s13);
  double c23 = sqrt(1 - s23 * s23);

  modUe1sq = c12*c12*c13*c13;
  modUe2sq = c13*c13*s12*s12;
  modUe3sq = s13*s13;
  modUmu1sq = s12*s12*c23*c23 + 2*s12*c23*c12*s23*s13*cos(deltaCP) + c12*c12*s23*s23*s13*s13;
  modUmu2sq = c12*c12*c23*c23 - 2*c12*c23*s12*s23*s13*cos(deltaCP) + s12*s12*s23*s23*s13*s13;
  modUmu3sq = s23*s23*c13*c13;
  modUtau1sq = s12*s12*s23*s23 - 2*s12*s23*s13*c12*c23*cos(deltaCP) + s13*s13*c12*c12*c23*c23;
  modUtau2sq = s23*s23*c12*c12 + 2*s23*c12*s12*c23*s13*cos(deltaCP) + s12*s12*c23*c23*s13*s13;
  modUtau3sq = c13*c13*c23*c23;

  Re12 = c13*c13*s12*c12*(-c12*c23*c23*s12 + s12*s13*s13*s23*s23*c12 + (s12*s12*s13*s23*c23 - c12*c23*s13*s23*c12) * cos(deltaCP));
  Im12 = -c13*c13*s12*c12*c23*s13*s23*sin(deltaCP);
  arg12 = atan2(Im12, Re12);
  Re13 = -s13*c13*c13*s23*s23*c12*s13*c12 - s13*c13*c13*s23*c23*s12*c12*cos(deltaCP);
  Im13 = s13*c13*c13*s23*c23*s12*c12*sin(deltaCP);
  arg13 = atan2(Im13, Re13);
  Re23 = -s13*s13*c13*c13*s23*s23*s12*s12 + s13*c13*c13*s23*s12*c12*c23*cos(deltaCP);
  Im23 = -s13*c13*c13*s23*s12*c12*c23*sin(deltaCP);
  arg23 = atan2(Im23, Re23);
  
  UMassUdaggerVacuum[0][0][re] = c13*c13*s12*s12*dmsq12 + s13*s13*dmsq13;
  UMassUdaggerVacuum[0][0][im] = 0.0;
  UMassUdaggerVacuum[0][1][re] = c13*s12*(c12*c23 - s12*s13*s23*cos(deltaCP))*dmsq12 + s13*c13*s23*cos(deltaCP)*dmsq13;
  UMassUdaggerVacuum[0][1][im] = c13*s12*s12*s13*s23*sin(deltaCP)*dmsq12 - s13*c13*s23*sin(deltaCP)*dmsq13;
  UMassUdaggerVacuum[0][2][re] = c13*s12*(-s23*c12 - s12*c23*s13*cos(deltaCP))*dmsq12 + s13*c13*c23*cos(deltaCP)*dmsq13;
  UMassUdaggerVacuum[0][2][im] = c13*s12*s12*c23*s13*sin(deltaCP)*dmsq12 - s13*c13*c23*sin(deltaCP)*dmsq13;
  UMassUdaggerVacuum[1][0][re] = c13*s12*(c12*c23 - s12*s13*s23*cos(deltaCP))*dmsq12 + s13*c13*s23*cos(deltaCP)*dmsq13;
  UMassUdaggerVacuum[1][0][im] = -c13*s12*s12*s13*s23*sin(deltaCP)*dmsq12 + s13*c13*s23*sin(deltaCP)*dmsq13;
  UMassUdaggerVacuum[1][1][re] = (c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2.0*c12*c23*s12*s13*s23*cos(deltaCP))*dmsq12
                               + c13*c13*s23*s23*dmsq13;
  UMassUdaggerVacuum[1][1][im] = 0.0;
  UMassUdaggerVacuum[1][2][re] = (-c12*c23*s23*c12 + s12*s12*c23*s13*s13*s23 + s23*c12*s12*s13*s23*cos(deltaCP)
			       -c12*c23*s12*c23*s13*cos(deltaCP)) * dmsq12 + c13*c13*s23*c23*dmsq13;
  UMassUdaggerVacuum[1][2][im] = (s23*c12*s12*s13*s23*sin(deltaCP) + c12*c23*s12*c23*s13*sin(deltaCP)) * dmsq12;
  UMassUdaggerVacuum[2][0][re] = c13*s12*(-s23*c12 - s12*c23*s13*cos(deltaCP))*dmsq12
                               + s13*c13*c23*cos(deltaCP)*dmsq13;
  UMassUdaggerVacuum[2][0][im] = -c13*s12*s12*c23*s13*sin(deltaCP)*dmsq12 + s13*c13*c23*sin(deltaCP)*dmsq13;
  UMassUdaggerVacuum[2][1][re] = (-c12*c23*s23*c12 + s12*s12*s13*s13*s23*c23 - c12*c23*s12*c23*s13*cos(deltaCP)
			       +s23*s23*c12*s12*s13*cos(deltaCP)) * dmsq12 + c13*c13*s23*c23*dmsq13;
  UMassUdaggerVacuum[2][1][im] = (-c12*c23*s12*c23*s13*sin(deltaCP) - s23*s23*c12*s12*s13*sin(deltaCP)) * dmsq12;
  UMassUdaggerVacuum[2][2][re] = (s23*s23*c12*c12 + s12*s12*c23*c23*s13*s13 + 2.0*s23*c12*s12*c23*s13*cos(deltaCP)) * dmsq12
                               + c13*c13*c23*c23*dmsq13;
  UMassUdaggerVacuum[2][2][im] = 0.0;
}
//________________________________________________________________________________________
double NuOsc3Flavour::CalcVacuumOscProb(int nu_start, int nu_end, double energy)
{
// This method calculates the 3-flavour oscillation probabilities in vacuum.
// The formula is found in many papers and theses, e.g. equation 11 in hep-ph 0506165.
// The neutrino codes are the PDG ones and the neutrino energy is in GeV.

  OscChannelId_t osch = util::GetOscChannel(nu_start, nu_end);

  if(osch == kOscUndefined) {
    cerr << "NuOsc3Flavour: Not allowing (neutrino = " << nu_start << ") --> "
         << "(neutrino = " << nu_end << ") oscillations" << endl;
    return 0.;
  }

  double oscProb = 0.0;
  
  double sinsqDelta23 = sin((kOscDimensFactor * dmsq23 * fBaseline) / energy) 
                      * sin((kOscDimensFactor * dmsq23 * fBaseline) / energy);
  double sinsqDelta12 = sin((kOscDimensFactor * dmsq12 * fBaseline) / energy) 
                      * sin((kOscDimensFactor * dmsq12 * fBaseline) / energy);
  double sinsqDelta13 = sin((kOscDimensFactor * dmsq13 * fBaseline) / energy) 
                      * sin((kOscDimensFactor * dmsq13 * fBaseline) / energy);

  double sinDelta23 = sin((2.0 * kOscDimensFactor * dmsq23 * fBaseline) / energy);
  double sinDelta12 = sin((2.0 * kOscDimensFactor * dmsq12 * fBaseline) / energy);
  double sinDelta13 = sin((2.0 * kOscDimensFactor * dmsq13 * fBaseline) / energy);

  // nue->nue, nuebar->nuebar
  if(osch == kOscNuE2NuE || osch == kOscNuEBar2NuEBar)
  {
    oscProb = 1 - (4 * (modUe2sq*modUe1sq*sinsqDelta12
            + modUe3sq*modUe1sq*sinsqDelta13 + modUe3sq*modUe2sq*sinsqDelta23));
  }    

  // numu->nue, nuebar->numubar
  else 
  if(osch == kOscNuMu2NuE || osch == kOscNuEBar2NuMuBar)
  { 
    oscProb = 2 * (Im12*sinDelta12 + Im13*sinDelta13 + Im23*sinDelta23) 
            - 4 * (Re12*sinsqDelta12 + Re13*sinsqDelta13 + Re23*sinsqDelta23);
  }

  // numubar->nuebar, nue->numu
  else 
  if(osch == kOscNuMuBar2NuEBar || osch == kOscNuE2NuMu)
  {
    oscProb = -2 * (Im12*sinDelta12 + Im13*sinDelta13 + Im23*sinDelta23)
            - 4 * (Re12*sinsqDelta12 + Re13*sinsqDelta13 + Re23*sinsqDelta23);
  }

  // numu->numu, numubar->numubar
  else 
  if(osch == kOscNuMu2NuMu || osch == kOscNuMuBar2NuMuBar)
  {
    oscProb = 1 - (4 * (modUmu2sq*modUmu1sq*sinsqDelta12
            + modUmu3sq*modUmu1sq*sinsqDelta13 + modUmu3sq*modUmu2sq*sinsqDelta23));
  }

  // nue->nutau, nutaubar->nuebar
  else 
  if(osch == kOscNuE2NuTau || osch == kOscNuTauBar2NuEBar)
  {
    oscProb += 1 - (4 * (modUe2sq*modUe1sq*sinsqDelta12
	       + modUe3sq*modUe1sq*sinsqDelta13 + modUe3sq*modUe2sq*sinsqDelta23));
    oscProb += -2 * (Im12*sinDelta12 + Im13*sinDelta13 + Im23*sinDelta23)
   	       - 4 * (Re12*sinsqDelta12 + Re13*sinsqDelta13 + Re23*sinsqDelta23);
    oscProb = 1.0 - oscProb; 
  }

  // nuebar->nutaubar, nutau->nue
  else 
  if(osch == kOscNuEBar2NuTauBar || osch == kOscNuTau2NuE)
  {
    oscProb += 1 - (4 * (modUe2sq*modUe1sq*sinsqDelta12
      		   + modUe3sq*modUe1sq*sinsqDelta13 + modUe3sq*modUe2sq*sinsqDelta23));
    oscProb += 2 * (Im12*sinDelta12 + Im13*sinDelta13 + Im23*sinDelta23)
	    - 4 * (Re12*sinsqDelta12 + Re13*sinsqDelta13 + Re23*sinsqDelta23);
    oscProb = 1.0 - oscProb;
  }

  // numu->nutau, nutaubar->numubar
  else 
  if(osch == kOscNuMu2NuTau || osch == kOscNuTauBar2NuMuBar)
  {
    oscProb += 2 * (Im12*sinDelta12 + Im13*sinDelta13 + Im23*sinDelta23)
  	- 4 * (Re12*sinsqDelta12 + Re13*sinsqDelta13 + Re23*sinsqDelta23);
    oscProb += 1 - (4 * (modUmu2sq*modUmu1sq*sinsqDelta12
	    + modUmu3sq*modUmu1sq*sinsqDelta13 + modUmu3sq*modUmu2sq*sinsqDelta23));
    oscProb = 1.0 - oscProb;
  }

  // numubar->nutaubar, nutau->numu
  else 
  if(osch == kOscNuMuBar2NuTauBar || osch == kOscNuTau2NuMu)
  {
    oscProb += -2 * (Im12*sinDelta12 + Im13*sinDelta13 + Im23*sinDelta23)
  	- 4 * (Re12*sinsqDelta12 + Re13*sinsqDelta13 + Re23*sinsqDelta23);
    oscProb += 1 - (4 * (modUmu2sq*modUmu1sq*sinsqDelta12
	    + modUmu3sq*modUmu1sq*sinsqDelta13 + modUmu3sq*modUmu2sq*sinsqDelta23));
    oscProb = 1.0 - oscProb;
  }  

  // nutau->nutau, nutaubar->nutaubar
  else 
  if(osch == kOscNuTau2NuTau || osch == kOscNuTauBar2NuTauBar)
  {
    oscProb = 1 - (4 * (modUtau2sq*modUtau1sq*sinsqDelta12
            + modUtau3sq*modUtau1sq*sinsqDelta13 + modUtau3sq*modUtau2sq*sinsqDelta23));
  }

  return oscProb;
}
//________________________________________________________________________________________
double NuOsc3Flavour::CalcVacuumOscProbPDG(int nu_start, int nu_end, double energy)
{
// This method calculates the 3-flavour oscillation probabilities in vacuum.
// The formula is in equations 13.13 and 13.14 in the PDG review at
// http://pdg.lbl.gov/2010/reviews/rpp2010-rev-neutrino-mixing.pdf
// The neutrino codes are the PDG ones and the neutrino energy is in GeV.
//

  OscChannelId_t osch = util::GetOscChannel(nu_start, nu_end);

  if(osch == kOscUndefined) {
    cerr << "NuOsc3Flavour: Not allowing (neutrino = " << nu_start << ") --> "
         << "(neutrino = " << nu_end << ") oscillations" << endl;
    return 0.;
  }

  double oscProb = 0.0;
  
  double cosDelta12;
  double cosDelta13;
  double cosDelta23;

  if( osch == kOscNuE2NuE         || 
      osch == kOscNuEBar2NuEBar   ||
      osch == kOscNuMu2NuMu       || 
      osch == kOscNuMuBar2NuMuBar ||
      osch == kOscNuTau2NuTau     || 
      osch == kOscNuTauBar2NuTauBar  )
  {
    cosDelta12 = cos((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy);
    cosDelta13 = cos((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy);
    cosDelta23 = cos((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy);
  }
  else 
  if( osch == kOscNuMu2NuE ||
      osch == kOscNuEBar2NuMuBar )
  {
    cosDelta12 = cos(((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy) - arg12);
    cosDelta13 = cos(((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy) - arg13);
    cosDelta23 = cos(((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy) - arg23);
  }
  else 
  if( osch == kOscNuMuBar2NuEBar ||
      osch == kOscNuE2NuMu )
  {
    cosDelta12 = cos(((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy) + arg12);
    cosDelta13 = cos(((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy) + arg13);
    cosDelta23 = cos(((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy) + arg23);
  }

  if( osch == kOscNuE2NuE ||
      osch == kOscNuEBar2NuEBar )
  {
    oscProb = modUe1sq*modUe1sq + modUe2sq*modUe2sq + modUe3sq*modUe3sq + 2*modUe2sq*modUe1sq*cosDelta12
            + 2*modUe3sq*modUe1sq*cosDelta13 + 2*modUe2sq*modUe3sq*cosDelta23;
  }
  else 
  if( osch == kOscNuMuBar2NuEBar ||
      osch == kOscNuE2NuMu       ||
      osch == kOscNuMu2NuE       ||
      osch == kOscNuEBar2NuMuBar )
  {
    oscProb = modUe1sq*modUmu1sq + modUe2sq*modUmu2sq + modUe3sq*modUmu3sq
            + 2*sqrt(Re12*Re12 + Im12*Im12)*cosDelta12 + 2*sqrt(Re13*Re13 + Im13*Im13)*cosDelta13
            + 2*sqrt(Re23*Re23 + Im23*Im23)*cosDelta23;
  }
  else 
  if( osch == kOscNuMu2NuMu ||
      osch == kOscNuMuBar2NuMuBar )
  {
    oscProb = modUmu1sq*modUmu1sq + modUmu2sq*modUmu2sq + modUmu3sq*modUmu3sq + 2*modUmu2sq*modUmu1sq*cosDelta12
            + 2*modUmu3sq*modUmu1sq*cosDelta13 + 2*modUmu2sq*modUmu3sq*cosDelta23;
  }
  else 
  if( osch == kOscNuTau2NuTau ||
      osch == kOscNuTauBar2NuTauBar )
  {
    oscProb = modUtau1sq*modUtau1sq + modUtau2sq*modUtau2sq + modUtau3sq*modUtau3sq + 2*modUtau2sq*modUtau1sq*cosDelta12
            + 2*modUtau3sq*modUtau1sq*cosDelta13 + 2*modUtau2sq*modUtau3sq*cosDelta23;
  }
  else 
  if( osch == kOscNuE2NuTau ||
      osch == kOscNuTauBar2NuEBar )
  {
    cosDelta12 = cos((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy);
    cosDelta13 = cos((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy);
    cosDelta23 = cos((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy);
      
    oscProb += modUe1sq*modUe1sq + modUe2sq*modUe2sq + modUe3sq*modUe3sq + 2*modUe2sq*modUe1sq*cosDelta12
            + 2*modUe3sq*modUe1sq*cosDelta13 + 2*modUe2sq*modUe3sq*cosDelta23;

    cosDelta12 = cos(((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy) + arg12);
    cosDelta13 = cos(((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy) + arg13);
    cosDelta23 = cos(((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy) + arg23);

    oscProb += modUe1sq*modUmu1sq + modUe2sq*modUmu2sq + modUe3sq*modUmu3sq
            + 2*sqrt(Re12*Re12 + Im12*Im12)*cosDelta12 + 2*sqrt(Re13*Re13 + Im13*Im13)*cosDelta13
            + 2*sqrt(Re23*Re23 + Im23*Im23)*cosDelta23;

    oscProb = 1.0 - oscProb;
  }
  else 
  if( osch == kOscNuEBar2NuTauBar ||
      osch == kOscNuTau2NuE )
  {
    cosDelta12 = cos((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy);
    cosDelta13 = cos((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy);
    cosDelta23 = cos((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy);

    oscProb += modUe1sq*modUe1sq + modUe2sq*modUe2sq 
             + modUe3sq*modUe3sq + 2*modUe2sq*modUe1sq*cosDelta12
             + 2*modUe3sq*modUe1sq*cosDelta13 + 2*modUe2sq*modUe3sq*cosDelta23;

    cosDelta12 = cos(((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy) - arg12);
    cosDelta13 = cos(((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy) - arg13);
    cosDelta23 = cos(((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy) - arg23);

    oscProb += modUe1sq*modUmu1sq + modUe2sq*modUmu2sq + modUe3sq*modUmu3sq
           + 2*sqrt(Re12*Re12 + Im12*Im12)*cosDelta12 + 2*sqrt(Re13*Re13 + Im13*Im13)*cosDelta13
           + 2*sqrt(Re23*Re23 + Im23*Im23)*cosDelta23;

    oscProb = 1.0 - oscProb;
  }
  else 
  if( osch == kOscNuMu2NuTau ||
      osch == kOscNuTauBar2NuMuBar )
  {
    cosDelta12 = cos(((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy) - arg12);
    cosDelta13 = cos(((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy) - arg13);
    cosDelta23 = cos(((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy) - arg23);

    oscProb += modUe1sq*modUmu1sq + modUe2sq*modUmu2sq + modUe3sq*modUmu3sq
            + 2*sqrt(Re12*Re12 + Im12*Im12)*cosDelta12 + 2*sqrt(Re13*Re13 + Im13*Im13)*cosDelta13
            + 2*sqrt(Re23*Re23 + Im23*Im23)*cosDelta23;
      
    cosDelta12 = cos((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy);
    cosDelta13 = cos((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy);
    cosDelta23 = cos((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy);

    oscProb += modUmu1sq*modUmu1sq + modUmu2sq*modUmu2sq + modUmu3sq*modUmu3sq + 2*modUmu2sq*modUmu1sq*cosDelta12
            + 2*modUmu3sq*modUmu1sq*cosDelta13 + 2*modUmu2sq*modUmu3sq*cosDelta23;

    oscProb = 1.0 - oscProb;
  }
  else 
  if( osch == kOscNuMuBar2NuTauBar ||
      osch == kOscNuTau2NuMu )
  {
    cosDelta12 = cos(((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy) + arg12);
    cosDelta13 = cos(((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy) + arg13);
    cosDelta23 = cos(((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy) + arg23);

    oscProb += modUe1sq*modUmu1sq + modUe2sq*modUmu2sq + modUe3sq*modUmu3sq
            + 2*sqrt(Re12*Re12 + Im12*Im12)*cosDelta12 + 2*sqrt(Re13*Re13 + Im13*Im13)*cosDelta13
            + 2*sqrt(Re23*Re23 + Im23*Im23)*cosDelta23;    


    cosDelta12 = cos((2 * kOscDimensFactor * dmsq12 * fBaseline) / energy);
    cosDelta13 = cos((2 * kOscDimensFactor * dmsq13 * fBaseline) / energy);
    cosDelta23 = cos((2 * kOscDimensFactor * dmsq23 * fBaseline) / energy);

    oscProb += modUmu1sq*modUmu1sq + modUmu2sq*modUmu2sq + 
               modUmu3sq*modUmu3sq + 2*modUmu2sq*modUmu1sq*cosDelta12 + 
               2*modUmu3sq*modUmu1sq*cosDelta13 + 2*modUmu2sq*modUmu3sq*cosDelta23;

    oscProb = 1.0 - oscProb;
  }

  return oscProb;
}
//________________________________________________________________________________________
double NuOsc3Flavour::CalcMatterOscProb(int nu_start, int nu_end, double energy)
{
// This method calculates the 3-flavour oscillation probabilities in constant 
// density matter.
// The neutrino codes are the PDG ones and the neutrino energy is in GeV.
//

  OscChannelId_t osch = util::GetOscChannel(nu_start, nu_end);

  if(osch == kOscUndefined) {
    cerr << "NuOsc3Flavour: Not allowing (neutrino = " << nu_start << ") --> "
         << "(neutrino = " << nu_end << ") oscillations" << endl;
    return 0.;
  }

  double oscProb = 0.;

  // The mass squared matrix has diagonal elements 0, dmsq12, dmsq13; 
  // all off-diagonal elements are zero.
  // UMassUdagger is the product of the PMNS matrix, the mass squared matrix, 
  // and the dagger of the PMNS matrix.
  double UMassUdagger[3][3][2];
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {    
      UMassUdagger[i][j][re] = UMassUdaggerVacuum[i][j][re];
      if( util::IsNeutrino(nu_start) )
      {
        UMassUdagger[i][j][im] = UMassUdaggerVacuum[i][j][im]; 
      }
      else 
      if( util::IsAntiNeutrino(nu_start) )
      {
	UMassUdagger[i][j][im] = -UMassUdaggerVacuum[i][j][im];
      }
    }//j
  }//i

  // Add potential due to constant density matter to first diagonal element of UMassUdagger.
  // The factor of 1.0e9 is to convert the energy from GeV to eV.
  if( util::IsNeutrino(nu_start) )
  {
    UMassUdagger[0][0][re] = UMassUdaggerVacuum[0][0][re] + 2.0*1.0e9*energy*sqrt(2.0)*kGF_SI*fNe;
  }
  else 
  if( util::IsAntiNeutrino(nu_start) )
  {
    UMassUdagger[0][0][re] = UMassUdaggerVacuum[0][0][re] - 2.0*1.0e9*energy*sqrt(2.0)*kGF_SI*fNe;
  }

  // Calculate the eigenvalues of UMassUdagger by solving the characteristic equation, which is a cubic.
  // Use the method of del Ferro, Tartaglia and Cardano as described in Kopp (arXiv:physics/0610206, section 3.1).
  double modUMassUdaggersq[3][3];
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
        {
            modUMassUdaggersq[i][j] = UMassUdagger[i][j][re] * UMassUdagger[i][j][re]
	                            + UMassUdagger[i][j][im] * UMassUdagger[i][j][im];
        }
    }

  // Equations 22-24 in Kopp
  double c0[2];
  c0[re] = UMassUdagger[0][0][re] * modUMassUdaggersq[1][2] + UMassUdagger[1][1][re] * modUMassUdaggersq[0][2]
         + UMassUdagger[2][2][re] * modUMassUdaggersq[0][1] 
         - (UMassUdagger[0][0][re] * UMassUdagger[1][1][re] * UMassUdagger[2][2][re]
         - UMassUdagger[0][0][im] * UMassUdagger[1][1][im] * UMassUdagger[2][2][re]
         - UMassUdagger[0][0][re] * UMassUdagger[1][1][im] * UMassUdagger[2][2][im]
         - UMassUdagger[0][0][im] * UMassUdagger[1][1][re] * UMassUdagger[2][2][im])
         - 2.0 * (UMassUdagger[0][2][re] * UMassUdagger[0][1][re] * UMassUdagger[1][2][re]
         + UMassUdagger[0][2][im] * UMassUdagger[0][1][im] * UMassUdagger[1][2][re]
         - UMassUdagger[0][2][re] * UMassUdagger[0][1][im] * UMassUdagger[1][2][im]
         + UMassUdagger[0][2][im] * UMassUdagger[0][1][re] * UMassUdagger[1][2][im]);
  c0[im] = UMassUdagger[0][0][im] * modUMassUdaggersq[1][2] + UMassUdagger[1][1][im] * modUMassUdaggersq[0][2]
         + UMassUdagger[2][2][im] * modUMassUdaggersq[0][1] 
         - (UMassUdagger[0][0][re] * UMassUdagger[1][1][re] * UMassUdagger[2][2][im]
         - UMassUdagger[0][0][im] * UMassUdagger[1][1][im] * UMassUdagger[2][2][im]
         + UMassUdagger[0][0][re] * UMassUdagger[1][1][im] * UMassUdagger[2][2][re]
         + UMassUdagger[0][0][im] * UMassUdagger[1][1][re] * UMassUdagger[2][2][re]);

  double c1[2];
  c1[re] = UMassUdagger[0][0][re] * UMassUdagger[1][1][re] - UMassUdagger[0][0][im] * UMassUdagger[1][1][im]
         + UMassUdagger[0][0][re] * UMassUdagger[2][2][re] - UMassUdagger[0][0][im] * UMassUdagger[2][2][im]
         + UMassUdagger[1][1][re] * UMassUdagger[2][2][re] - UMassUdagger[1][1][im] * UMassUdagger[2][2][im]
         - modUMassUdaggersq[0][1] - modUMassUdaggersq[0][2] - modUMassUdaggersq[1][2]; 
  c1[im] = UMassUdagger[0][0][re] * UMassUdagger[1][1][im] + UMassUdagger[0][0][im] * UMassUdagger[1][1][re]
         + UMassUdagger[0][0][re] * UMassUdagger[2][2][im] + UMassUdagger[0][0][im] * UMassUdagger[2][2][re]
         + UMassUdagger[1][1][re] * UMassUdagger[2][2][im] + UMassUdagger[1][1][im] * UMassUdagger[2][2][re];

  double c2[2];
  c2[re] = -UMassUdagger[0][0][re] - UMassUdagger[1][1][re] - UMassUdagger[2][2][re];
  c2[im] = -UMassUdagger[0][0][im] - UMassUdagger[1][1][im] - UMassUdagger[2][2][im];

  // c0, c1 and c2 are always real, and their values are unaffected by the value of deltaCP.
  // Equations 26-27 in Kopp.
  double p = c2[re]*c2[re] - 3.0*c1[re];
  double q = -0.5*27.0*c0[re] - c2[re]*c2[re]*c2[re] + 0.5*9.0*c2[re]*c1[re];
  
  // Equation 32 in Kopp.
  double phi = atan2(sqrt(27.0 * (0.25*c1[re]*c1[re]*(p - c1[re]) + c0[re]*(q + (27.0/4.0)*c0[re]))), q) / 3.0;

  double eigenvalues[3];

  // Equations 33-34 in Kopp.
  // The eigenvalues of a Hermitian matrix are always real.
  if(dmsq23 > 0)
  {
      eigenvalues[0] = (sqrt(p) * 2.0 * cos(phi + 2.0*3.1415927/3.0)) / 3.0 - c2[re] / 3.0;
      eigenvalues[1] = (sqrt(p) * 2.0 * cos(phi - 2.0*3.1415927/3.0)) / 3.0 - c2[re] / 3.0;
      eigenvalues[2] = (sqrt(p) * 2.0 * cos(phi)) / 3.0 - c2[re] / 3.0;
  }
  else
  {
      eigenvalues[0] = (sqrt(p) * 2.0 * cos(phi - 2.0*3.1415927/3.0)) / 3.0 - c2[re] / 3.0;
      eigenvalues[1] = (sqrt(p) * 2.0 * cos(phi)) / 3.0 - c2[re] / 3.0;
      eigenvalues[2] = (sqrt(p) * 2.0 * cos(phi + 2.0*3.1415927/3.0)) / 3.0 - c2[re] / 3.0;
  }

  // Check that the sum of the eigenvalues is equal to the trace of the matrix
  double eigenvalueSum = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
  double trace = UMassUdagger[0][0][re] + UMassUdagger[1][1][re] + UMassUdagger[2][2][re];

  if(fabs((eigenvalueSum - trace) / trace) > 1.0e-4)
    {
    cerr<<"NuOsc3Flavour::CalcMatterOscProb() : eigenvalues may not be correct."<<endl;
    cerr<<"Fractional difference of sum of eigenvalues and trace of matrix is "
	<<fabs((eigenvalueSum - trace) / trace)<<"."<<endl;
    //return 0.0;
    }


  double UMassUdaggerMinusEv[3][3][2];
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<2; k++)
      {
	 UMassUdaggerMinusEv[i][j][k] = UMassUdagger[i][j][k];
      }//k
    }//j
  }//i

  double UMatter[3][3][2];
  double normalisation;
  double eigenvectorCheck[3][2];
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<2; j++)
    {
      eigenvectorCheck[i][j] = 0.0;
    }//j     
  }//i

  double numerator1[2];
  double numerator3[2];
  double denominator[2];

  for(int i=0; i<3; i++)
  {
    // UMassUdaggerMinusEv is UMassUdagger with the eigenvalue subtracted from all diagonal real elements
    for(int j=0; j<3; j++) {
      UMassUdaggerMinusEv[j][j][re] = UMassUdagger[j][j][re] - eigenvalues[i];
    }

    // Calculate the complex eigenvectors of UMassUdagger using an algebraic method.
    // Take the second component of each eigenvector to be 1.0 (real).
    // Calculate the first and third components in terms of the second component.
    // This calculation uses the Hermitian nature of UMassUdagger 
    numerator1[re]  = UMassUdaggerMinusEv[1][1][re]   * UMassUdaggerMinusEv[0][2][re]
                    - UMassUdaggerMinusEv[1][1][im]   * UMassUdaggerMinusEv[0][2][im]
                    - (UMassUdaggerMinusEv[1][2][re]  * UMassUdaggerMinusEv[0][1][re]
                    - UMassUdaggerMinusEv[1][2][im]   * UMassUdaggerMinusEv[0][1][im]);
    numerator1[im]  = UMassUdaggerMinusEv[1][1][re]   * UMassUdaggerMinusEv[0][2][im]
                    + UMassUdaggerMinusEv[1][1][im]   * UMassUdaggerMinusEv[0][2][re]
                    - (UMassUdaggerMinusEv[1][2][re]  * UMassUdaggerMinusEv[0][1][im]
                    + UMassUdaggerMinusEv[1][2][im]   * UMassUdaggerMinusEv[0][1][re]);
    numerator3[re]  = UMassUdaggerMinusEv[0][1][re]   * UMassUdaggerMinusEv[0][1][re]
                    + UMassUdaggerMinusEv[0][1][im]   * UMassUdaggerMinusEv[0][1][im]
                    - (UMassUdaggerMinusEv[0][0][re]  * UMassUdaggerMinusEv[1][1][re]
                    - UMassUdaggerMinusEv[0][0][im]   * UMassUdaggerMinusEv[1][1][im]);
    numerator3[im]  = -(UMassUdaggerMinusEv[0][0][re] * UMassUdaggerMinusEv[1][1][im]
	            + UMassUdaggerMinusEv[0][0][im]   * UMassUdaggerMinusEv[1][1][re]);
    denominator[re] = UMassUdaggerMinusEv[0][0][re]   * UMassUdaggerMinusEv[1][2][re]
                    - UMassUdaggerMinusEv[0][0][im]   * UMassUdaggerMinusEv[1][2][im]
                    - (UMassUdaggerMinusEv[0][2][re]  * UMassUdaggerMinusEv[0][1][re]
                    + UMassUdaggerMinusEv[0][2][im]   * UMassUdaggerMinusEv[0][1][im]);
    denominator[im] = UMassUdaggerMinusEv[0][0][re]   * UMassUdaggerMinusEv[1][2][im]
                    + UMassUdaggerMinusEv[0][0][im]   * UMassUdaggerMinusEv[1][2][re]
                    - (-UMassUdaggerMinusEv[0][2][re] * UMassUdaggerMinusEv[0][1][im]
                    + UMassUdaggerMinusEv[0][2][im]   * UMassUdaggerMinusEv[0][1][re]);

    // The columns of UMatter are the eigenvectors of UMassUdagger
    UMatter[0][i][re] = (numerator1[re]  * denominator[re] + numerator1[im]  * denominator[im])
                      / (denominator[re] * denominator[re] + denominator[im] * denominator[im]); 
    UMatter[0][i][im] = (-numerator1[re] * denominator[im] + numerator1[im]  * denominator[re])
                      / (denominator[re] * denominator[re] + denominator[im] * denominator[im]);
    UMatter[1][i][re] = 1.0;
    UMatter[1][i][im] = 0.0;
    UMatter[2][i][re] = (numerator3[re]  * denominator[re] + numerator3[im]  * denominator[im])
                      / (denominator[re] * denominator[re] + denominator[im] * denominator[im]);
    UMatter[2][i][im] = (-numerator3[re] * denominator[im] + numerator3[im]  * denominator[re])
                      / (denominator[re] * denominator[re] + denominator[im] * denominator[im]);

    // Now normalise the eigenvector to unit length 
    normalisation = sqrt(UMatter[0][i][re] * UMatter[0][i][re] + UMatter[0][i][im] * UMatter[0][i][im]
                       + UMatter[1][i][re] * UMatter[1][i][re] + UMatter[1][i][im] * UMatter[1][i][im]  
                       + UMatter[2][i][re] * UMatter[2][i][re] + UMatter[2][i][im] * UMatter[2][i][im]);
    for(int k=0; k<3; k++)
    {
       UMatter[k][i][re] = UMatter[k][i][re] / normalisation;
       UMatter[k][i][im] = UMatter[k][i][im] / normalisation;
    }
  
    // Check the eigenvector is correct.
    // Multiplying any row of UMassUdaggerMinusEv by the column eigenvector should give zero   
    for(int m=0; m<3; m++)
    {
      eigenvectorCheck[i][re] += UMassUdaggerMinusEv[i][m][re] * UMatter[m][i][re]
                              - UMassUdaggerMinusEv[i][m][im] * UMatter[m][i][im];
      eigenvectorCheck[i][im] += UMassUdaggerMinusEv[i][m][re] * UMatter[m][i][im]
                              + UMassUdaggerMinusEv[i][m][im] * UMatter[m][i][re];
    }

 }//end of for(int i=0; i<3; i++)

 for(int i=0; i<3; i++)
   {
   if(fabs(eigenvectorCheck[i][re]) > 1.0e-10)
     {
     cerr<<"NuOsc3Flavour::CalcMatterOscProb() : eigenvector may not be correct."<<endl;
     cerr<<"Real part of check is "<<eigenvectorCheck[i][re]<<"."<<endl;
     //return 0.0;
     }
   if(fabs(eigenvectorCheck[i][im]) > 1.0e-10)
     {
     cerr<<"NuOsc3Flavour::CalcMatterOscProb() : eigenvector may not be correct."<<endl;
     cerr<<"Imaginary part of check is "<<eigenvectorCheck[i][im]<<"."<<endl;
     //return 0.0;
     }
   }//end of for(int i=0; i<3; i++)


 // Check that the complex eigenvectors are orthogonal  
 double scalarProduct[3][2];
 for(int i=0; i<3; i++)
 {
    for(int j=0; j<2; j++)
    {
      scalarProduct[i][j] = 0.0;
    }      
 }

 // For complex vectors, the scalar product is defined as the complex product 
 // of the first vector and the complex conjugate of the second vector
 for(int i=0; i<2; i++)
 {
    for(int j=0; j<3; j++)
    {
      scalarProduct[i][re] +=  UMatter[j][i][re] * UMatter[j][i+1][re]
	                   +   UMatter[j][i][im] * UMatter[j][i+1][im];
      scalarProduct[i][im] += -UMatter[j][i][re] * UMatter[j][i+1][im]
	                   +   UMatter[j][i][im] * UMatter[j][i+1][re];
    }
 }
  
 for(int j=0; j<3; j++)
 {
    scalarProduct[2][re] +=  UMatter[j][0][re] * UMatter[j][2][re]
                          +  UMatter[j][0][im] * UMatter[j][2][im];
    scalarProduct[2][im] += -UMatter[j][0][re] * UMatter[j][2][im]
                           + UMatter[j][0][im] * UMatter[j][2][re];
 }

 for(int i=0; i<3; i++)
   {
   if(fabs(scalarProduct[i][re]) > 1.0e-4)
     {
     cerr<<"NuOsc3Flavour::CalcMatterOscProb() : eigenvectors may not be orthogonal."<<endl;
     cerr<<"Real part of scalar product is "<<scalarProduct[i][re]<<"."<<endl;
     //return 0.0;
     }        
   if(fabs(scalarProduct[i][im]) > 1.0e-4)
     {
     cerr<<"NuOsc3Flavour::CalcMatterOscProb() : eigenvectors may not be orthogonal."<<endl;
     cerr<<"Imaginary part of scalar product is "<<scalarProduct[i][im]<<"."<<endl;
     //return 0.0;
     }
   }//end of for(int i=0; i<3; i++)

 double UMatterDagger[3][3][2];

 for(int i=0; i<3; i++)
 {
    for(int j=0; j<3; j++)
    {
      UMatterDagger[i][j][re] = UMatter[j][i][re];
      UMatterDagger[i][j][im] = -UMatter[j][i][im];
    }
 }

 //check that the Hermitian conjugate of UMatter is equal to its inverse
 double inverseCheck[3][3][2];
 for(int i=0; i<3; i++)
   {
   for(int j=0; j<3; j++)
     {
     inverseCheck[i][j][re] = 0.0;
     inverseCheck[i][j][im] = 0.0;
     }
   }

 for(int i=0; i<3; i++)
   {
   for(int j=0; j<3; j++)
     {
     for(int k=0; k<3; k++)
       {
       inverseCheck[i][j][re] += UMatter[i][k][re] * UMatterDagger[k][j][re]
                              - UMatter[i][k][im] * UMatterDagger[k][j][im];
       inverseCheck[i][j][im] += UMatter[i][k][re] * UMatterDagger[k][j][im]
                              + UMatter[i][k][im] * UMatterDagger[k][j][re];
       }
     }
   }

 for(int i=0; i<3; i++)
   {
   for(int j=0; j<3; j++)
     {
     if(i == j)
       {
       if(fabs(inverseCheck[i][j][re] - 1.0) > 1.0e-4)
	 {
	 cerr<<"NuOsc3Flavour::CalcMatterOscProb() : Hermitian conjugate of matter mixing matrix may not be its inverse."<<endl;
	 cerr<<"Real part of diagonal element of product - 1.0 = "<<fabs(inverseCheck[i][j][re] - 1.0)<<"."<<endl;
	 //return 0.0;
	 }
       if(fabs(inverseCheck[i][j][im]) > 1.0e-4)
	 {
	 cerr<<"NuOsc3Flavour::CalcMatterOscProb() : Hermitian conjugate of matter mixing matrix may not be its inverse."<<endl;
	 cerr<<"Imaginary part of diagonal element of product = "<<inverseCheck[i][j][im]<<"."<<endl;
	 //return 0.0;
	 }      
       }
     else
       {
       if(fabs(inverseCheck[i][j][re]) > 1.0e-4)
	 {
	 cerr<<"NuOsc3Flavour::CalcMatterOscProb() : Hermitian conjugate of matter mixing matrix may not be its inverse."<<endl;
	 cerr<<"Real part of off-diagonal element of product = "<<inverseCheck[i][j][re]<<"."<<endl;
	 //return 0.0;
	 }
       if(fabs(inverseCheck[i][j][im]) > 1.0e-4)
	 {
	 cerr<<"NuOsc3Flavour::CalcMatterOscProb() : Hermitian conjugate of matter mixing matrix may not be its inverse."<<endl; 
	 cerr<<"Imaginary part of off-diagonal element of product = "<<inverseCheck[i][j][im]<<"."<<endl;
	 //return 0.0;
	 }              
       }  
     }//end of for(int j=0; j<3; j++)
   }//end of for(int i=0; i<3; i++)

 double initialFlavourState[3][2];
 for(int i=0; i<3; i++)
 {
    initialFlavourState[i][re] = 0.0;
    initialFlavourState[i][im] = 0.0;
 }

 // The neutrino is formed in a definite flavour state
 if(nu_start == kPdgNuE || nu_start == kPdgNuEBar)
 {
    initialFlavourState[0][re] = 1.0;
 }
 else 
 if(nu_start == kPdgNuMu || nu_start == kPdgNuMuBar)
 {
    initialFlavourState[1][re] = 1.0;
 }
 else 
 if(nu_start == kPdgNuTau || nu_start == kPdgNuTauBar)
 {
    initialFlavourState[2][re] = 1.0;
 }

 double initialMassState[3][2];
 for(int i=0; i<3; i++)
 {
    for(int j=0; j<2; j++)
    {
      initialMassState[i][j] = 0.0;
    }
 }

 // Convert the initial flavour state to initial mass states by multiplying by UMatterDagger
 for(int i=0; i<3; i++)
 {
    for(int j=0; j<3; j++)
    {
      initialMassState[i][re] += UMatterDagger[i][j][re] * initialFlavourState[j][re]
	                       - UMatterDagger[i][j][im] * initialFlavourState[j][im];
      initialMassState[i][im] += UMatterDagger[i][j][re] * initialFlavourState[j][im]
                               + UMatterDagger[i][j][im] * initialFlavourState[j][re];
    }
 }

 // Now propagate the mass eigenstates to the far detector
 double finalMassState[3][2];

 // For the first mass eigenstate, the difference in the squared masses is zero.
 // Therefore only retain terms multiplied by cos(0), and drop those multiplied by sin(0).
 finalMassState[0][re] = initialMassState[0][re];
 finalMassState[0][im] = initialMassState[0][im];
 finalMassState[1][re] = initialMassState[1][re] *  cos((2.0 * kOscDimensFactor * (eigenvalues[1] - eigenvalues[0]) * fBaseline) / energy)
                       - initialMassState[1][im] * -sin((2.0 * kOscDimensFactor * (eigenvalues[1] - eigenvalues[0]) * fBaseline) / energy);  
 finalMassState[1][im] = initialMassState[1][re] * -sin((2.0 * kOscDimensFactor * (eigenvalues[1] - eigenvalues[0]) * fBaseline) / energy)
                       + initialMassState[1][im] *  cos((2.0 * kOscDimensFactor * (eigenvalues[1] - eigenvalues[0]) * fBaseline) / energy);
 finalMassState[2][re] = initialMassState[2][re] *  cos((2.0 * kOscDimensFactor * (eigenvalues[2] - eigenvalues[0]) * fBaseline) / energy)
                       - initialMassState[2][im] * -sin((2.0 * kOscDimensFactor * (eigenvalues[2] - eigenvalues[0]) * fBaseline) / energy);
 finalMassState[2][im] = initialMassState[2][re] * -sin((2.0 * kOscDimensFactor * (eigenvalues[2] - eigenvalues[0]) * fBaseline) / energy)
                       + initialMassState[2][im] *  cos((2.0 * kOscDimensFactor * (eigenvalues[2] - eigenvalues[0]) * fBaseline) / energy);

 double finalFlavourState[3][2];
 for(int i=0; i<3; i++)
 {
    for(int j=0; j<2; j++)
    {
      finalFlavourState[i][j] = 0.0;
    }
 }

 // Convert final mass states to flavour states by multiplying by UMatter 
 for(int i=0; i<3; i++)
 {
    for(int j=0; j<3; j++)
    {
      finalFlavourState[i][re] += UMatter[i][j][re] * finalMassState[j][re]
	                       - UMatter[i][j][im] * finalMassState[j][im];
      finalFlavourState[i][im] += UMatter[i][j][re] * finalMassState[j][im]
	                       + UMatter[i][j][im] * finalMassState[j][re];
    }
 }

 // Each component of finalFlavourState is the amplitude for the neutrino to interact 
 // in the far detector as that flavour.
 // The oscillation probabilities are the moduli squared of the amplitudes
 if(nu_end == kPdgNuE || nu_end == kPdgNuEBar)
 {
    oscProb = finalFlavourState[0][re]*finalFlavourState[0][re] + 
              finalFlavourState[0][im]*finalFlavourState[0][im];
 }
 else 
 if(nu_end == kPdgNuMu || nu_end == kPdgNuMuBar)
 {
    oscProb = finalFlavourState[1][re]*finalFlavourState[1][re] + 
              finalFlavourState[1][im]*finalFlavourState[1][im];
 }
 else 
 if(nu_end == kPdgNuTau || nu_end == kPdgNuTauBar)
 {
    oscProb = finalFlavourState[2][re]*finalFlavourState[2][re] + 
              finalFlavourState[2][im]*finalFlavourState[2][im];
 }

 return oscProb;
}
//________________________________________________________________________________________
