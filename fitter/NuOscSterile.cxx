
//________________________________________________________________________________________
/*

 >> ChangeLog:


*/
//________________________________________________________________________________________

//#define __check_diagonalisation__

#include <TMath.h>

#include "NuOscSterile.h"
#include "Conventions.h"
#include "Utilities.h"
#include "Opt.h"

#ifdef __GSL_ENABLED__ 
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#endif

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

//#define DEBUG_INFO

//________________________________________________________________________________________
NuOscSterile::NuOscSterile()
{
  this->Init();
}
//________________________________________________________________________________________
NuOscSterile::~NuOscSterile()
{

}
//________________________________________________________________________________________
void NuOscSterile::Init()
{
  // set default baseline
  this->SetBaseline(kBaseline);

  // set the default mass hierarchy 
  NuMassHierarchy_t nmh = Opt::Instance()->gMassHierarchyAssumed;
  this->SetMassHierarchy(nmh);

  // set default density 
  this->SetDensity(kCrustDensity);

  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      sinMixingAngle[i][j] = 0.0;
      cosMixingAngle[i][j] = 1.0;
      }
    }

  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      dmsq[i][j] = 0.0;
      }
    }  

  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      deltaCP[i][j] = 0.0;
      }
    }

  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      for (int k=0; k<nFlav; k++)
	{
	for (int l=0; l<nFlav; l++)
	  {
	  if(k == l)
	    {
	    rotationMatrix[i][j][k][l][re] = 1.0;
	    rotationMatrix[i][j][k][l][im] = 0.0;
	    }
	  else
	    {
	    rotationMatrix[i][j][k][l][re] = 0.0;
	    rotationMatrix[i][j][k][l][im] = 0.0;
	    }
	  }
	}
      }
    }

  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      UVacuum[i][j][re] = 0.0;
      UVacuum[i][j][im] = 0.0;
      }
    }

}
//________________________________________________________________________________________
void NuOscSterile::SetBaseline(double length)
{
  fBaseline = length;
}
//________________________________________________________________________________________
void NuOscSterile::SetDensity(double density)
{
  // Convert input Earth crust density (g cm^-3) to number density of electrons Ne (cm^-3)
  // Ne = density * N_A * Z / A
  // where N_A = Avogadro number (mol^-1) , Z = atomic number, A = atomic weight (g/mol)
  // take Z/A = 1/2 since most elements in Earth's crust are light 

  fNe = density * kNA * 0.5;
}
//________________________________________________________________________________________
void NuOscSterile::SetMassHierarchy(NuMassHierarchy_t h)
{
#ifdef DEBUG_INFO
  cout << "NuOscSterile >> Mass hierarchy assumed: " << util::AsString(h) << endl;
#endif
  fMassHierarchy = h;
}
//________________________________________________________________________________________
void NuOscSterile::SetMixingAngles(double sinsqtheta12, double sinsqtheta13, double sinsqtheta14, double sinsqtheta15,
				   double sinsqtheta23, double sinsqtheta24, double sinsqtheta25, double sinsqtheta34, double sinsqtheta35)
{
  //Neglect oscillations between sterile states 4 and 5 since they give no observable effects, i.e. theta45 = 0.

  // Code returns NaN if sin^2(mixing angle) is unphysical. 
  // Don't let sin^2(mixing angle) be below 0 or above 1. 
  sinsqtheta12 = TMath::Min(1.,sinsqtheta12);
  sinsqtheta12 = TMath::Max(0.,sinsqtheta12);
  sinsqtheta13 = TMath::Min(1.,sinsqtheta13);
  sinsqtheta13 = TMath::Max(0.,sinsqtheta13);
  sinsqtheta14 = TMath::Min(1.,sinsqtheta14);
  sinsqtheta14 = TMath::Max(0.,sinsqtheta14);
  sinsqtheta15 = TMath::Min(1.,sinsqtheta15);
  sinsqtheta15 = TMath::Max(0.,sinsqtheta15);
  sinsqtheta23 = TMath::Min(1.,sinsqtheta23);
  sinsqtheta23 = TMath::Max(0.,sinsqtheta23);
  sinsqtheta24 = TMath::Min(1.,sinsqtheta24);
  sinsqtheta24 = TMath::Max(0.,sinsqtheta24);
  sinsqtheta25 = TMath::Min(1.,sinsqtheta25);
  sinsqtheta25 = TMath::Max(0.,sinsqtheta25);
  sinsqtheta34 = TMath::Min(1.,sinsqtheta34);
  sinsqtheta34 = TMath::Max(0.,sinsqtheta34);
  sinsqtheta35 = TMath::Min(1.,sinsqtheta35);
  sinsqtheta35 = TMath::Max(0.,sinsqtheta35);

  // Calculate sin and cos(mixingAngle) from sin^2(mixingAngle) 
  sinMixingAngle[0][1] = sqrt(sinsqtheta12);
  sinMixingAngle[0][2] = sqrt(sinsqtheta13);
  sinMixingAngle[0][3] = sqrt(sinsqtheta14);
  sinMixingAngle[0][4] = sqrt(sinsqtheta15);
  sinMixingAngle[1][2] = sqrt(sinsqtheta23);
  sinMixingAngle[1][3] = sqrt(sinsqtheta24);
  sinMixingAngle[1][4] = sqrt(sinsqtheta25);
  sinMixingAngle[2][3] = sqrt(sinsqtheta34);
  sinMixingAngle[2][4] = sqrt(sinsqtheta35);
  cosMixingAngle[0][1] = sqrt(1.0 - sinsqtheta12);
  cosMixingAngle[0][2] = sqrt(1.0 - sinsqtheta13);
  cosMixingAngle[0][3] = sqrt(1.0 - sinsqtheta14);
  cosMixingAngle[0][4] = sqrt(1.0 - sinsqtheta15);
  cosMixingAngle[1][2] = sqrt(1.0 - sinsqtheta23);
  cosMixingAngle[1][3] = sqrt(1.0 - sinsqtheta24);
  cosMixingAngle[1][4] = sqrt(1.0 - sinsqtheta25);
  cosMixingAngle[2][3] = sqrt(1.0 - sinsqtheta34);
  cosMixingAngle[2][4] = sqrt(1.0 - sinsqtheta35);

  usingOscParams = true;

}
//________________________________________________________________________________________
void NuOscSterile::SetMixingDoubleAnglesOct1(double sinsq2theta12, double sinsq2theta13, double sinsq2theta14, double sinsq2theta15,
                                   double sinsq2theta23, double sinsq2theta24, double sinsq2theta25, double sinsq2theta34, double sinsq2theta35)
{
  //Neglect oscillations between sterile states 4 and 5 since they give no observable effects, i.e. theta45 = 0.

  // Code returns NaN if sin^2(2*mixing angle) is unphysical.
  // Don't let sin^2(2*mixing angle) be below 0 or above 1.
  sinsq2theta12 = TMath::Min(1.,sinsq2theta12);
  sinsq2theta12 = TMath::Max(0.,sinsq2theta12);
  sinsq2theta13 = TMath::Min(1.,sinsq2theta13);
  sinsq2theta13 = TMath::Max(0.,sinsq2theta13);
  sinsq2theta14 = TMath::Min(1.,sinsq2theta14);
  sinsq2theta14 = TMath::Max(0.,sinsq2theta14);
  sinsq2theta15 = TMath::Min(1.,sinsq2theta15);
  sinsq2theta15 = TMath::Max(0.,sinsq2theta15);
  sinsq2theta23 = TMath::Min(1.,sinsq2theta23);
  sinsq2theta23 = TMath::Max(0.,sinsq2theta23);
  sinsq2theta24 = TMath::Min(1.,sinsq2theta24);
  sinsq2theta24 = TMath::Max(0.,sinsq2theta24);
  sinsq2theta25 = TMath::Min(1.,sinsq2theta25);
  sinsq2theta25 = TMath::Max(0.,sinsq2theta25);
  sinsq2theta34 = TMath::Min(1.,sinsq2theta34);
  sinsq2theta34 = TMath::Max(0.,sinsq2theta34);
  sinsq2theta35 = TMath::Min(1.,sinsq2theta35);
  sinsq2theta35 = TMath::Max(0.,sinsq2theta35);

  //Calculate sin and cos(mixingAngle) from sin^2(2*mixingAngle)
  sinMixingAngle[0][1] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta12)));
  sinMixingAngle[0][2] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta13)));
  sinMixingAngle[0][3] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta14)));
  sinMixingAngle[0][4] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta15)));
  //theta23 is in first octant
  sinMixingAngle[1][2] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta23)));
  sinMixingAngle[1][3] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta24)));
  sinMixingAngle[1][4] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta25)));
  sinMixingAngle[2][3] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta34)));
  sinMixingAngle[2][4] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta35)));
  cosMixingAngle[0][1] = sqrt(1.0 - sinMixingAngle[0][1] * sinMixingAngle[0][1]);
  cosMixingAngle[0][2] = sqrt(1.0 - sinMixingAngle[0][2] * sinMixingAngle[0][2]);
  cosMixingAngle[0][3] = sqrt(1.0 - sinMixingAngle[0][3] * sinMixingAngle[0][3]);
  cosMixingAngle[0][4] = sqrt(1.0 - sinMixingAngle[0][4] * sinMixingAngle[0][4]);
  cosMixingAngle[1][2] = sqrt(1.0 - sinMixingAngle[1][2] * sinMixingAngle[1][2]);
  cosMixingAngle[1][3] = sqrt(1.0 - sinMixingAngle[1][3] * sinMixingAngle[1][3]);
  cosMixingAngle[1][4] = sqrt(1.0 - sinMixingAngle[1][4] * sinMixingAngle[1][4]);
  cosMixingAngle[2][3] = sqrt(1.0 - sinMixingAngle[2][3] * sinMixingAngle[2][3]);
  cosMixingAngle[2][4] = sqrt(1.0 - sinMixingAngle[2][4] * sinMixingAngle[2][4]);

  usingOscParams = true;

}
//________________________________________________________________________________________
void NuOscSterile::SetMixingDoubleAnglesOct2(double sinsq2theta12, double sinsq2theta13, double sinsq2theta14, double sinsq2theta15,
					     double sinsq2theta23, double sinsq2theta24, double sinsq2theta25, double sinsq2theta34, double sinsq2theta35)
{
  //Neglect oscillations between sterile states 4 and 5 since they give no observable effects, i.e. theta45 = 0.

  //Code returns NaN if sin^2(2*mixing angle) is unphysical.
  //Don't let sin^2(2*mixing angle) be below 0 or above 1.

  sinsq2theta12 = TMath::Min(1.,sinsq2theta12);
  sinsq2theta12 = TMath::Max(0.,sinsq2theta12);
  sinsq2theta13 = TMath::Min(1.,sinsq2theta13);
  sinsq2theta13 = TMath::Max(0.,sinsq2theta13);
  sinsq2theta14 = TMath::Min(1.,sinsq2theta14);
  sinsq2theta14 = TMath::Max(0.,sinsq2theta14);
  sinsq2theta15 = TMath::Min(1.,sinsq2theta15);
  sinsq2theta15 = TMath::Max(0.,sinsq2theta15);
  sinsq2theta23 = TMath::Min(1.,sinsq2theta23);
  sinsq2theta23 = TMath::Max(0.,sinsq2theta23);
  sinsq2theta24 = TMath::Min(1.,sinsq2theta24);
  sinsq2theta24 = TMath::Max(0.,sinsq2theta24);
  sinsq2theta25 = TMath::Min(1.,sinsq2theta25);
  sinsq2theta25 = TMath::Max(0.,sinsq2theta25);
  sinsq2theta34 = TMath::Min(1.,sinsq2theta34);
  sinsq2theta34 = TMath::Max(0.,sinsq2theta34);
  sinsq2theta35 = TMath::Min(1.,sinsq2theta35);
  sinsq2theta35 = TMath::Max(0.,sinsq2theta35);

  //Calculate sin and cos(mixingAngle) from sin^2(2*mixingAngle)
  sinMixingAngle[0][1] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta12)));
  sinMixingAngle[0][2] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta13)));
  sinMixingAngle[0][3] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta14)));
  sinMixingAngle[0][4] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta15)));
  //theta23 is in second octant
  sinMixingAngle[1][2] = sqrt(1 - (0.5 * (1 - sqrt(1 - sinsq2theta23))));
  sinMixingAngle[1][3] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta24)));
  sinMixingAngle[1][4] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta25)));
  sinMixingAngle[2][3] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta34)));
  sinMixingAngle[2][4] = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta35)));
  cosMixingAngle[0][1] = sqrt(1.0 - sinMixingAngle[0][1] * sinMixingAngle[0][1]);
  cosMixingAngle[0][2] = sqrt(1.0 - sinMixingAngle[0][2] * sinMixingAngle[0][2]);
  cosMixingAngle[0][3] = sqrt(1.0 - sinMixingAngle[0][3] * sinMixingAngle[0][3]);
  cosMixingAngle[0][4] = sqrt(1.0 - sinMixingAngle[0][4] * sinMixingAngle[0][4]);
  cosMixingAngle[1][2] = sqrt(1.0 - sinMixingAngle[1][2] * sinMixingAngle[1][2]);
  cosMixingAngle[1][3] = sqrt(1.0 - sinMixingAngle[1][3] * sinMixingAngle[1][3]);
  cosMixingAngle[1][4] = sqrt(1.0 - sinMixingAngle[1][4] * sinMixingAngle[1][4]);
  cosMixingAngle[2][3] = sqrt(1.0 - sinMixingAngle[2][3] * sinMixingAngle[2][3]);
  cosMixingAngle[2][4] = sqrt(1.0 - sinMixingAngle[2][4] * sinMixingAngle[2][4]);

  usingOscParams = true;
}
//________________________________________________________________________________________
void NuOscSterile::SetDmsq(double dm21sq, double dm31sq, double dm41sq, double dm51sq)
{
  //This assumes that m4^2 and m5^2 are both > m1^2, m2^2 and m3^2.
  dmsq[1][0] = dm21sq;
  dmsq[2][0] = dm31sq;
  if(fMassHierarchy == kNmhInverted) {
    dmsq[2][0] = -dmsq[2][0];
  }
  dmsq[3][0] = dm41sq;
  dmsq[4][0] = dm51sq;
  dmsq[2][1] = dmsq[2][0] - dmsq[1][0];
  dmsq[3][1] = dm41sq - dm21sq;
  dmsq[4][1] = dm51sq - dm21sq;
  if(fMassHierarchy == kNmhNormal) 
    {
    dmsq[3][2] = dm41sq - dm31sq;
    dmsq[4][2] = dm51sq - dm31sq;
    }
  else if(fMassHierarchy == kNmhInverted)  
    {
    dmsq[3][2] = dm41sq + dm31sq;
    dmsq[4][2] = dm51sq + dm31sq;
    }
  dmsq[4][3] = dm51sq - dm41sq;

}
//________________________________________________________________________________________
void NuOscSterile::SetDeltaCP(double dcp12, double dcp13, double dcp15, double dcp34, double dcp35)
{
  deltaCP[0][1] = dcp12;
  deltaCP[0][2] = dcp13;
  deltaCP[0][4] = dcp15;
  deltaCP[2][3] = dcp34;
  deltaCP[2][4] = dcp35;

  usingOscParams = true;

}
//________________________________________________________________________________________
void NuOscSterile::SetVacuumPMNSReal
(double Ue1, double Ue2, double Ue3, double Ue4, double Ue5,
 double Umu1, double Umu2, double Umu3, double Umu4, double Umu5,
 double Utau1, double Utau2, double Utau3, double Utau4, double Utau5,
 double U41, double U42, double U43, double U44, double U45,
 double U51, double U52, double U53, double U54, double U55
)
{
  UVacuum[0][0][re] = Ue1;
  UVacuum[0][1][re] = Ue2;
  UVacuum[0][2][re] = Ue3;
  UVacuum[0][3][re] = Ue4;
  UVacuum[0][4][re] = Ue5;
  UVacuum[1][0][re] = Umu1;
  UVacuum[1][1][re] = Umu2;
  UVacuum[1][2][re] = Umu3;
  UVacuum[1][3][re] = Umu4;
  UVacuum[1][4][re] = Umu5;
  UVacuum[2][0][re] = Utau1;
  UVacuum[2][1][re] = Utau2;
  UVacuum[2][2][re] = Utau3;
  UVacuum[2][3][re] = Utau4;
  UVacuum[2][4][re] = Utau5;
  UVacuum[3][0][re] = U41;
  UVacuum[3][1][re] = U42;
  UVacuum[3][2][re] = U43;
  UVacuum[3][3][re] = U44;
  UVacuum[3][4][re] = U45;
  UVacuum[4][0][re] = U51;
  UVacuum[4][1][re] = U52;
  UVacuum[4][2][re] = U53;
  UVacuum[4][3][re] = U54;
  UVacuum[4][4][re] = U55;

  usingOscParams = false;

}
//________________________________________________________________________________________

void NuOscSterile::SetVacuumPMNSIm
(double Ue1, double Ue2, double Ue3, double Ue4, double Ue5,
 double Umu1, double Umu2, double Umu3, double Umu4, double Umu5,
 double Utau1, double Utau2, double Utau3, double Utau4, double Utau5,
 double U41, double U42, double U43, double U44, double U45,
 double U51, double U52, double U53, double U54, double U55
)
{
  UVacuum[0][0][im] = Ue1;
  UVacuum[0][1][im] = Ue2;
  UVacuum[0][2][im] = Ue3;
  UVacuum[0][3][im] = Ue4;
  UVacuum[0][4][im] = Ue5;
  UVacuum[1][0][im] = Umu1;
  UVacuum[1][1][im] = Umu2;
  UVacuum[1][2][im] = Umu3;
  UVacuum[1][3][im] = Umu4;
  UVacuum[1][4][im] = Umu5;
  UVacuum[2][0][im] = Utau1;
  UVacuum[2][1][im] = Utau2;
  UVacuum[2][2][im] = Utau3;
  UVacuum[2][3][im] = Utau4;
  UVacuum[2][4][im] = Utau5;
  UVacuum[3][0][im] = U41;
  UVacuum[3][1][im] = U42;
  UVacuum[3][2][im] = U43;
  UVacuum[3][3][im] = U44;
  UVacuum[3][4][im] = U45;
  UVacuum[4][0][im] = U51;
  UVacuum[4][1][im] = U52;
  UVacuum[4][2][im] = U53;
  UVacuum[4][3][im] = U54;
  UVacuum[4][4][im] = U55;

  usingOscParams = false;

}
//________________________________________________________________________________________
void NuOscSterile::SetVacuumPMNSImToZero()
{
  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      UVacuum[i][j][im] = 0.0;
      }
    }

  usingOscParams = false;

}
//________________________________________________________________________________________
void NuOscSterile::SetVacuumPMNSElementByIndex(int iflav, int imass, double real, double imag)
{
  UVacuum[iflav][imass][re] = real;
  UVacuum[iflav][imass][im] = imag;
  usingOscParams = false;
  return;
}
//________________________________________________________________________________________
void NuOscSterile::SetVacuumPMNSElementByType(int nu_pdg, int nmass, double real, double imag)
{
  int imass = nmass-1;
  int iflav;
  switch (abs(nu_pdg)) {
    case kPdgNuE  : iflav=0; break;
    case kPdgNuMu : iflav=1; break;
    case kPdgNuTau: iflav=2; break;
    case kPdgNuS4 : iflav=3; break;
    case kPdgNuS5 : iflav=4; break;
    default: cerr<<"NuOscSterile/SetVacuumPMNSElementByType >> Unrecognised PDG: "<<nu_pdg<<endl;
  }
  UVacuum[iflav][imass][re] = real;
  UVacuum[iflav][imass][im] = imag;
  usingOscParams = false;
  return;
}
//________________________________________________________________________________________ 
void NuOscSterile::CalcVacuumPMNS()
{
  //This calculation for 5 neutrinos is based on appendix A in arXiv 1111.4225.
  //It neglects oscillations between the two sterile states 4 and 5 since they produce no observable effects. 
  
  //Define 2D rotation matrices, where rotationMatrix[i][j] is an identity matrix plus 
  //element[i][i] = cos(mixingAngle[i][j]
  //element[i][j] = sin(mixingAngle[i][j] 
  //element[j][i] = -sin(mixingAngle[i][j] 
  //element[j][j] = cos(mixingAngle[i][j]
  //In matrices 12, 13, 15, 34 and 35, element[i][j] is multplied by exp(-deltaCP[i][j]) and element[j][i] by exp(deltaCP[i][j])
  //Indices 0 and 1 are those of matrix
  //Indices 2 and 3 are elements of that matrix
  //Index 4: 0 is real part, 1 is imaginary part of matrix element
  for (int i=0; i<nFlav; i++)
    {
    for (int j=0; j<nFlav; j++)
      {
      for (int k=0; k<nFlav; k++)
	{
	for (int l=0; l<nFlav; l++)
	  {
	  if(k == l)
	    {
	    rotationMatrix[i][j][k][l][re] = 1.0;
            rotationMatrix[i][j][k][l][im] = 0.0;
	    }
	  if(k == i && l == i)
	    {
	    rotationMatrix[i][j][k][l][re] = cosMixingAngle[i][j];
	    rotationMatrix[i][j][k][l][im] = 0.0; 
	    }
	  else if(k == i && l == j)
	    {
	    if(nFlav == 5)
	      {
	      if((i == 0 && (j == 1 || j == 2 || j == 4)) || (i == 2 && (j == 3 || j == 4)))
	        {
		rotationMatrix[i][j][k][l][re] = sinMixingAngle[i][j] * cos(deltaCP[i][j]);
		rotationMatrix[i][j][k][l][im] = sinMixingAngle[i][j] * -sin(deltaCP[i][j]);
		}
	      else
		{
		rotationMatrix[i][j][k][l][re] = sinMixingAngle[i][j];
		rotationMatrix[i][j][k][l][im] = 0.0;
		}
	      }
	    }
	  else if(k == j && l == i)
	    {
	    if(nFlav == 5)
              {
	      if((i == 0 && (j == 1 || j == 2 || j == 4)) || (i == 2 && (j == 3 || j == 4)))
		{
		rotationMatrix[i][j][k][l][re] = -sinMixingAngle[i][j] * cos(deltaCP[i][j]);
		rotationMatrix[i][j][k][l][im] = -sinMixingAngle[i][j] * sin(deltaCP[i][j]);
		}
	      else
		{
		rotationMatrix[i][j][k][l][re] = -sinMixingAngle[i][j];
		rotationMatrix[i][j][k][l][im] = 0.0;
		}
	      }
	    }
	  else if(k == j && l == j)
	    {
	    rotationMatrix[i][j][k][l][re] = cosMixingAngle[i][j];
	    rotationMatrix[i][j][k][l][im] = 0.0;
	    }
	  else
	    {
	    if(k != l)
	      {
	      rotationMatrix[i][j][k][l][re] = 0.0;
	      rotationMatrix[i][j][k][l][im] = 0.0;
	      }
	    }
	  }
	}
      }
    }

  //Set UVacuum to be an identity matrix.
  for(int i=0; i<nFlav; i++)
    {
      for(int j=0; j<nFlav; j++)
	{
	  if(i == j)
	    UVacuum[i][j][re] = 1.0;
	  else
	    UVacuum[i][j][re] = 0.0;
	  UVacuum[i][j][im] = 0.0;
	}
    }

  double UVacuumCopy[nFlav][nFlav][2];

  //Multiply rotation matrices 12, 13, 14, 15, 23, 24, 25, 34 and 35 in that order.
  //The product of these matrices is the vacuum PMNS matrix UVacuum.
  //In general, rotation matrix 45 would be added to the product, but this is neglected
  //since oscillations from neutrino 4 to 5 or vice versa produce no observable effects. 
  for(int r=0; r<3; r++)
    {
    for(int s=r+1; s<nFlav; s++)
      {
	  
      for(int i=0; i<nFlav; i++)
	{
	for(int j=0; j<nFlav; j++)
	  {
	  UVacuumCopy[i][j][re] = UVacuum[i][j][re];
	  UVacuumCopy[i][j][im] = UVacuum[i][j][im];
          UVacuum[i][j][re] = 0.0;
	  UVacuum[i][j][im] = 0.0;
	  }
	}

      for(int i=0; i<nFlav; i++)
	{
	for(int j=0; j<nFlav; j++)
	  {
	  for(int k=0; k<nFlav; k++)
	    {
            UVacuum[i][j][re] += rotationMatrix[r][s][i][k][re] * UVacuumCopy[k][j][re]
                               - rotationMatrix[r][s][i][k][im] * UVacuumCopy[k][j][im];
            UVacuum[i][j][im] += rotationMatrix[r][s][i][k][re] * UVacuumCopy[k][j][im]
                               + rotationMatrix[r][s][i][k][im] * UVacuumCopy[k][j][re];
	    }
	  }
	}

      }//end of for(int s=r+1; s<nFlav; s++)
    }//end of for(int r=0; r<3; r++)

}
//________________________________________________________________________________________
bool NuOscSterile::Check3FlavPMNS()
{
  //Check the first 3 rows and columns of the PMNS matrix against the standard 3-flavour parametrisation
  //To use this method, all mixing angles with indices 4 and 5 must be set to zero

  double PMNSMatrix3Flav[3][3][2];

  PMNSMatrix3Flav[0][0][re] = cosMixingAngle[0][1]*cosMixingAngle[0][2];
  PMNSMatrix3Flav[0][0][im] = 0.0;
  PMNSMatrix3Flav[0][1][re] = cosMixingAngle[0][2]*sinMixingAngle[0][1];
  PMNSMatrix3Flav[0][1][im] = 0.0;
  PMNSMatrix3Flav[0][2][re] = sinMixingAngle[0][2]*cos(deltaCP[0][2]);
  PMNSMatrix3Flav[0][2][im] = sinMixingAngle[0][2]*-sin(deltaCP[0][2]);
  PMNSMatrix3Flav[1][0][re] = -cosMixingAngle[1][2]*sinMixingAngle[0][1]-sinMixingAngle[0][2]*sinMixingAngle[1][2]*cosMixingAngle[0][1]*cos(deltaCP[0][2]);
  PMNSMatrix3Flav[1][0][im] = -sinMixingAngle[0][2]*sinMixingAngle[1][2]*cosMixingAngle[0][1]*sin(deltaCP[0][2]);
  PMNSMatrix3Flav[1][1][re] = cosMixingAngle[0][1]*cosMixingAngle[1][2]-sinMixingAngle[0][1]*sinMixingAngle[0][2]*sinMixingAngle[1][2]*cos(deltaCP[0][2]);
  PMNSMatrix3Flav[1][1][im] = -sinMixingAngle[0][1]*sinMixingAngle[0][2]*sinMixingAngle[1][2]*sin(deltaCP[0][2]);
  PMNSMatrix3Flav[1][2][re] = cosMixingAngle[0][2]*sinMixingAngle[1][2];
  PMNSMatrix3Flav[1][2][im] = 0.0;
  PMNSMatrix3Flav[2][0][re] = sinMixingAngle[0][1]*sinMixingAngle[1][2]-sinMixingAngle[0][2]*cosMixingAngle[0][1]*cosMixingAngle[1][2]*cos(deltaCP[0][2]);
  PMNSMatrix3Flav[2][0][im] = -sinMixingAngle[0][2]*cosMixingAngle[0][1]*cosMixingAngle[1][2]*sin(deltaCP[0][2]);
  PMNSMatrix3Flav[2][1][re] = -sinMixingAngle[1][2]*cosMixingAngle[0][1]-sinMixingAngle[0][1]*cosMixingAngle[1][2]*sinMixingAngle[0][2]*cos(deltaCP[0][2]);
  PMNSMatrix3Flav[2][1][im] = -sinMixingAngle[0][1]*cosMixingAngle[1][2]*sinMixingAngle[0][2]*sin(deltaCP[0][2]);
  PMNSMatrix3Flav[2][2][re] = cosMixingAngle[0][2]*cosMixingAngle[1][2];
  PMNSMatrix3Flav[2][2][im] = 0.0;

  this->CalcVacuumPMNS();

  for(int i=0; i<3; i++)
    {
    for(int j=0; j<3; j++)
      {
      if(fabs(UVacuum[i][j][re] - PMNSMatrix3Flav[i][j][re]) > 1E-12)
	{
	cout<<"NuOscSterile::Check3FlavPMNS() check failed for PMNS matrix element "<<i<<" "<<j<<" re; ";
	cout<<"difference between UVacuum and PMNSMatrix3Flav is "<<fabs(UVacuum[i][j][re] - PMNSMatrix3Flav[i][j][re])<<endl; 
        return false;
	}
      if(fabs(UVacuum[i][j][im] - PMNSMatrix3Flav[i][j][im]) > 1E-12)
        {
        cout<<"NuOscSterile::Check3FlavPMNS() check failed for PMNS matrix element "<<i<<" "<<j<<" im; ";
        cout<<"difference between UVacuum and PMNSMatrix3Flav is "<<fabs(UVacuum[i][j][im] - PMNSMatrix3Flav[i][j][im])<<endl;
        return false;
        }
      }
    }

  return true;
}
//________________________________________________________________________________________
void NuOscSterile::CalcMEProduct(int nu_start, int nu_end, int mass1, int mass2)
{

  double MEProduct12[2];
  double MEProduct34[2];

  MEProduct12[re] = UVacuum[nu_start][mass2][re]  *  UVacuum[nu_end][mass2][re]
                - (-UVacuum[nu_start][mass2][im]) *  UVacuum[nu_end][mass2][im];
  MEProduct12[im] = UVacuum[nu_start][mass2][re]  *  UVacuum[nu_end][mass2][im]
                + (-UVacuum[nu_start][mass2][im]) *  UVacuum[nu_end][mass2][re];

  MEProduct34[re] = UVacuum[nu_start][mass1][re]  *  UVacuum[nu_end][mass1][re]
                 - (UVacuum[nu_start][mass1][im]) * -UVacuum[nu_end][mass1][im];
  MEProduct34[im] = UVacuum[nu_start][mass1][re]  * -UVacuum[nu_end][mass1][im]
                 + (UVacuum[nu_start][mass1][im]) *  UVacuum[nu_end][mass1][re];


  MEProduct[re] = MEProduct12[re] * MEProduct34[re] - MEProduct12[im] * MEProduct34[im];
  MEProduct[im] = MEProduct12[re] * MEProduct34[im] + MEProduct12[im] * MEProduct34[re]; 

}
//________________________________________________________________________________________
double NuOscSterile::CalcVacuumOscProb(int nu_start, int nu_end, double energy)
{
  //This method calculates the oscillation probabilities in vacuum by extending the 
  //standard 3-flavour formula to 4 or 5 flavours.
  //Oscillations between the sterile states 4 and 5 are neglected as they give no observable effects.
  //The 3-flavour formula is found in many papers and theses, e.g. equation 11 in hep-ph 0506165.
  //The neutrino energy is in GeV.

  OscChannelId_t osch = util::GetOscChannel(nu_start, nu_end);

  if(osch == kOscUndefined) {
    cerr << "NuOscSterile: Not allowing (neutrino = " << nu_start << ") --> "
         << "(neutrino = " << nu_end << ") oscillations" << endl;
    return 0.;
  }

  double oscProb = 0.0;

  if(usingOscParams)
    this->CalcVacuumPMNS();

  //**************************** OSCILLATIONS BETWEEN ACTIVE NEUTRINOS ******************************

  //nue->nue, nuebar->nuebar
  if(osch == kOscNuE2NuE || osch == kOscNuEBar2NuEBar)
    {
    oscProb = 1.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
	this->CalcMEProduct(0, 0, i, j);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    } 

  //nue->numu, numubar->nuebar
  else if(osch == kOscNuE2NuMu || osch == kOscNuMuBar2NuEBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(0, 1, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nue->nutau, nutaubar->nuebar
  else if(osch == kOscNuE2NuTau || osch == kOscNuTauBar2NuEBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(0, 2, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //numu->nue, nuebar->numubar 
  else if(osch == kOscNuMu2NuE || osch == kOscNuEBar2NuMuBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(0, 1, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);	
        }
      } 
    }

  //numu->numu, numubar->numubar
  else if(osch == kOscNuMu2NuMu || osch == kOscNuMuBar2NuMuBar)
    {
    oscProb = 1.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(1, 1, i, j);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
	        * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //numu->nutau, nutaubar->numubar
  else if(osch == kOscNuMu2NuTau || osch == kOscNuTauBar2NuMuBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(1, 2, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nutau->nue, nuebar->nutaubar
  else if(osch == kOscNuTau2NuE || osch == kOscNuEBar2NuTauBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(0, 2, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nutau->numu, numubar->nutaubar
  else if(osch == kOscNuTau2NuMu || osch == kOscNuMuBar2NuTauBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(1, 2, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
       }
     }
   }

  //nutau->nutau, nutaubar->nutaubar
  else if(osch == kOscNuTau2NuTau || osch == kOscNuTauBar2NuTauBar)
    {
    oscProb = 1.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
	this->CalcMEProduct(2, 2, i, j);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
	         * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //********************** OSCILLATIONS BETWEEN ACTIVE AND STERILE NEUTRINOS ***************************

  //nue->nus4, nus4bar->nuebar
  else if(osch ==  kOscNuE2NuS4 || osch == kOscNuS4Bar2NuEBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(0, 3, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nue->nus5, nus5bar->nuebar
  else if(osch ==  kOscNuE2NuS5 || osch == kOscNuS5Bar2NuEBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
	this->CalcMEProduct(0, 4, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //numu->nus4, nus4bar->numubar
  else if(osch ==  kOscNuMu2NuS4 || osch == kOscNuS4Bar2NuMuBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(1, 3, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //numu->nus5, nus5bar->numubar
  else if(osch ==  kOscNuMu2NuS5 || osch == kOscNuS5Bar2NuMuBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
	this->CalcMEProduct(1, 4, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nutau->nus4, nus4bar->nutaubar
  else if(osch ==  kOscNuTau2NuS4 || osch == kOscNuS4Bar2NuTauBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
        this->CalcMEProduct(2, 3, i, j);
        oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        }
      }
    }

  //nutau->nus5, nus5bar->nutaubar
  else if(osch ==  kOscNuTau2NuS5 || osch == kOscNuS5Bar2NuTauBar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
	this->CalcMEProduct(2, 4, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nus4->nue, nuebar->nus4bar
  else if(osch ==  kOscNuS42NuE || osch == kOscNuEBar2NuS4Bar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
        this->CalcMEProduct(0, 3, i, j);
        oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        }
      }
    }

  //nus4->numu, numubar->nus4bar
  else if(osch ==  kOscNuS42NuMu || osch == kOscNuMuBar2NuS4Bar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
        this->CalcMEProduct(1, 3, i, j);
        oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        }
      }
    }

  //nus4->nutau, nutaubar->nus4bar
  else if(osch ==  kOscNuS42NuTau || osch == kOscNuTauBar2NuS4Bar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(2, 3, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  //nus5->nue, nuebar->nus5bar
  else if(osch ==  kOscNuS52NuE || osch == kOscNuEBar2NuS5Bar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
        this->CalcMEProduct(0, 4, i, j);
        oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        }
      }
    }

  //nus5->numu, numubar->nus5bar
  else if(osch ==  kOscNuS52NuMu || osch == kOscNuMuBar2NuS5Bar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
        {
        this->CalcMEProduct(1, 4, i, j);
        oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
        }
      }
    }

  //nus5->nutau, nutaubar->nus5bar
  else if(osch ==  kOscNuS52NuTau || osch == kOscNuTauBar2NuS5Bar)
    {
    oscProb = 0.0;
    for(int i=0; i<nFlav-1; i++)
      {
      for(int j=i+1; j<nFlav; j++)
	{
	this->CalcMEProduct(2, 4, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
                * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  return oscProb;
  
}
//________________________________________________________________________________________
double NuOscSterile::CalcMatterOscProb(int nu_start, int nu_end, double energy)
{
  //This method calculates the oscillation probabilities in constant density matter.
  //The neutrino energy is in GeV.

  OscChannelId_t osch = util::GetOscChannel(nu_start, nu_end);

  if(osch == kOscUndefined) {
    cerr << "NuOscSterile: Not allowing (neutrino = " << nu_start << ") --> "
         << "(neutrino = " << nu_end << ") oscillations" << endl;
    return 0.;
  }

  double oscProb = 0.0;

  if(usingOscParams)
    this->CalcVacuumPMNS();

  //UMassUdagger is the product of the PMNS matrix, the mass squared matrix and the
  //dagger of the PMNS matrix.
  double UMassUdagger[nFlav][nFlav][2] = {};

  //The mass squared matrix has diagonal elements 0, dm21sq, dm31sq, dm41sq, dm51sq,
  //and all off-diagonal elements are zero.
  double massSqMatrix[nFlav][nFlav][2] = {};

  //UMass is the product of UVacuum and the mass squared matrix
  double UMass[nFlav][nFlav][2] = {};

  massSqMatrix[1][1][re] = dmsq[1][0];
  massSqMatrix[2][2][re] = dmsq[2][0];
  massSqMatrix[3][3][re] = dmsq[3][0];
  massSqMatrix[4][4][re] = dmsq[4][0];
  
  for(int i=0; i<nFlav; i++)
    {
    for(int j=0; j<nFlav; j++)
      {
      for(int k=0; k<nFlav; k++)
        {
        UMass[i][j][re] += UVacuum[i][k][re] * massSqMatrix[k][j][re]
	                 - UVacuum[i][k][im] * massSqMatrix[k][j][im];
        UMass[i][j][im] += UVacuum[i][k][re] * massSqMatrix[k][j][im]
	                 + UVacuum[i][k][im] * massSqMatrix[k][j][re];
        }
      }
    }

  for(int i=0; i<nFlav; i++)
    {
    for(int j=0; j<nFlav; j++)
      {
      for(int k=0; k<nFlav; k++)
	{
        UMassUdagger[i][j][re] += UMass[i][k][re] *  UVacuum[j][k][re]
	                        - UMass[i][k][im] * -UVacuum[j][k][im];
        UMassUdagger[i][j][im] += UMass[i][k][re] * -UVacuum[j][k][im]
	                        + UMass[i][k][im] *  UVacuum[j][k][re];
	}
      }
    }

  //For antineutrinos, use the complex conjugate of UMassUdagger
  if(util::IsAntiNeutrino(nu_start))
    {
    for(int i=0; i<nFlav; i++)
      {
      for(int j=0; j<nFlav; j++)
        {
	UMassUdagger[i][j][im] = -UMassUdagger[i][j][im];
        }
      }
    }

  //Define number density of neutrons in matter and assume it is the same as that of electrons
  double Nn = fNe;

  // Add potential due to constant density matter to UMassUdagger.
  // The factor of 1.0e9 is to convert the energy from GeV to eV.
  if(util::IsNeutrino(nu_start))
    {
    //First term is for CC interactions between nue and electrons in matter 
    //Second term is for NC interactions between any active neutrino and protons, neutrons or electrons in matter
    UMassUdagger[0][0][re] = UMassUdagger[0][0][re] + 2.0*1.0e9*energy*sqrt(2.0)*kGF_SI*fNe - 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    //Only apply potential for NC interactions since only nue can have CC interactions with electrons in matter
    UMassUdagger[1][1][re] -= 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    UMassUdagger[2][2][re] -= 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    }
  else if(util::IsAntiNeutrino(nu_start))
    {
    UMassUdagger[0][0][re] = UMassUdagger[0][0][re] - 2.0*1.0e9*energy*sqrt(2.0)*kGF_SI*fNe + 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    UMassUdagger[1][1][re] += 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    UMassUdagger[2][2][re] += 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    }
  
  double eigenvalues[nFlav] = {};
  double UMatter[nFlav][nFlav][2] = {};
  double UMatterDagger[nFlav][nFlav][2] = {};
  
  //Diagonalise UMassUdagger, put eigenvalues into array with same name and enter eigenvectors as  
  //columns of UMatter. 
  //Use complex form of symmetric bidiagonalisation and QR reduction method in GSL.
#ifdef __GSL_ENABLED__ 
  util::DiagonaliseGSL(nFlav, &UMassUdagger[0][0][0], &UMatter[0][0][0], &eigenvalues[0]);
#else
  std::cerr<<"Need GSL"<<std::endl;
  exit(-1);
#endif

#ifdef __check_diagonalisation__
  // Check that the sum of the eigenvalues is equal to the trace of the matrix
  double eigenvalueSum = eigenvalues[0] + eigenvalues[1] + eigenvalues[2] + eigenvalues[3] + eigenvalues[4];
  double trace = UMassUdagger[0][0][re] + UMassUdagger[1][1][re] + UMassUdagger[2][2][re] + UMassUdagger[3][3][re] + UMassUdagger[4][4][re];

  cout<<"eigenvalueSum / trace =   "<<eigenvalueSum<<"    "<<trace<<endl;
  cout<<"Fractional difference =   "<<fabs((eigenvalueSum - trace) / trace)<<endl; 

  if(fabs((eigenvalueSum - trace) / trace) > 1.0e-10) 
    {
      cerr << "Warning in NuOscSterile::CalcMatterOscProb() : fractional difference between sum of eigenvalues and trace of matrix = "
	   << fabs((eigenvalueSum - trace) / trace) << "." << endl;
    }

  double UMassUdaggerMinusEv[nFlav][nFlav][2];
  for(int i=0; i<nFlav; i++)
    {
    for(int j=0; j<nFlav; j++)
      {
      for(int k=0; k<2; k++)
	{
	UMassUdaggerMinusEv[i][j][k] = UMassUdagger[i][j][k];
	}//k
      }//j
    }//i

  double eigenvectorCheck[nFlav][2] = {};

  for(int i=0; i<nFlav; i++)
    {
    // UMassUdaggerMinusEv is UMassUdagger with the eigenvalue subtracted from all diagonal real elements
    for(int j=0; j<nFlav; j++) 
      UMassUdaggerMinusEv[j][j][re] = UMassUdagger[j][j][re] - eigenvalues[i];
    
    // Check the eigenvector is correct.
    // Multiplying any row of UMassUdaggerMinusEv by the column eigenvector should give zero   
    for(int m=0; m<nFlav; m++)
      {
      eigenvectorCheck[i][re] += UMassUdaggerMinusEv[i][m][re] * UMatter[m][i][re]
	                       - UMassUdaggerMinusEv[i][m][im] * UMatter[m][i][im];
      eigenvectorCheck[i][im] += UMassUdaggerMinusEv[i][m][re] * UMatter[m][i][im]
	                       + UMassUdaggerMinusEv[i][m][im] * UMatter[m][i][re];
      }

    }//end of for(int i=0; i<nFlav; i++)

  for(int i=0; i<nFlav; i++) 
    {
      cout<<"eigenvectorCheck   "<<eigenvectorCheck[i][re]<<"    "<<eigenvectorCheck[i][im]<<endl;
      if(fabs(eigenvectorCheck[i][re]) > 1.0e-10) 
	{
	  cerr << "Warning in NuOscSterile::CalcMatterOscProb() : real part of eigenvector check = "
	       << eigenvectorCheck[i][re] <<"." << endl;
	}
      if(fabs(eigenvectorCheck[i][im]) > 1.0e-10) 
	{
	  cerr << "Warning in NuOscSterile::CalcMatterOscProb() : imaginary part of eigenvector check = "
	       << eigenvectorCheck[i][im] << "." << endl;
	}
    }//end of for(int i=0; i<nFlav; i++)

#endif //__check_diagonalisation__

  //UMatterDagger is the Hermitian conjugate (complex conjugate of transpose) of UMatter
  for(int i=0; i<nFlav; i++)
    {
    for(int j=0; j<nFlav; j++)
      {
      UMatterDagger[i][j][re] =  UMatter[j][i][re];
      UMatterDagger[i][j][im] = -UMatter[j][i][im];
      }
    }

  double initialFlavourState[nFlav][2];
  for(int i=0; i<nFlav; i++)
    {
    initialFlavourState[i][re] = 0.0;
    initialFlavourState[i][im] = 0.0;
    }

  // The neutrino is formed in a definite flavour state
  if(nu_start == kPdgNuE || nu_start == kPdgNuEBar)
    initialFlavourState[0][re] = 1.0;
  else if(nu_start == kPdgNuMu || nu_start == kPdgNuMuBar)
    initialFlavourState[1][re] = 1.0;
  else if(nu_start == kPdgNuTau || nu_start == kPdgNuTauBar)
    initialFlavourState[2][re] = 1.0;
  else if(nu_start == kPdgNuS4 || nu_start == kPdgNuS4Bar)
    initialFlavourState[3][re] = 1.0;
  else if(nu_start == kPdgNuS5 || nu_start == kPdgNuS5Bar)
    initialFlavourState[4][re] = 1.0;

  double initialMassState[nFlav][2];
  for(int i=0; i<nFlav; i++)
  {
    initialMassState[i][re] = 0.0;
    initialMassState[i][im] = 0.0;
  }

  // Convert the initial flavour state to initial mass states by multiplying by UMatterDagger
  for(int i=0; i<nFlav; i++)
  {
    for(int j=0; j<nFlav; j++)
    {
      initialMassState[i][re] += UMatterDagger[i][j][re] * initialFlavourState[j][re]
                               - UMatterDagger[i][j][im] * initialFlavourState[j][im];
      initialMassState[i][im] += UMatterDagger[i][j][re] * initialFlavourState[j][im]
                               + UMatterDagger[i][j][im] * initialFlavourState[j][re];
    }
  }

  // Now propagate the mass eigenstates to the far detector
  double finalMassState[nFlav][2];
  for(int i=0; i<nFlav; i++)
  {
    finalMassState[i][re] = 0.0;
    finalMassState[i][im] = 0.0;
  }

  // For the first mass eigenstate, the difference in the squared masses is zero.
  // Therefore only retain terms multiplied by cos(0), and drop those multiplied by sin(0).
  finalMassState[0][re] = initialMassState[0][re];
  finalMassState[0][im] = initialMassState[0][im];

  for(int i=1; i<nFlav; i++)
  {
    finalMassState[i][re] = initialMassState[i][re] *  cos((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy)
                          - initialMassState[i][im] * -sin((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy);  
    finalMassState[i][im] = initialMassState[i][re] * -sin((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy)
                          + initialMassState[i][im] *  cos((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy);
  }

  double finalFlavourState[nFlav][2];
  for(int i=0; i<nFlav; i++)
  {
    finalFlavourState[i][re] = 0.0;
    finalFlavourState[i][im] = 0.0;
  }

  // Convert final mass states to flavour states by multiplying by UMatter 
  for(int i=0; i<nFlav; i++)
  {
    for(int j=0; j<nFlav; j++)
    {
      finalFlavourState[i][re] += UMatter[i][j][re] * finalMassState[j][re]
				- UMatter[i][j][im] * finalMassState[j][im];
      finalFlavourState[i][im] += UMatter[i][j][re] * finalMassState[j][im]
				+ UMatter[i][j][im] * finalMassState[j][re];
    }
  }

  //Each component of finalFlavourState is the amplitude for the neutrino to interact 
  //in the far detector as that flavour.
  //The oscillation probabilities are the moduli squared of the amplitudes
  if(nu_end == kPdgNuE || nu_end == kPdgNuEBar)
    {
    oscProb = finalFlavourState[0][re]*finalFlavourState[0][re] 
            + finalFlavourState[0][im]*finalFlavourState[0][im];
    }
  else if(nu_end == kPdgNuMu || nu_end == kPdgNuMuBar)
    {
    oscProb = finalFlavourState[1][re]*finalFlavourState[1][re] 
            + finalFlavourState[1][im]*finalFlavourState[1][im];
    }
  else if(nu_end == kPdgNuTau || nu_end == kPdgNuTauBar)
    {
    oscProb = finalFlavourState[2][re]*finalFlavourState[2][re] 
            + finalFlavourState[2][im]*finalFlavourState[2][im];
    }
  else if(nu_end == kPdgNuS4 || nu_end == kPdgNuS4Bar)
    {
    oscProb = finalFlavourState[3][re]*finalFlavourState[3][re]
	    + finalFlavourState[3][im]*finalFlavourState[3][im];
    }
  else if(nu_end == kPdgNuS5|| nu_end == kPdgNuS5Bar)
    {
    oscProb = finalFlavourState[4][re]*finalFlavourState[4][re]
	    + finalFlavourState[4][im]*finalFlavourState[4][im];
    }

  return oscProb;

}
//________________________________________________________________________________________
void NuOscSterile::CalcMatterOscProbAll(bool isAntineutrino, double energy, double * osc_probs_out)
{
  //This method calculates all the oscillation probabilities in constant density matter.
  //The neutrino energy is in GeV.

  if(usingOscParams)
    this->CalcVacuumPMNS();

  //UMassUdagger is the product of the PMNS matrix, the mass squared matrix and the
  //dagger of the PMNS matrix.
  double UMassUdagger[nFlav][nFlav][2] = {};

  //The mass squared matrix has diagonal elements 0, dm21sq, dm31sq, dm41sq, dm51sq,
  //and all off-diagonal elements are zero.
  double massSqMatrix[nFlav][nFlav][2] = {};

  //UMass is the product of UVacuum and the mass squared matrix
  double UMass[nFlav][nFlav][2] = {};

  massSqMatrix[1][1][re] = dmsq[1][0];
  massSqMatrix[2][2][re] = dmsq[2][0];
  massSqMatrix[3][3][re] = dmsq[3][0];
  massSqMatrix[4][4][re] = dmsq[4][0];
  
  for(int i=0; i<nFlav; i++) {
    for(int j=0; j<nFlav; j++) {
      for(int k=0; k<nFlav; k++) {
        UMass[i][j][re] += UVacuum[i][k][re] * massSqMatrix[k][j][re]
	                 - UVacuum[i][k][im] * massSqMatrix[k][j][im];
        UMass[i][j][im] += UVacuum[i][k][re] * massSqMatrix[k][j][im]
	                 + UVacuum[i][k][im] * massSqMatrix[k][j][re];
      }
    }
  }

  for(int i=0; i<nFlav; i++) {
    for(int j=0; j<nFlav; j++) {
      for(int k=0; k<nFlav; k++) {
        UMassUdagger[i][j][re] += UMass[i][k][re] *  UVacuum[j][k][re]
	                        - UMass[i][k][im] * -UVacuum[j][k][im];
        UMassUdagger[i][j][im] += UMass[i][k][re] * -UVacuum[j][k][im]
	                        + UMass[i][k][im] *  UVacuum[j][k][re];
      }
    }
  }

  //For antineutrinos, use the complex conjugate of UMassUdagger
  if(isAntineutrino) {
    for(int i=0; i<nFlav; i++) {
      for(int j=0; j<nFlav; j++) {
	UMassUdagger[i][j][im] = -UMassUdagger[i][j][im];
      }
    }
  }

  //Define number density of neutrons in matter and assume it is the same as that of electrons
  double Nn = fNe;

  // Add potential due to constant density matter to UMassUdagger.
  // The factor of 1.0e9 is to convert the energy from GeV to eV.
  if(!isAntineutrino) {
    //First term is for CC interactions between nue and electrons in matter 
    //Second term is for NC interactions between any active neutrino and protons, neutrons or electrons in matter
    UMassUdagger[0][0][re] = UMassUdagger[0][0][re] + 2.0*1.0e9*energy*sqrt(2.0)*kGF_SI*fNe - 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    //Only apply potential for NC interactions since only nue can have CC interactions with electrons in matter
    UMassUdagger[1][1][re] -= 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    UMassUdagger[2][2][re] -= 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
  }
  else {
    UMassUdagger[0][0][re] = UMassUdagger[0][0][re] - 2.0*1.0e9*energy*sqrt(2.0)*kGF_SI*fNe + 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    UMassUdagger[1][1][re] += 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
    UMassUdagger[2][2][re] += 2.0*1.0e9*energy*0.5*sqrt(2.0)*kGF_SI*Nn;
  }
  
  double eigenvalues[nFlav] = {};
  double UMatter[nFlav][nFlav][2] = {};
  double UMatterDagger[nFlav][nFlav][2] = {};
  
  //Diagonalise UMassUdagger, put eigenvalues into array with same name and enter eigenvectors as  
  //columns of UMatter. 
  //Use complex form of symmetric bidiagonalisation and QR reduction method in GSL.
#ifdef __GSL_ENABLED__ 
  util::DiagonaliseGSL(nFlav, &UMassUdagger[0][0][0], &UMatter[0][0][0], &eigenvalues[0]);
#else
  std::cerr<<"Need GSL"<<std::endl;
  exit(-1);
#endif

#ifdef __check_diagonalisation__
  // Check that the sum of the eigenvalues is equal to the trace of the matrix
  double eigenvalueSum = eigenvalues[0] + eigenvalues[1] + eigenvalues[2] + eigenvalues[3] + eigenvalues[4];
  double trace = UMassUdagger[0][0][re] + UMassUdagger[1][1][re] + UMassUdagger[2][2][re] + UMassUdagger[3][3][re] + UMassUdagger[4][4][re];

  cout<<"eigenvalueSum / trace =   "<<eigenvalueSum<<"    "<<trace<<endl;
  cout<<"Fractional difference =   "<<fabs((eigenvalueSum - trace) / trace)<<endl; 

  if(fabs((eigenvalueSum - trace) / trace) > 1.0e-10) {
    cerr << "Warning in NuOscSterile::CalcMatterOscProb() : fractional difference between sum of eigenvalues and trace of matrix = "
	 << fabs((eigenvalueSum - trace) / trace) << "." << endl;
  }

  double UMassUdaggerMinusEv[nFlav][nFlav][2];
  for(int i=0; i<nFlav; i++) {
    for(int j=0; j<nFlav; j++) {
      for(int k=0; k<2; k++) {
	UMassUdaggerMinusEv[i][j][k] = UMassUdagger[i][j][k];
      }//k
    }//j
  }//i

  double eigenvectorCheck[nFlav][2] = {};

  for(int i=0; i<nFlav; i++) {
    // UMassUdaggerMinusEv is UMassUdagger with the eigenvalue subtracted from all diagonal real elements
    for(int j=0; j<nFlav; j++) 
      UMassUdaggerMinusEv[j][j][re] = UMassUdagger[j][j][re] - eigenvalues[i];
    
    // Check the eigenvector is correct.
    // Multiplying any row of UMassUdaggerMinusEv by the column eigenvector should give zero   
    for(int m=0; m<nFlav; m++) {
      eigenvectorCheck[i][re] += UMassUdaggerMinusEv[i][m][re] * UMatter[m][i][re]
	                       - UMassUdaggerMinusEv[i][m][im] * UMatter[m][i][im];
      eigenvectorCheck[i][im] += UMassUdaggerMinusEv[i][m][re] * UMatter[m][i][im]
	                       + UMassUdaggerMinusEv[i][m][im] * UMatter[m][i][re];
    }
  }//end of for(int i=0; i<nFlav; i++)

  for(int i=0; i<nFlav; i++) {
    cout<<"eigenvectorCheck   "<<eigenvectorCheck[i][re]<<"    "<<eigenvectorCheck[i][im]<<endl;
    if(fabs(eigenvectorCheck[i][re]) > 1.0e-10) {
      cerr << "Warning in NuOscSterile::CalcMatterOscProb() : real part of eigenvector check = "
	   << eigenvectorCheck[i][re] <<"." << endl;
    }
    if(fabs(eigenvectorCheck[i][im]) > 1.0e-10) {
      cerr << "Warning in NuOscSterile::CalcMatterOscProb() : imaginary part of eigenvector check = "
	   << eigenvectorCheck[i][im] << "." << endl;
    }
  }//end of for(int i=0; i<nFlav; i++)

#endif //__check_diagonalisation__

  //UMatterDagger is the Hermitian conjugate (complex conjugate of transpose) of UMatter
  for(int i=0; i<nFlav; i++) {
    for(int j=0; j<nFlav; j++) {
      UMatterDagger[i][j][re] =  UMatter[j][i][re];
      UMatterDagger[i][j][im] = -UMatter[j][i][im];
    }
  }
  
  // Treat that block of memory as a 2D array
  double(&osc_probs)[nFlav][nFlav] = *reinterpret_cast< double (*) [nFlav][nFlav]> (osc_probs_out);
  
  for (int i_initflav = 0 ; i_initflav < nFlav ; i_initflav++) {
    double initialFlavourState[nFlav][2] = {};
    double initialMassState   [nFlav][2] = {};
    double finalMassState     [nFlav][2] = {};
    double finalFlavourState  [nFlav][2] = {};
    
    initialFlavourState[i_initflav][re] = 1.0;

    // Convert the initial flavour state to initial mass states by multiplying by UMatterDagger
    for(int i=0; i<nFlav; i++) {
      for(int j=0; j<nFlav; j++) {
	initialMassState[i][re] += UMatterDagger[i][j][re] * initialFlavourState[j][re]
				 - UMatterDagger[i][j][im] * initialFlavourState[j][im];
	initialMassState[i][im] += UMatterDagger[i][j][re] * initialFlavourState[j][im]
				 + UMatterDagger[i][j][im] * initialFlavourState[j][re];
      }
    }
    
    // Now propagate the mass eigenstates to the far detector
    // For the first mass eigenstate, the difference in the squared masses is zero.
    // Therefore only retain terms multiplied by cos(0), and drop those multiplied by sin(0).
    finalMassState[0][re] = initialMassState[0][re];
    finalMassState[0][im] = initialMassState[0][im];

    for(int i=1; i<nFlav; i++) {
      finalMassState[i][re] = initialMassState[i][re] *  cos((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy)
			    - initialMassState[i][im] * -sin((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy);  
      finalMassState[i][im] = initialMassState[i][re] * -sin((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy)
			    + initialMassState[i][im] *  cos((2.0 * kOscDimensFactor * (eigenvalues[i] - eigenvalues[0]) * fBaseline) / energy);
    }

    // Convert final mass states to flavour states by multiplying by UMatter 
    for(int i=0; i<nFlav; i++) {
      for(int j=0; j<nFlav; j++) {
	finalFlavourState[i][re] += UMatter[i][j][re] * finalMassState[j][re]
				  - UMatter[i][j][im] * finalMassState[j][im];
	finalFlavourState[i][im] += UMatter[i][j][re] * finalMassState[j][im]
				  + UMatter[i][j][im] * finalMassState[j][re];
      }
    }

    //Each component of finalFlavourState is the amplitude for the neutrino to interact 
    //in the far detector as that flavour.
    //The oscillation probabilities are the moduli squared of the amplitudes
    for (int i_finalflav = 0 ; i_finalflav < nFlav ; i_finalflav++ ) {
      osc_probs[i_initflav][i_finalflav] = finalFlavourState[i_finalflav][re]*finalFlavourState[i_finalflav][re] 
					 + finalFlavourState[i_finalflav][im]*finalFlavourState[i_finalflav][im];
    }
    
  } //i_initflav
  
  return;
}
