//________________________________________________________________________________________
/*

 >> ChangeLog:


*/
//________________________________________________________________________________________

#include <TMath.h>

#include "NuOscSterile.h"
#include "Conventions.h"
#include "Utilities.h"
#include "Opt.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

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

  for (int i=0; i<nFlav+1; i++)
    {
    for (int j=0; j<nFlav+1; j++)
      {
      sinMixingAngle[i][j] = 0.0;
      cosMixingAngle[i][j] = 1.0;
      }
    }

  for (int i=0; i<nFlav+1; i++)
    {
    for (int j=0; j<nFlav+1; j++)
      {
      dmsq[i][j] = 0.0;
      }
    }  

  for (int i=0; i<nFlav+1; i++)
    {
    for (int j=0; j<nFlav+1; j++)
      {
      deltaCP[i][j] = 0.0;
      }
    }

  for (int i=0; i<nFlav+1; i++)
    {
    for (int j=0; j<nFlav+1; j++)
      {
      for (int k=0; k<nFlav+1; k++)
	{
	for (int l=0; l<nFlav+1; l++)
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

  for (int i=0; i<nFlav+1; i++)
    {
    for (int j=0; j<nFlav+1; j++)
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
  cout << "NuOscSterile >> Mass hierarchy assumed: " << util::AsString(h) << endl;

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
  sinMixingAngle[1][2] = sqrt(sinsqtheta12);
  sinMixingAngle[1][3] = sqrt(sinsqtheta13);
  sinMixingAngle[1][4] = sqrt(sinsqtheta14);
  sinMixingAngle[1][5] = sqrt(sinsqtheta15);
  sinMixingAngle[2][3] = sqrt(sinsqtheta23);
  sinMixingAngle[2][4] = sqrt(sinsqtheta24);
  sinMixingAngle[2][5] = sqrt(sinsqtheta25);
  sinMixingAngle[3][4] = sqrt(sinsqtheta34);
  sinMixingAngle[3][5] = sqrt(sinsqtheta35);
  cosMixingAngle[1][2] = sqrt(1.0 - sinsqtheta12);
  cosMixingAngle[1][3] = sqrt(1.0 - sinsqtheta13);
  cosMixingAngle[1][4] = sqrt(1.0 - sinsqtheta14);
  cosMixingAngle[1][5] = sqrt(1.0 - sinsqtheta15);
  cosMixingAngle[2][3] = sqrt(1.0 - sinsqtheta23);
  cosMixingAngle[2][4] = sqrt(1.0 - sinsqtheta24);
  cosMixingAngle[2][5] = sqrt(1.0 - sinsqtheta25);
  cosMixingAngle[3][4] = sqrt(1.0 - sinsqtheta34);
  cosMixingAngle[3][5] = sqrt(1.0 - sinsqtheta35);

  usingOscParams = true;

}
//________________________________________________________________________________________
void NuOscSterile::SetDmsq(double dm21sq, double dm31sq, double dm41sq, double dm51sq)
{
  //This assumes that m4^2 and m5^2 are both > m1^2, m2^2 and m3^2.
  dmsq[2][1] = dm21sq;
  dmsq[3][1] = dm31sq;
  if(fMassHierarchy == kNmhInverted) {
    dmsq[3][1] = -dmsq[3][1];
  }
  dmsq[4][1] = dm41sq;
  dmsq[5][1] = dm51sq;
  dmsq[3][2] = dmsq[3][1] - dmsq[2][1];
  dmsq[4][2] = dm41sq - dm21sq;
  dmsq[5][2] = dm51sq - dm21sq;
  if(fMassHierarchy == kNmhNormal) 
    {
    dmsq[4][3] = dm41sq - dm31sq;
    dmsq[5][3] = dm51sq - dm31sq;
    }
  else if(fMassHierarchy == kNmhInverted)  
    {
    dmsq[4][3] = dm41sq + dm31sq;
    dmsq[5][3] = dm51sq + dm31sq;
    }

  usingOscParams = true;
}
//________________________________________________________________________________________
void NuOscSterile::SetDeltaCP(double dcp12, double dcp13, double dcp15, double dcp34, double dcp35)
{
  deltaCP[1][2] = dcp12;
  deltaCP[1][3] = dcp13;
  deltaCP[1][5] = dcp15;
  deltaCP[3][4] = dcp34;
  deltaCP[3][5] = dcp35;

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
  for(int j=0; j<=3; j++)
    {
    UVacuum[j][0][re] = 0.0;
    UVacuum[j][0][im] = 0.0;
    UVacuum[0][j][re] = 0.0;
    UVacuum[0][j][im] = 0.0;
    }

  UVacuum[1][1][re] = Ue1;
  UVacuum[1][2][re] = Ue2;
  UVacuum[1][3][re] = Ue3;
  UVacuum[1][4][re] = Ue4;
  UVacuum[1][5][re] = Ue5;
  UVacuum[2][1][re] = Umu1;
  UVacuum[2][2][re] = Umu2;
  UVacuum[2][3][re] = Umu3;
  UVacuum[2][4][re] = Umu4;
  UVacuum[2][5][re] = Umu5;
  UVacuum[3][1][re] = Utau1;
  UVacuum[3][2][re] = Utau2;
  UVacuum[3][3][re] = Utau3;
  UVacuum[3][4][re] = Utau4;
  UVacuum[3][5][re] = Utau5;
  UVacuum[4][1][re] = U41;
  UVacuum[4][2][re] = U42;
  UVacuum[4][3][re] = U43;
  UVacuum[4][4][re] = U44;
  UVacuum[4][5][re] = U45;
  UVacuum[5][1][re] = U51;
  UVacuum[5][2][re] = U52;
  UVacuum[5][3][re] = U53;
  UVacuum[5][4][re] = U54;
  UVacuum[5][5][re] = U55;

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
  for(int j=0; j<=3; j++)
    {
    UVacuum[j][0][re] = 0.0;
    UVacuum[j][0][im] = 0.0;
    UVacuum[0][j][re] = 0.0;
    UVacuum[0][j][im] = 0.0;
    }
  UVacuum[1][1][im] = Ue1;
  UVacuum[1][2][im] = Ue2;
  UVacuum[1][3][im] = Ue3;
  UVacuum[1][4][im] = Ue4;
  UVacuum[1][5][im] = Ue5;
  UVacuum[2][1][im] = Umu1;
  UVacuum[2][2][im] = Umu2;
  UVacuum[2][3][im] = Umu3;
  UVacuum[2][4][im] = Umu4;
  UVacuum[2][5][im] = Umu5;
  UVacuum[3][1][im] = Utau1;
  UVacuum[3][2][im] = Utau2;
  UVacuum[3][3][im] = Utau3;
  UVacuum[3][4][im] = Utau4;
  UVacuum[3][5][im] = Utau5;
  UVacuum[4][1][im] = U41;
  UVacuum[4][2][im] = U42;
  UVacuum[4][3][im] = U43;
  UVacuum[4][4][im] = U44;
  UVacuum[4][5][im] = U45;
  UVacuum[5][1][im] = U51;
  UVacuum[5][2][im] = U52;
  UVacuum[5][3][im] = U53;
  UVacuum[5][4][im] = U54;
  UVacuum[5][5][im] = U55;

  usingOscParams = false;

}
//________________________________________________________________________________________
void NuOscSterile::SetVacuumPMNSImToZero()
{
  for (int i=0; i<=nFlav; i++)
    {
    for (int j=0; j<=nFlav; j++)
      {
      UVacuum[i][j][im] = 0.0;
      }
    }

  usingOscParams = false;

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
  //Index 0 is not used and is set to 0 in init()
  //Indices 1 and 2 are those of matrix
  //Indices 3 and 4 are elements of that matrix
  //Index 5: 0 is real part, 1 is imaginary part of matrix element
  for (int i=1; i<=nFlav; i++)
    {
    for (int j=1; j<=nFlav; j++)
      {
      for (int k=1; k<=nFlav; k++)
	{
	for (int l=1; l<=nFlav; l++)
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
	    if(nFlav == 5 && ((i == 1 && (j == 2) || j == 3 || j == 5) || (i == 3 && (j == 4 || j == 5))))
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
	  else if(k == j && l == i)
	    {
	    if(nFlav == 5 && ((i == 1 && (j == 2) || j == 3 || j == 5) || (i == 3 && (j == 4 || j == 5))))
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
  for(int i=1; i<=nFlav; i++)
    {
      for(int j=1; j<=nFlav; j++)
	{
	  if(i == j)
	    UVacuum[i][j][re] = 1.0;
	  else
	    UVacuum[i][j][re] = 0.0;
	  UVacuum[i][j][im] = 0.0;
	}
    }

  double UVacuumCopy[nFlav+1][nFlav+1][2];

  //Multiply rotation matrices 12, 13, 14, 15, 23, 24, 25, 34 and 35 in that order.
  //The product of these matrices is the vacuum PMNS matrix UVacuum.
  //In general, rotation matrix 45 would be added to the product, but this is neglected
  //since oscillations from neutrino 4 to 5 or vice versa produce no observable effects. 
  for(int r=1; r<=3; r++)
    {
    for(int s=r+1; s<=nFlav; s++)
      {
	  
      for(int i=1; i<=nFlav; i++)
	{
	for(int j=1; j<=nFlav; j++)
	  {
	  UVacuumCopy[i][j][re] = UVacuum[i][j][re];
	  UVacuumCopy[i][j][im] = UVacuum[i][j][im];
          UVacuum[i][j][re] = 0.0;
	  UVacuum[i][j][im] = 0.0;
	  }
	}

      for(int i=1; i<=nFlav; i++)
	{
	for(int j=1; j<=nFlav; j++)
	  {
	  for(int k=1; k<=nFlav; k++)
	    {
            UVacuum[i][j][re] += rotationMatrix[r][s][i][k][re] * UVacuumCopy[k][j][re]
                               - rotationMatrix[r][s][i][k][im] * UVacuumCopy[k][j][im];
            UVacuum[i][j][im] += rotationMatrix[r][s][i][k][re] * UVacuumCopy[k][j][im]
                               + rotationMatrix[r][s][i][k][im] * UVacuumCopy[k][j][re];
	    }
	  }
	}

      }//end of for(int s=r+1; s<=nFlav; s++)
    }//end of for(int r=1; r<=3; r++)

}
//________________________________________________________________________________________
bool NuOscSterile::Check3FlavPMNS()
{
  //Check the first 3 rows and columns of the PMNS matrix against the standard 3-flavour parametrisation
  //To use this method, all mixing angles with indices 4 and 5 must be set to zero

  //Index 0 is not used and is set to 0.
  //Indices 1 - 3 are used for the 3 rows and columns of the 3-flavor PMNS matrix.
  double PMNSMatrix3Flav[4][4][2];
  for(int j=0; j<=3; j++)
    {
    PMNSMatrix3Flav[j][0][re] = 0.0;
    PMNSMatrix3Flav[j][0][im] = 0.0;
    PMNSMatrix3Flav[0][j][re] = 0.0;
    PMNSMatrix3Flav[0][j][im] = 0.0;
    }

  PMNSMatrix3Flav[1][1][re] = cosMixingAngle[1][2]*cosMixingAngle[1][3];
  PMNSMatrix3Flav[1][1][im] = 0.0;
  PMNSMatrix3Flav[1][2][re] = cosMixingAngle[1][3]*sinMixingAngle[1][2];
  PMNSMatrix3Flav[1][2][im] = 0.0;
  PMNSMatrix3Flav[1][3][re] = sinMixingAngle[1][3]*cos(deltaCP[1][3]);
  PMNSMatrix3Flav[1][3][im] = sinMixingAngle[1][3]*-sin(deltaCP[1][3]);
  PMNSMatrix3Flav[2][1][re] = -cosMixingAngle[2][3]*sinMixingAngle[1][2]-sinMixingAngle[1][3]*sinMixingAngle[2][3]*cosMixingAngle[1][2]*cos(deltaCP[1][3]);
  PMNSMatrix3Flav[2][1][im] = -sinMixingAngle[1][3]*sinMixingAngle[2][3]*cosMixingAngle[1][2]*sin(deltaCP[1][3]);
  PMNSMatrix3Flav[2][2][re] = cosMixingAngle[1][2]*cosMixingAngle[2][3]-sinMixingAngle[1][2]*sinMixingAngle[1][3]*sinMixingAngle[2][3]*cos(deltaCP[1][3]);
  PMNSMatrix3Flav[2][2][im] = -sinMixingAngle[1][2]*sinMixingAngle[1][3]*sinMixingAngle[2][3]*sin(deltaCP[1][3]);
  PMNSMatrix3Flav[2][3][re] = cosMixingAngle[1][3]*sinMixingAngle[2][3];
  PMNSMatrix3Flav[2][3][im] = 0.0;
  PMNSMatrix3Flav[3][1][re] = sinMixingAngle[1][2]*sinMixingAngle[2][3]-sinMixingAngle[1][3]*cosMixingAngle[1][2]*cosMixingAngle[2][3]*cos(deltaCP[1][3]);
  PMNSMatrix3Flav[3][1][im] = -sinMixingAngle[1][3]*cosMixingAngle[1][2]*cosMixingAngle[2][3]*sin(deltaCP[1][3]);
  PMNSMatrix3Flav[3][2][re] = -sinMixingAngle[2][3]*cosMixingAngle[1][2]-sinMixingAngle[1][2]*cosMixingAngle[2][3]*sinMixingAngle[1][3]*cos(deltaCP[1][3]);
  PMNSMatrix3Flav[3][2][im] = -sinMixingAngle[1][2]*cosMixingAngle[2][3]*sinMixingAngle[1][3]*sin(deltaCP[1][3]);
  PMNSMatrix3Flav[3][3][re] = cosMixingAngle[1][3]*cosMixingAngle[2][3];
  PMNSMatrix3Flav[3][3][im] = 0.0;

  this->CalcVacuumPMNS();

  for(int i=1; i<=3; i++)
    {
    for(int j=1; j<=3; j++)
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
  //The 3-flavour formula is found in many papers and theses, e.g. equation 11 in hep-ph 0506165.
  //The neutrino codes are the PDG ones and the neutrino energy is in GeV.

  OscChannelId_t osch = util::GetOscChannel(nu_start, nu_end);

  if(osch == kOscUndefined) {
    cerr << "NuOsc3Flavour: Not allowing (neutrino = " << nu_start << ") --> "
         << "(neutrino = " << nu_end << ") oscillations" << endl;
    return 0.;
  }

  double oscProb;

  if(usingOscParams)
    this->CalcVacuumPMNS();

  // nue->nue, nuebar->nuebar
  if(osch == kOscNuE2NuE || osch == kOscNuEBar2NuEBar)
    {
    oscProb = 1.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
        {
	this->CalcMEProduct(1, 1, i, j);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    } 

  // numu->nue, nuebar->numubar 
  else if(osch == kOscNuMu2NuE || osch == kOscNuEBar2NuMuBar)
    { 
    oscProb = 0.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
	{
	this->CalcMEProduct(2, 1, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);	
        }
      } 
    }

  // numubar->nuebar, nue->numu
  else if(osch == kOscNuMuBar2NuEBar || osch == kOscNuE2NuMu)
    {
    oscProb = 0.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
	{
	this->CalcMEProduct(2, 1, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  // numu->numu, numubar->numubar
  else if(osch == kOscNuMu2NuMu || osch == kOscNuMuBar2NuMuBar)
    {
    oscProb = 1.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
	{
	this->CalcMEProduct(2, 2, i, j);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
	        * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  // nue->nutau, nutaubar->nuebar
  else if(osch == kOscNuE2NuTau || osch == kOscNuTauBar2NuEBar)
    {
    oscProb = 0.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
	{
	this->CalcMEProduct(3, 1, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  // nuebar->nutaubar, nutau->nue
  else if(osch == kOscNuEBar2NuTauBar || osch == kOscNuTau2NuE)
    {
    oscProb = 0.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
	{
	this->CalcMEProduct(3, 1, i, j);
	oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  // numu->nutau, nutaubar->numubar
  else if(osch == kOscNuMu2NuTau || osch == kOscNuTauBar2NuMuBar)
     {
     oscProb = 0.0;
     for(int i=1; i<=3; i++)
       {
       for(int j=i+1; j<=nFlav; j++)
	 {
	 this->CalcMEProduct(2, 3, i, j);
	 oscProb += 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	 oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		 * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	 }
       }
     }

  // numubar->nutaubar, nutau->numu
  else if(osch == kOscNuMuBar2NuTauBar || osch == kOscNuTau2NuMu)
    {
    oscProb = 0.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
	{
	this->CalcMEProduct(2, 3, i, j);
	oscProb -= 2 * MEProduct[im] * sin((2.0 * kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
		* sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
       }
     }
   }

  // nutau->nutau, nutaubar->nutaubar
  else if(osch == kOscNuTau2NuTau || osch == kOscNuTauBar2NuTauBar)
    {
    oscProb = 1.0;
    for(int i=1; i<=3; i++)
      {
      for(int j=i+1; j<=nFlav; j++)
        {
	this->CalcMEProduct(3, 3, i, j);
	oscProb -= 4.0 * MEProduct[re] * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy)
	         * sin((kOscDimensFactor * dmsq[j][i] * fBaseline) / energy);
	}
      }
    }

  return oscProb;
  
}
