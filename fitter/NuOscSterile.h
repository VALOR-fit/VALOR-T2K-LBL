
//____________________________________________________________________________
/*!

\class    t2kana::NuOscSterile

\brief    Neutrino-oscillation probabilities (incl. matter effects) for 
          1 or 2 sterile neutrinos in addition to the 3 active neutrinos.

\author   Nick Grant <n.grant@lancaster.ac.uk>

\created  14th January 2014

*/
//____________________________________________________________________________

#ifndef _NuOscSterile_
#define _NuOscSterile_

#include <iostream>
#include <math.h>

#include "Conventions.h"

namespace t2kana {

class NuOscSterile
{
 public:     
  NuOscSterile();
  ~NuOscSterile();

  static int const nFlav = 5;

  //
  // Set baseline, earth density, mass hierarchy, squared mass-splittings and
  // mixing angles / phases.
  //

  void SetBaseline      (double length);
  void SetDensity       (double density);
  void SetMassHierarchy (t2kana::conv::NuMassHierarchy_t h);

  void SetMixingAngles  (double sinsqtheta12, 
                         double sinsqtheta13, 
                         double sinsqtheta14, 
			 double sinsqtheta15,
                         double sinsqtheta23,
                         double sinsqtheta24,
			 double sinsqtheta25,
                         double sinsqtheta34,
                         double sinsqtheta35
			 );

  void SetDmsq          (double dm21sq,
			 double dm31sq,
			 double dm41sq,
			 double dm51sq
			 );

  void SetDeltaCP       (double dcp12,
		         double dcp13,
		         double dcp15,
		         double dcp34,
		         double dcp35
			 );

  void SetVacuumPMNSReal(double Ue1, double Ue2, double Ue3, double Ue4, double Ue5,
                         double Umu1, double Umu2, double Umu3, double Umu4, double Umu5,
                         double Utau1, double Utau2, double Utau3, double Utau4, double Utau5,
                         double U41, double U42, double U43, double U44, double U45,
                         double U51, double U52, double U53, double U54, double U55
                         );

  void SetVacuumPMNSIm  (double Ue1, double Ue2, double Ue3, double Ue4, double Ue5,
                         double Umu1, double Umu2, double Umu3, double Umu4, double Umu5,
                         double Utau1, double Utau2, double Utau3, double Utau4, double Utau5,
                         double U41, double U42, double U43, double U44, double U45,
                         double U51, double U52, double U53, double U54, double U55
			 );

  void SetVacuumPMNSImToZero();
  void CalcVacuumPMNS   ();
  bool Check3FlavPMNS   ();
  void CalcMEProduct    (int nu_start, 
                        int nu_end,
                        int mass1,
                        int mass2
		        ); 

  double CalcVacuumOscProb(int nu_start,
			   int nu_end,
                           double energy);

  double CalcMatterOscProb(int nu_start,
                           int nu_end,
                           double energy);
			   
  // Calculate all matter probabilities at once.
  // Osc_probs_out should point to a block of memory nflav*nflav
  // It fills this with an array[i][j] of P(nu_i)->P(nu_j)
  void CalcMatterOscProbAll(bool isAntineutrino, double energy, double * osc_probs_out);

 private:

  enum {
    re,
    im
  };

  void Init(void);

  t2kana::conv::NuMassHierarchy_t fMassHierarchy; ///< mass hierarchy
  double fBaseline; ///< baseline (km)
  double fNe;       ///< number density of electrons in Earth's crust (cm^-3)

  bool usingOscParams; 

  double sinMixingAngle[nFlav][nFlav];
  double cosMixingAngle[nFlav][nFlav];
  double dmsq[nFlav][nFlav];
  double deltaCP[nFlav][nFlav];

  double rotationMatrix[nFlav][nFlav][nFlav][nFlav][2];
  double UVacuum[nFlav][nFlav][2];
  double MEProduct[2];

};

}//end of namespace t2kana

#endif
