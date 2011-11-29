//____________________________________________________________________________
/*!

\class    t2kana::NuOsc3Flavour

\brief    Neutrino-oscillation probabilities (incl. matter effects) in 
          computed in a 3-flavour framework.

\author   Nick Grant <n.grant@lancaster.ac.uk>

\created  July 05, 2011

*/
//____________________________________________________________________________

#ifndef _NuOsc3Flavour_
#define _NuOsc3Flavour_

#include <iostream>
#include <math.h>

#include "Conventions.h"

namespace t2kana {

class NuOsc3Flavour
{
public:     
  NuOsc3Flavour();
 ~NuOsc3Flavour();
 
  //
  // Set baseline, earth density, mass hierarchy, squared mass-splittings and
  // mixing angles / phases.
  //

  void SetBaseline      (double length);
  void SetDensity       (double density);
  void SetTheta23Octant (int octant);
  void SetMassHierarchy (t2kana::conv::NuMassHierarchy_t h); 
  void SetOscParameters (double sinsq2theta12, /* sin^2(2*theta_{12}) */
                         double sinsq2theta13, /* sin^2(2*theta_{13}) */
                         double sinsq2theta23, /* sin^2(2*theta_{23}) */
                         double dmsq12,        /*  Dm^2_{12} in eV^2/c^4 */
                         double dmsq23,        /* |Dm^2_{23}| (always positive) in eV^2/c^4 */
                         double deltaCP        /* delta_{CP} phase in radians */
                        );
  void SetOscParametersSinSqT23 (double sinsq2theta12, /* sin^2(2*theta_{12}) */
                                 double sinsq2theta13, /* sin^2(2*theta_{13}) */
                                 double sinsqtheta23,  /* sin^2(theta_{23}) */
                                 double dmsq12,        /*  Dm^2_{12} in eV^2/c^4 */
                                 double dmsq23,        /* |Dm^2_{23}| (always positive) in eV^2/c^4 */
                                 double deltaCP        /* delta_{CP} phase in radians */
                                 );

 //         
 // Calculate 3-flavour oscillation probabilities in vacuum and in matter.
 // The neutrino codes are the PDG ones and the neutrino energy is in GeV.
 //         

 double CalcVacuumOscProb      (int nu_start, int nu_end, double energy);                
 double CalcVacuumOscProbPDG   (int nu_start, int nu_end, double energy);
 double CalcMatterOscProb      (int nu_start, int nu_end, double energy);
 double Calc2FlavMatterOscProb (int nu_start, int nu_end, double energy, int index);
 
private:

  enum {
    re,
    im
  };

 void Init(void);

 t2kana::conv::NuMassHierarchy_t fMassHierarchy; ///< mass hierarchy
 double fBaseline; ///< baseline (km)
 double fNe;       ///< number density of electrons in Earth's crust (cm^-3)

 // Squared mass splittings
 double dmsq23;    //< eV^2
 double dmsq12;    //< eV^2
 double dmsq13;    //< eV^2 

 // Moduli squared of PMNS matrix elements
 double modUe1sq;
 double modUe2sq;
 double modUe3sq;
 double modUmu1sq;
 double modUmu2sq;
 double modUmu3sq;
 double modUtau1sq;
 double modUtau2sq;
 double modUtau3sq;

 double Re12;
 double Im12;
 double arg12;
 double Re13;
 double Im13;
 double arg13;
 double Re23;
 double Im23;
 double arg23;

 // The mass squared matrix has diagonal elements 0, dmsq12, dmsq13
 // and all off-diagonal elements are zero
 // UMassUdaggerVacuum is the product of the PMNS matrix, the mass 
 // squared matrix and the dagger of the PMNS matrix in vacuum.
 double UMassUdaggerVacuum[3][3][2];		

 bool useEquations12;
 int nFlavourIndex;
 int theta23Octant;
};

}//end of namespace t2kana

#endif

