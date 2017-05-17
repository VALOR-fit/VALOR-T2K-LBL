//________________________________________________________________________________________
/*

   Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

21-7-2014 Davide
Added samples for multiring forward horn current

23-06-2016 Chris
Added check for parameter maxima and minima since these are decoupled from the 
FitParamList

*/
//________________________________________________________________________________________

#include <cassert>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <TSystem.h>

#include "valor/Messenger/Messenger.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/Registry/Registry.h"
#include "valor/Utilities/UtlSystem.h"

#include "UtlFitParamListFactory.h"
#include "Enumerations.h"
#include "Conventions.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"
#include "Utilities.h"
#include "T2KPhysicalParameterLimits.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace valor;
using namespace valor::conv;

FitParamList * util::NamedParameterList(string name)
{

  //
  // 3-flavour oscillation fit with/without systematics included 
  //

  FitParamList * fp = 0;

  if(name == "nuosc3f") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }
  if(name == "nuosc3f+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);

    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }
  if(name == "nuosc3f_12constrain+escale") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    int np = fp->GetNumberParameters();
    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = fp->WhichParam(ip);
      if(!fp->IsNuisance(p)) continue; 
      if(p != kNuisSuperkEnergyScale &&
          p != kPhypMixing12 &&
          p != kPhypDm12sq ) 
        fp -> FixParameter(p);
    }
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
  }

  //
  // 3-flavour oscillation fit with various oscillation param subsets either fixed or
  // treated as nuisance paramers.
  //

  if(name == "nuosc3f_12constrain") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }
  if(name == "nuosc3f_12constrain+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }
  if(name == "nuosc3f_12fix") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_12fix+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }


  // NEW FHC+RHC JOINT
  if(name == "nuosc3f_fixall_fitdcp") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);    
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }
  if(name == "nuosc3f_fixall_fitdcp+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);    
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }
  //

  if(name == "nuosc3f_numudisapp_fix4") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_numudisapp_nuis3fixdcp") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance (kPhypMixing13, true);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq, true);
    fp -> FixParameter (kPhypCPV);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3f_numudisapp_fix4+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);

    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //#ifdef __FHC_RUNNING__
    //~ fp -> FixParameter (kNuisXSecRPA ); 
    //#endif
    //#ifdef __RHC_RUNNING__    
    //fp -> FixParameter (kRHCNuisXSecRPA );  // xsec syst are fully RHC-FHC correlated
    //#endif
#endif

  }

  if(name == "nuosc3f_numudisapp_nuis3fixdcp+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing13, true);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq, true);
    fp -> FixParameter (kPhypCPV);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3f_nueapp_fix4") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }
  if(name == "nuosc3f_nueapp_fix4+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3f_dcpsensitivity_fix5") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_dcpsensitivity_fix3") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(4);//use future osc errors
  }

  if(name == "nuosc3f_dcpsensitivity") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(3);//use current osc errors
  }

  if(name == "nuosc3f_dcpsensitivity_oscnuis") {
    fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(5);//use future osc errors
  }

  if(name == "nuosc3f_dcpsensitivity_fix3+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(5);//use future osc errors
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3f_dcpsensitivity+nuisance") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(6);//use current osc errors
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3f_dcpsensitivity_allsyst") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    //        fp -> SetNuisance  (kPhypDm23sq,   true);
    //        fp -> SetNuisance  (kPhypMixing13, true);
    fp -> FixParameter (kPhypMixing12      );
    fp -> FixParameter (kPhypDm12sq        );
    //fp -> FixParameter (kPhypCPV);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);//future oscillation
  }

  if(name == "nuosc3fbar_numudissap_allsyst") {
    fp = FitParamList::NuOsc3FlavourBarDissapFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 8);
    fp -> FixParameter (kPhypCPV);
    fp -> SetNuisance  (kPhypMixing13 , true);
    fp -> SetNuisance  (kPhypMixing12 , true);
    fp -> SetNuisance  (kPhypDm12sq   , true);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_numudisapp_constrainall") {
    fp = FitParamList::NuOsc3FlavourBarFit();
    fp -> SetNuisance  (kPhypMixing13 , true);
    fp -> SetNuisance  (kPhypMixing23 , true);
    fp -> SetNuisance  (kPhypMixing12 , true);
    fp -> SetNuisance  (kPhypCPV      , true);
    fp -> SetNuisance  (kPhypDm12sq   , true);
    fp -> SetNuisance  (kPhypDm23sq   , true);

    fp -> SetNuisance  (kPhypMixing13Bar , true);
    fp -> SetNuisance  (kPhypMixing12Bar , true);
    fp -> SetNuisance  (kPhypCPVBar      , true);
    fp -> SetNuisance  (kPhypDm12sqBar   , true);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_numudisapp_constrainall+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);
    fp -> SetNuisance  (kPhypMixing13 , true);
    fp -> SetNuisance  (kPhypMixing23 , true);
    fp -> SetNuisance  (kPhypMixing12 , true);
    fp -> SetNuisance  (kPhypCPV      , true);
    fp -> SetNuisance  (kPhypDm12sq   , true);
    fp -> SetNuisance  (kPhypDm23sq   , true);

    fp -> SetNuisance  (kPhypMixing13Bar , true);
    fp -> SetNuisance  (kPhypMixing12Bar , true);
    fp -> SetNuisance  (kPhypCPVBar      , true);
    fp -> SetNuisance  (kPhypDm12sqBar   , true);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);

    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //

  }

  if(name == "nuosc3fbar_numudisapp+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);

#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA );
#endif
  }


  if(name == "nuosc3fbar_numudisapp_fixnonbar_fix4") {
    fp = FitParamList::NuOsc3FlavourBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);

    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //

  }



  if(name == "nuosc3fbar_numudisapp_fixnonbar_fix5+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);

    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //

  }




  if(name == "nuosc3fbar_fixnonbar_12fix") {
    fp = FitParamList::NuOsc3FlavourBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar  );

    //fp -> FixParameter (kPhypCPVBar  );
    //fp -> SetNuisance  (kPhypMixing23Bar , true);
    //fp -> SetNuisance  (kPhypDm23sqBar   , true);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_fixnonbar_12fix+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar  );

    //fp -> FixParameter (kPhypMixing13Bar);
    //fp -> FixParameter (kPhypCPVBar  );
    //fp -> SetNuisance  (kPhypMixing23Bar , true);
    //fp -> SetNuisance  (kPhypDm23sqBar   , true);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3fbar_fixnonbar_fix5") {
    fp = FitParamList::NuOsc3FlavourBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar  );
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypDm23sqBar);
    fp -> FixParameter (kPhypCPVBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_fixnonbar_fix5+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar  );
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypDm23sqBar);
    fp -> FixParameter (kPhypCPVBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }


  if(name == "nuosc3fbar_nueapp_fixnonbar_fix4") {
    fp = FitParamList::NuOsc3FlavourBarFit();
    //util::AppendAllSystematicsTo(fp, 12);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar  );
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_nueapp_fixnonbar_fix4+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 12);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar  );
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }





  // NEW BETABAR

  if(name == "nuosc3fbar_normbar_fixall") {
    fp = FitParamList::NuOsc3FlavourBarNormBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> FixParameter (kPhypBetaBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }


  if(name == "nuosc3fbar_normbar_fixall+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 13);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> FixParameter (kPhypBetaBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3fbar_normbar_fixosc") {
    fp = FitParamList::NuOsc3FlavourBarNormBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_normbar_fixosc+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 13);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypCPVBar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3fbar_normbar_fitdcpbar") {
    fp = FitParamList::NuOsc3FlavourBarNormBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3fbar_normbar_fitdcpbar+nuisance") {
    fp = FitParamList::NuOsc3FlavourBarNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 13);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> FixParameter (kPhypMixing13Bar);
    fp -> FixParameter (kPhypMixing23Bar);
    fp -> FixParameter (kPhypMixing12Bar);
    fp -> FixParameter (kPhypDm12sqBar);
    fp -> FixParameter (kPhypDm23sqBar);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  //

  // NEW BETABAR 2

  if(name == "nuosc3f_normbar_fitdcp") {
    fp = FitParamList::NuOsc3FlavourNormBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_normbar_fitdcp+nuisance") {
    fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 7);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }


  if(name == "nuosc3f_normbar_fit23_fitdcp") {
    fp = FitParamList::NuOsc3FlavourNormBarFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_normbar_fit23_fitdcp+nuisance") {
    fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 7);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }


  if(name == "nuosc3f_normbar_fix12") {
    fp = FitParamList::NuOsc3FlavourNormBarFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_normbar_fix12+nuisance") {
    fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 7);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif
    //
  }

  if(name == "nuosc3f_normbar_fixall+nuisance") {
    fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 7);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);

    fp -> SetParListName(name);
    fp -> SetSystMode(0);
  }

  if(name == "nuosc3f_normbar+nuisance") {
    fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);

    util::AppendAllSystematicsTo(fp, 7);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
#if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
#endif

    //
  }
  if(name == "nuosc3f+nuisance+fix_12_23") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);
  }
  if(name == "nuosc3f+nuisance+fix_noncp") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);
  }
  if(name == "nuosc3f+nuisance+fix_non23") {
    fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypCPV);
  }

  if (!fp) {
    cerr<<"FitParamListFactory/NamedParamList >> Unable to make: "<<name<<endl;
    util::Die(-1); return fp;
  }
  
  // Do a quick check to ensure our parameter list is full
  if (! fp->IsComplete() ) {
    LOG("UtlFitParamListFactory",pERROR)<<"Fit parameter list is not complete";
    fp->Print(/*LOG("UtlFitParamListFactory",pERROR)*/);
    util::Die(-1);
  }
  // Make sure that the parameter limits that are set aren't in unphysical regions
  T2KPhysicalParameterLimits::Instance()->SetLimits(fp);
  return fp;
}

//________________________________________________________________________________________
void util::AppendAllSystematicsTo(FitParamList * fparm, int start)
{
// Add nuisance params
//


#if defined(__NUEONLY__) // nue appearance
LOG("UtlFitParamListFactory",pINFO) << "Determine Append";
  #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
    
    #if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__)
    // fhc all + skdetfsi multiring
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    #if defined(__MULTIRING_FHC_NUE__)
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRe);
    #elif defined(__MULTIRING_FHC_NUMU__)
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRmu);
    #else
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing);
    #endif
    
    #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
    
      #if defined(__FULLY_CORRELATED__)
      // fhc+rhc fully correlated 
      util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
      util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
      util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
    
      // NEW FHC+RHC JOINT
      #else
      // fhc+rhc w/ correlations 
      util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
      util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfFHCpRHCParamsInBANFFfit);
      util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfFHCpRHCParamsInBANFFfit+1+kNumOfFHCpRHCSKDetFSIParams);    
      #endif
      //

    #elif defined(__FHC_RUNNING__)
    // fhc only 
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);

    #elif defined(__RHC_RUNNING__)
    // rhc only 
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfRHCParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfRHCParamsInBANFFfit+1+kNumOfRHCSKDetFSIParams);
    //

    #endif
    
  #else

    #ifdef __MULTIRING_FHC_RUNNING__
    // multiring fhc only 
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    #if defined(__MULTIRING_FHC_NUE__)
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRe);
    #elif defined(__MULTIRING_FHC_NUMU__)
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRmu);
    #else
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing);
    #endif
    
    #endif

  #endif // sample analysis option

#else // joint fit and NumuDissap

#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)

    #if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__)

    // fhc all + skdetfsi multiring 
    util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
    util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit);
    #if defined(__MULTIRING_FHC_NUE__)
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRe);
    #elif defined(__MULTIRING_FHC_NUMU__)
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRmu);
    #else
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing);
    #endif
    
    #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

      #if defined(__FULLY_CORRELATED__)
      // fhc+rhc fully correlated 
      util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
      util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit); 
      util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);

      // NEW FHC+RHC JOINT
      #else
      // fhc+rhc w/ correlations 
      util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
      util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfFHCpRHCParamsInBANFFfit); 
      util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfFHCpRHCParamsInBANFFfit+1+kNumOfFHCpRHCSKDetFSIParams);    
      #endif
      //

    #elif defined(__FHC_RUNNING__)
    // fhc only 
    util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
    util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit); 
    util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
  
    #elif defined(__RHC_RUNNING__)
    // rhc only 
    util::AppendAllCorrelatedBANFFSystematics  (fparm, start);
    util::AppendAllSuperkDetectorFSISystematics(fparm, start+kNumOfRHCParamsInBANFFfit);
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfRHCParamsInBANFFfit+1+kNumOfRHCSKDetFSIParams);
    //
    #endif

  #else
    #ifdef __MULTIRING_FHC_RUNNING__
    // multiring fhc only 
    util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
    util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit);
    #if defined(__MULTIRING_FHC_NUE__)
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRe);
    #elif defined(__MULTIRING_FHC_NUMU__)
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRmu);
    #else
    util::AppendAllUncorrelatedXSecSystematics (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing);
    #endif

    #endif

  #endif // sample analysis option  
#endif // oscillation analysis option

}
//________________________________________________________________________________________

//________________________________________________________________________________________
void util::AppendAllUncorrelatedXSecSystematics(
#if defined(__USE_NIWG2013__)
  FitParamList * fparm, int start)
#else
  FitParamList * /*fparm*/, int /*start*/)
#endif
{
// Add nuisance params
//


  // NEW XSEC --> comment all old non-banff xsec systematics --> Not needed any longer

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)
  int ip = start; // first unallocated element in the param list
  
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom));
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc));
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape));
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth));
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta));
  #ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE));
  #endif
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh));
  #ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNC1piC));
  #endif
  #ifndef __NUMUONLY__
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP*0.3,  false, true, AsString(kNuisXSecNormNCcoh));
  #endif
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth));
  #ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu));
  #endif
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu));

#elif defined(__USE_NIWG2015__)
  // NEW XSEC
  //~ fparm->AddParam( ip++, kNuisXSecRPA,           0., 0.,                1.,              SYST_STEP,  false, true, AsString(kNuisXSecRPA) );
  //~ fparm->AddParam( ip++, kRHCNuisXSecRPA,           0., 0.,                1.,              SYST_STEP,  false, true, AsString(kRHCNuisXSecRPA) );
#endif

}

#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsApp(
#if defined(__USE_NIWG2013__)
  FitParamList * fparm, int start)
#else
  FitParamList * /*fparm*/, int /*start*/)
#endif
{
// Add nuisance params
//

  // NEW XSEC --> comment all old non-banff xsec systematics

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)
  int ip = start; // first unallocated element in the param list

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom) );
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc) );
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape) );
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth) );
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta) );
  //  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE) );
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth) );
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu) );
  //  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu) );

#elif defined(__USE_NIWG2015__)
  // NEW XSEC
  //fparm->AddParam( ip++, kNuisXSecRPA,           0., 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecRPA) );
#endif

}
#endif
//________________________________________________________________________________________
void util::AppendAllSuperkDetectorFSISystematics(FitParamList * fparm, int /*start*/)
{
  // Add nuisance params
  //
  const Registry & reg = Opt::Instance()->GlobalRegistry();
  
  const SuperkDetectorFSI * skdetfsi = dynamic_cast<const SuperkDetectorFSI*>(reg.GetObject("skdetfsi"));
  
  if (!skdetfsi) {
    LOG("UtilFitParamListFactory",pFATAL)<<"Failed to get skdetfsi from Opt::GlobalRegistry";
    util::Die(-1); return;
  }
  
  skdetfsi->AppendParamsTo(*fparm);
}
//________________________________________________________________________________________
void util::AppendAllCorrelatedBANFFSystematics(FitParamList * fparm, int start)
{
  // Add nuisance params
  //
  const Registry & reg = Opt::Instance()->GlobalRegistry();
  
  const BANFF * banff = dynamic_cast<const BANFF*>(reg.GetObject("banff"));
  
  if (!banff) {
    LOG("UtilFitParamListFactory",pFATAL)<<"Failed to get banff from Opt::GlobalRegistry";
    util::Die(-1); return;
  }
  
  banff->AppendParamsTo(*fparm);
}
//________________________________________________________________________________________
