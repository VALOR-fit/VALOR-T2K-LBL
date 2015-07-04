//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

21-7-2014 Davide
Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

#include <cassert>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <TSystem.h>

#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"

#include "UtlFitParamListFactory.h"
#include "FitCore/FitParamList.h"

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

  if(name == "nuosc3f") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    // NEW XSEC
    #if defined(__USE_NIWG2015__)
    //~ fp -> FixParameter (kNuisXSecRPA ); 
    #endif
    //
    return fp;
  }
  if(name == "nuosc3f_12constrain+escale") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
  }
  
  //
  // 3-flavour oscillation fit with various oscillation param subsets either fixed or
  // treated as nuisance paramers.
  //

  if(name == "nuosc3f_12constrain") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_12constrain+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
  }
  if(name == "nuosc3f_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  if(name == "nuosc3f_12fix+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
   }

  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_nuis3fixdcp") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
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
    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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

    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_nuis3fixdcp+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
  }

  if(name == "nuosc3f_nueapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix5") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(4);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(3);//use current osc errors
    return fp;
  }
  
  if(name == "nuosc3f_dcpsensitivity_oscnuis") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(5);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
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
    return fp;
  }
  
  if(name == "nuosc3f_dcpsensitivity_allsyst") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    //        fp -> SetNuisance  (kPhypDm23sq,   true);
    //        fp -> SetNuisance  (kPhypMixing13, true);
    fp -> FixParameter (kPhypMixing12      );
    fp -> FixParameter (kPhypDm12sq        );
    //fp -> FixParameter (kPhypCPV);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);//future oscillation
    return fp;
  }

  if(name == "nuosc3fbar_numudisapp_constrainall") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit();
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
    return fp;
  }

  if(name == "nuosc3fbar_numudisapp_constrainall+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
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

    return fp;
  }
  
  if(name == "nuosc3fbar_numudisapp_fixnonbar_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit();
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
    
    return fp;
  }

  if(name == "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
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

    return fp;
  }



  if(name == "nuosc3fbar_numudisapp_fixnonbar_fix5+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
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

    return fp;
  }




  if(name == "nuosc3fbar_fixnonbar_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit();
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
    return fp;
  }
 
  if(name == "nuosc3fbar_fixnonbar_12fix+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
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
    return fp;
  }
  
  




  if(name == "nuosc3fbar_fixnonbar_fix5") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit();
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
    return fp;
  }
 
  if(name == "nuosc3fbar_fixnonbar_fix5+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
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
    return fp;
  }
  

  if(name == "nuosc3fbar_nueapp_fixnonbar_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit();
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
    return fp;
  }

  if(name == "nuosc3fbar_nueapp_fixnonbar_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarFit(kNumOfNuisanceParams);
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
    return fp;
  }
  
  
  
  
  
  // NEW BETABAR
  
  if(name == "nuosc3fbar_normbar_fixall") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarNormBarFit();
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
    return fp;
  }
  
  if(name == "nuosc3fbar_normbar_fixall+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarNormBarFit(kNumOfNuisanceParams);
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
    return fp;
  }
  
  if(name == "nuosc3fbar_normbar_fixosc") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarNormBarFit();
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
    return fp;
  }
    
    if(name == "nuosc3fbar_normbar_fixosc+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarNormBarFit(kNumOfNuisanceParams);
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
    return fp;
  }

  if(name == "nuosc3fbar_normbar_fitdcpbar") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarNormBarFit();
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
    return fp;
  }
    
    if(name == "nuosc3fbar_normbar_fitdcpbar+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourBarNormBarFit(kNumOfNuisanceParams);
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
    return fp;
  }

    //
    
    // NEW BETABAR 2
    
    if(name == "nuosc3f_normbar_fitdcp") {
      FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit();
      fp -> FixParameter (kPhypMixing13);
      fp -> FixParameter (kPhypMixing23);
      fp -> FixParameter (kPhypMixing12);
      fp -> FixParameter (kPhypDm12sq);
      fp -> FixParameter (kPhypDm23sq);
            
      fp -> SetParListName(name);
      fp -> SetSystMode(0);
      return fp;
    }
    
    if(name == "nuosc3f_normbar_fitdcp+nuisance") {
      FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
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
      return fp;
    }


    if(name == "nuosc3f_normbar_fit23_fitdcp") {
      FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit();
      fp -> FixParameter (kPhypMixing13);
      fp -> FixParameter (kPhypMixing12);
      fp -> FixParameter (kPhypDm12sq);
            
      fp -> SetParListName(name);
      fp -> SetSystMode(0);
      return fp;
    }
    
    if(name == "nuosc3f_normbar_fit23_fitdcp+nuisance") {
      FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
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
      return fp;
    }
    

    if(name == "nuosc3f_normbar_fix12") {
      FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit();
      fp -> FixParameter (kPhypMixing12);
      fp -> FixParameter (kPhypDm12sq);
            
      fp -> SetParListName(name);
      fp -> SetSystMode(0);
      return fp;
    }
    
    if(name == "nuosc3f_normbar_fix12+nuisance") {
      FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
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
      return fp;
    }



 if(name == "nuosc3f_normbar_fixall+nuisance") {
   FitParamList * fp = FitParamList::NuOsc3FlavourNormBarFit(kNumOfNuisanceParams);
   util::AppendAllSystematicsTo(fp, 7);
   fp -> FixParameter (kPhypMixing13);
   fp -> FixParameter (kPhypMixing23);
   fp -> FixParameter (kPhypMixing12);
   fp -> FixParameter (kPhypCPV);
   fp -> FixParameter (kPhypDm12sq);
   fp -> FixParameter (kPhypDm23sq);
      
   fp -> SetParListName(name);
   fp -> SetSystMode(0);
   return fp;
 }
 
 //  
  
 cerr<<"FitParamListFactory/NamedParamList >> Unable to make: "<<name<<endl;
 exit(-1); 
}
    
//________________________________________________________________________________________
void util::AppendAllSystematicsTo(FitParamList * fparm, int start)
{
// Add nuisance params
//

#if defined(__NUMUONLY__) // numu disappearance only

  #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  

    #if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__)
    // fhc all + skdetfsi multiring
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing+1);
    
    #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING) 

      #if defined(__FULLY_CORRELATED__) 
      // rhc+fhc fully correlated 
      util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
      util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
      util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
    
      // NEW FHC+RHC JOINT
      #else
      // rhc+fhc w/ correlations 
      util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
      util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfRHCpFHCParamsInBANFFfit);
      util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfRHCpFHCParamsInBANFFfit+1+kNumOfFHCpRHCSKDetFSIParams);
      #endif
      //
    
    #elif defined(__FHC_RUNNING__)
    // fhc only 
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);

    #elif defined(__RHC_RUNNING__)
    // rhc only 
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfRHCParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfRHCParamsInBANFFfit+1+kNumOfRHCSKDetFSIParams);
    
    #endif

  #else
    
   #ifdef __MULTIRING_FHC_RUNNING__
    // multiring fhc only 
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParamsMultiRing);
    // #endif
    #endif // sample analysis option

  #endif

#elif defined(__NUEONLY__) // nue appearance

  #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
    
    #if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__)
    // fhc all + skdetfsi multiring
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing+1);
    
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
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParamsMultiRing);
    #endif

  #endif // sample analysis option

#else // joint fit

  #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
    #if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__)

    // fhc all + skdetfsi multiring 
    util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
    util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit);          
    util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing+1);  

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
    util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParamsMultiRing); 
    #endif

  #endif // sample analysis option  
#endif // oscillation analysis option

}
//________________________________________________________________________________________

//________________________________________________________________________________________
void util::AppendAllUncorrelatedXSecSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // NEW XSEC --> comment all old non-banff xsec systematics --> Not needed any longer

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)
  
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
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  // NEW XSEC --> comment all old non-banff xsec systematics
  
  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom) );
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc) );
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape) );
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth) );
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta) );
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE) );
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNC1piC) );
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth) );
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu) );
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu) );

#elif defined(__USE_NIWG2015__)
  // NEW XSEC
  //fparm->AddParam( ip++, kNuisXSecRPA,           0., 0.,               1.,             SYST_STEP,  false, true, AsString(kNuisXSecRPA) );
#endif

}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  // NEW XSEC --> comment all old non-banff xsec systematics

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)

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

void util::AppendAllSuperkDetectorFSISystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // Get the covariance matrices & best fits of FSI & BANFF from xml, 
  // to have something sensible to initialise to
  // Parse the input XML configuration file
  
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();  

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "util/AppendAllSuperkDetectorFSISystematics >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathSKDetFSICovMatr.c_str());
  
  if(!SKDetFSICovMtrx) {
    cerr << "util/AppendAllSuperkDetectorFSISystematics >> " << "No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //~ xml->FreeNode(root_node);
  xml->FreeDoc(xmldoc);
  delete xml;
 
 
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  

  // fhc all + skdetfsi multiring
  #if defined (__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__) 
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );   // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale) );   // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRing; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }

  // NEW FLUX
  // fhc+rhc fully correlated   
  #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

    #if defined(__FULLY_CORRELATED__) 
    fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );   // FHC SK escale
    for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    
    // NEW FHC+RHC JOINT
    
    // fhc+rhc w/ correlations
    #else

    fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale)    );   // FHC SK escale
    for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    
    fparm->AddParam( ip++, kRHCNuisSuperkEnergyScale,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kRHCNuisSuperkEnergyScale) );   // RHC SK escale
    for(int i = 0; i < kNumOfRHCSKDetFSIParams; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    // 
    #endif

    //

  // fhc only
  #elif defined(__FHC_RUNNING__)    
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );   // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }

  // rhc only  
  #elif defined(__RHC_RUNNING__)  
  fparm->AddParam( ip++, kRHCNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kRHCNuisSuperkEnergyScale) );   // RHC SK escale
  for(int i = 0; i < kNumOfRHCSKDetFSIParams; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  //

  #endif

#else
  
  #ifdef __MULTIRING_FHC_RUNNING__ // multiring fhc only 
  
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale) );   // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRing; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + i);    
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #endif
#endif // sample analysis option
///

  if (SKDetFSICovMtrx) delete SKDetFSICovMtrx;

}

#ifdef __NUMUONLY__

//________________________________________________________________________________________
void util::AppendSuperkDetectorFSISystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file


  //const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "util/AppendSuperkDetectorFSISystematicsDisap >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
  //TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathSKDetFSICovMatr.c_str());

  if(!SKDetFSICovMtrx) {
    cerr << "util/AppendSuperkDetectorFSISystematicsDisap >> " << "No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //~ xml->FreeNode(root_node);
  xml->FreeDoc(xmldoc);
  delete xml;

  
  // FSI
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  

  // fhc all + skdetfsi multiring 
  #if defined (__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__) 
    fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
    for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
    for(int i = 0; i < kNumOfSKDetFSIParamsMultiRmu; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    
    // fhc+rhc 
    #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

      // fhc+rhc fully correlated
      #if defined(__FULLY_CORRELATED__) 
      fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
      for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
        // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
        PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
        double err = sqrt(SKDetFSICov[i][i]);
        double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
        fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
      }
      
      
      // NEW FHC+RHC JOINT
      
      // fhc+rhc fully w/ correlations
      #else
      fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
      for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
	double err = sqrt(SKDetFSICov[i][i]);
	double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
	fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
      }      
      fparm->AddParam( ip++, kRHCNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kRHCNuisSuperkEnergyScale));  // RHC SK escale
      for(int i = 0; i < kNumOfRHCSKDetFSIParams1Rmu; i++) {
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + i);
	double err = sqrt(SKDetFSICov[i][i]);
	double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
	fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
      }      
      #endif
      //
    
    // fhc only 
   #elif defined(__FHC_RUNNING__) 
    fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
    for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    
    // rhc only 
   #elif defined(__RHC_RUNNING__) 
    fparm->AddParam( ip++, kRHCNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kRHCNuisSuperkEnergyScale));  // RHC SK escale
    for(int i = 0; i < kNumOfRHCSKDetFSIParams1Rmu; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    //

   #endif

#else
  #ifdef __MULTIRING_FHC_RUNNING__
  // multiring fhc only 
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRmu; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #endif
#endif // sample analysis option
  ///
  if (SKDetFSICovMtrx) delete SKDetFSICovMtrx;
}
#endif // #ifdef __NUMUONLY__

#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendSuperkDetectorFSISystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  
  // Get the covariance matrices & best fits of FSI & BANFF from xml, 
  // to have something sensible to initialise to
  // Parse the input XML configuration file

  //const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "util/AppendSuperkDetectorFSISystematicsApp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
  //TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathSKDetFSICovMatr.c_str());
  
  if(!SKDetFSICovMtrx) {
    cerr << "util/AppendSuperkDetectorFSISystematicsApp >> " << "No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //~ xml->FreeNode(root_node);
  xml->FreeDoc(xmldoc);
  delete xml;
 
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)
  
  // fhc all + skdetfsi multiring
  #if defined (__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__) 
    fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
    for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
    for(int i = 0; i < kNumOfSKDetFSIParamsMultiRe; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    
    // fhc+rhc 
  #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

    // fhc+rhc fully correlated
    #if defined(__FULLY_CORRELATED__)  
      fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
      for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
	double err = sqrt(SKDetFSICov[i][i]);
	double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
	fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
      }
      


      // NEW FHC+RHC JOINT
      
      // fhc+rhc w/ correlations
    #else
      fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
      for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
	double err = sqrt(SKDetFSICov[i][i]);
	double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
	fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
      }
      fparm->AddParam( ip++, kRHCNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kRHCNuisSuperkEnergyScale));  // FHC SK escale
      for(int i = 0; i < kNumOfRHCSKDetFSIParams1Re; i++) {
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + i);
	double err = sqrt(SKDetFSICov[i][i]);
	double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
	fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
      }
    #endif
      //
    
  #elif defined(__FHC_RUNNING__) // fhc only 
    fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
    for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    
  #elif defined(__RHC_RUNNING__) // rhc only 
    fparm->AddParam( ip++, kRHCNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kRHCNuisSuperkEnergyScale));  // FHC SK escale
    for(int i = 0; i < kNumOfRHCSKDetFSIParams1Re; i++) {
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + i);
      double err = sqrt(SKDetFSICov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
    }
    //

#endif

#else
  #ifdef __MULTIRING_FHC_RUNNING__ // multiring fhc only 
  
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRe; i++) {
    // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #endif
#endif // sample analysis option

  if (SKDetFSICovMtrx) delete SKDetFSICovMtrx;
}
#endif // #ifdef __NUEONLY__




//________________________________________________________________________________________
void util::AppendAllCorrelatedBANFFSystematics(FitParamList * fparm, int start)
{

// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // Get the covariance matrices & best fits of FSI & BANFF from xml, 
  // to have something sensible to initialise to
  // Parse the input XML configuration file
  
  //const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "util/AppendAllCorrelatedBANFFSystematics >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathBanffParam.c_str());

  if(!BANFFBestFitValues) {
    cerr << "util/AppendAllCorrelatedBANFFSystematics >> " << " No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());

  if(!BANFFCovMtrx) {
    cerr << "util/AppendAllCorrelatedBANFFSystematics >> " << " No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //~ xml->FreeNode(root_node);
  xml->FreeDoc(xmldoc);
  delete xml;

  //
  // OLD
  //
  //   for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
  // #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__)  
  //     PhysParamId_t p = kListOfBANFFParams[i];
  // #elif defined(__FHC_RUNNING__) 
  //     PhysParamId_t p = kListOfBANFFParams[i];
  // #elif defined(__RHC_RUNNING__)
  //     PhysParamId_t p = kListOfRHCBANFFParams[i];
  // #endif



  // NEW FHC+RHC JOINT

  // Take the number of banff parameters
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)  
  #if defined(__FULLY_CORRELATED__)  
  const int nparams = kNumOfParamsInBANFFfit;
  #else
  const int nparams = kNumOfFHCpRHCParamsInBANFFfit;
  #endif
#elif defined(__FHC_RUNNING__)
  const int nparams = kNumOfParamsInBANFFfit;
#elif defined(__RHC_RUNNING__)
  const int nparams = kNumOfRHCParamsInBANFFfit;
#endif

  // get parameter
  for(int i = 0; i < nparams; i++) {  
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
    #if defined(__FULLY_CORRELATED__)  
    PhysParamId_t p = kListOfBANFFParams[i];
    #else
    PhysParamId_t p = kListOfFHCpRHCBANFFParams[i];
    #endif
#elif defined(__FHC_RUNNING__) 
    PhysParamId_t p = kListOfBANFFParams[i];
#elif defined(__RHC_RUNNING__)
    PhysParamId_t p = kListOfRHCBANFFParams[i];
#endif
    //
    
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    //double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    double lowerr=0;
    if(p==kNuisBANFFShapeCCoth || p==kRHCNuisBANFFShapeCCoth){
      lowerr = bf-(err*SYST_NSIGMA);
    }
    else lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  } // loop over nparams 
  
  if(BANFFBestFitValues ) delete BANFFBestFitValues;
  if(BANFFCovMtrx       ) delete BANFFCovMtrx;
}
//________________________________________________________________________________________


#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendCorrelatedBANFFSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  
  //const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "util/AppendCorrelatedBANFFSystematicsDisap >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathBanffParam.c_str());


  if(!BANFFBestFitValues) {
    cerr << "util/AppendCorrelatedBANFFSystematicsDisap >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());

  if(!BANFFCovMtrx) {
    cerr << "util/AppendCorrelatedBANFFSystematicsDisap >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //~ xml->FreeNode(root_node);
  xml->FreeDoc(xmldoc);
  delete xml;

  // 
  // OLD
  //
  //   for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {    
  // #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) 
  //     PhysParamId_t p = kListOfBANFFParams[i];
  // #elif defined(__FHC_RUNNING__) 
  //     PhysParamId_t p = kListOfBANFFParams[i];
  // #elif defined(__RHC_RUNNING__)
  //     PhysParamId_t p = kListOfRHCBANFFParams[i];
  // #endif
  //

  // NEW FHC+RHC JOINT 
  
  // Take the number of banff parameters
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)  
  #if defined(__FULLY_CORRELATED__)  
  const int nparams = kNumOfParamsInBANFFfit;
  #else
  const int nparams = kNumOfFHCpRHCParamsInBANFFfit;
  #endif
#elif defined(__FHC_RUNNING__)
  const int nparams = kNumOfParamsInBANFFfit;
#elif defined(__RHC_RUNNING__)
  const int nparams = kNumOfRHCParamsInBANFFfit
#endif

   // get parameter
  for(int i = 0; i < nparams; i++) {  
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
    #if defined(__FULLY_CORRELATED__)  
    PhysParamId_t p = kListOfBANFFParams[i];
    #else
    PhysParamId_t p = kListOfFHCpRHCBANFFParams[i];
    #endif
#elif defined(__FHC_RUNNING__) 
    PhysParamId_t p = kListOfBANFFParams[i];
#elif defined(__RHC_RUNNING__)
    PhysParamId_t p = kListOfRHCBANFFParams[i];
#endif
   
    //

    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    //double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    double lowerr=0;
    if(p==kNuisBANFFShapeCCoth || p==kRHCNuisBANFFShapeCCoth){
      lowerr = bf-(err*SYST_NSIGMA);
    }
    else lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));

    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  } // loop over nparams
  
  if(BANFFBestFitValues ) delete BANFFBestFitValues;
  if(BANFFCovMtrx       ) delete BANFFCovMtrx;
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendCorrelatedBANFFSystematicsApp(FitParamList * fparm, int start)
{

// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  
  //const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "util/AppendCorrelatedBANFFSystematicsApp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathBanffParam.c_str());

  if(!BANFFBestFitValues) {
    cerr << "util/AppendCorrelatedBANFFSystematicsApp >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());

  if(!BANFFCovMtrx) {
    cerr << "util/AppendCorrelatedBANFFSystematicsApp >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //~ xml->FreeNode(root_node);
  xml->FreeDoc(xmldoc);
  delete xml;

  // 
  // OLD
  //
  //   for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
  // #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__)  
  //     PhysParamId_t p = kListOfBANFFParams[i];
  // #elif defined(__FHC_RUNNING__) 
  //     PhysParamId_t p = kListOfBANFFParams[i];
  // #elif defined(__RHC_RUNNING__)
  //     PhysParamId_t p = kListOfRHCBANFFParams[i];
  // #endif
 

  // NEW FHC+RHC JOINT

  // Take the number of banff parameters
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)  
  #if defined(__FULLY_CORRELATED__)  
  const int nparams = kNumOfParamsInBANFFfit;
  #else
  const int nparams = kNumOfFHCpRHCParamsInBANFFfit;
  #endif
#elif defined(__FHC_RUNNING__)
  const int nparams = kNumOfParamsInBANFFfit;
#elif defined(__RHC_RUNNING__)
  const int nparams = kNumOfRHCParamsInBANFFfit;
#endif

   // get parameter
  for(int i = 0; i < nparams; i++) {  
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
    #if defined(__FULLY_CORRELATED__)  
    PhysParamId_t p = kListOfBANFFParams[i];
    #else
    PhysParamId_t p = kListOfFHCpRHCBANFFParams[i];
    #endif
#elif defined(__FHC_RUNNING__) 
    PhysParamId_t p = kListOfBANFFParams[i];
#elif defined(__RHC_RUNNING__)
    PhysParamId_t p = kListOfRHCBANFFParams[i];
#endif
    //

    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    //double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    double lowerr=0;
    if(p==kNuisBANFFShapeCCoth || p==kRHCNuisBANFFShapeCCoth){
      lowerr = bf-(err*SYST_NSIGMA);
    }
    else lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));

    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  } // loop over nparams
  
  if(BANFFBestFitValues ) delete BANFFBestFitValues;
  if(BANFFCovMtrx       ) delete BANFFCovMtrx;
}
#endif
