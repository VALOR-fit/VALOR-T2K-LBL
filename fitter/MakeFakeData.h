//____________________________________________________________________________
/*!
 
 \class     valor::utils
 
 \brief     Utilities for formatting and style
 
 \author    Raj Shah <raj.shah \at stfc.ac.uk>
 University of Oxford & STFC Rutherford Appleton Laboratory
 
 \created   Oct 10, 2015
 
 \copyright (c) 2010-2015, The VALOR group.
 
 */
//____________________________________________________________________________

#ifndef _valor_MakeFakeData_
#define _valor_MakeFakeData_

#include <iostream>
#include <string>
#include "Conventions.h"
#include "TH3.h"
#include "TH2.h"
#include "TGraph.h"

#include <TXMLEngine.h>


using std::vector;

namespace valor {
  namespace util {
    
    class  UtlMakeFakeData {
      
    public:
      static  UtlMakeFakeData * Instance();
      bool    Read(TXMLEngine* xml, XMLNodePointer_t root_node);
      double  GetWeight();
      void    SetEventInfo(float _fBrPtrue[10],float _fBrDirNu[10][3],float _fBrDir[10][7][3],float _Etrue, conv::SuperkMode_t _m, double _ERec, conv::SuperkSample_t _s, float _fBrAmome[10], float _fBrAmomm[10]);
      void    SetReweightStyle(string _s = "error");
      double  GetEreco();
      bool    IsActive(){return fActive;};
      void    Activate(){ fActive = true;};
      void    Deactivate(){ fActive = false;};
      string  GetReweightStyle(){return fReweightString;};
      
      enum ReweightType
      {
        FDUndefined   = -1,
        Martini     = 0,
        MECPDD2D    = 1,
        MECPDD3D    = 2,
        Nieves1P1H  = 3,
        MECPDD3DFlatFlux = 4,
	// new 2017
	2p2hPDD2D_2017 = 5,
	2p2hnonPDD2D_2017 = 6
	//
      };
      
    private:
      UtlMakeFakeData();
      UtlMakeFakeData(const UtlMakeFakeData &) {}
      ~UtlMakeFakeData() {}
      
      TGraph* Martini_to_Nieves_MEC_ratio_Nu;
      TGraph* Martini_to_Nieves_MEC_ratio_NuBar;
      
      TH2D* NievesMomentumOffset;
      
      TH2D* PDD_ratio_nu_2D;
      TH2D* PDD_ratio_nubar_2D;
      TH3D* PDD_ratio_nu_3D;
      TH3D* PDD_ratio_nubar_3D;

      TH3D* PDD_ratio_numu;
      TH3D* PDD_ratio_numubar;
      TH3D* PDD_ratio_nue;
      TH3D* PDD_ratio_nuebar;

      // new 2017
      TH2D* 2p2hPDD_2017;
      TH2D* 2p2hnonPDD_2017;
      double Weight2p2hPDD_2017();
      double Weight2p2hnonPDD_2017();
      //

      double WeightMECPDD3DFlatFlux();
      double WeightMECPDD3D();
      double WeightMECPDD2D();
      double WeightNieves1P1H();
      double WeightMartini();
      double ErecoFunc_Nu(double Mom,double costh,double mass_lep);
      double ErecoFunc_NuBar(double Mom,double costh,double mass_lep);
      double ErecoFuncCCOnePi_Nu(double Mom,double costh,double mass_lep);
      void   Init(void);
      
      double EReco1P1H();
      static UtlMakeFakeData * fInstance;
      bool   fActive;
      string fReweightString;
      float      fBrAmomm[10];                  ///< branch var for reading input MC mini-DST (1 mu-like ring)
      float      fBrAmome[10];                  ///< branch var for reading input MC mini-DST (1 e-like ring)
      
      //Event info
      
      float Etrue;
      double ERec;
      conv::SuperkMode_t m;
      conv::SuperkSample_t s;
      float fBrDirNu[10][3];
      float fBrDir[10][7][3];
      float fBrPtrue[10];
      
      
      ReweightType reweight_option;
      
      struct Cleaner {
        void DummyMethodAndSilentCompiler() { }
        ~Cleaner() {
          if (UtlMakeFakeData::fInstance !=0) {
            delete UtlMakeFakeData::fInstance;
            UtlMakeFakeData::fInstance = 0;
          }
        }
      };
      friend struct Cleaner;
      
      //      ClassDef(UtlMakeFakeData,1);
    }; //
  }//util
} // valor

#endif



