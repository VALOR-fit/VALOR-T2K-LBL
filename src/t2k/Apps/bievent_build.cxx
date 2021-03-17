//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip>
//root 

#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMarker.h"
// ToDo:check all includes
// fmwk
#include "valor/Conventions/Opt.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Utilities/UtlSystem.h" //for config path 
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/FitCore/CorrelatedSystematics.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/CompoundFitModel.h"
#include "valor/FitCore/UtlOscillation.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/FitCore/UtlFit.h"


#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"

#include "valor/Utilities/Env.h"
#include "valor/Utilities/VThreadPool.h"

#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/SingleSampleKineBinning.h"
#include "valor/Utilities/OscParamValues.h"

#include "valor/Ntuple/StdFitNtp.h"
#include "valor/KineBinning/VALORBinning.h"


//temp 
#include "valor/Conventions/ParameterPool.h"
#include "valor/FitCore/PhysParamList.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

using std::string;
using std::cout;
using std::cerr;
using std::endl;
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");
void Bievent_build(CmdLnArgList& arg);

void Set_MH(valor::FitModelABC* model);
void Get_events( TObjArray *array, valor::FitModelABC* model, double& Nfhc, double& Nrhc);
void RandomizeOsc(valor::FitParamList *fp);
bool OscParamsPhysical(double mix12val, double mix13val, double dm12val);

int main(int argc, char** argv)
{
   gROOT->SetBatch(1);

   CmdLnArgList arg;
   arg.Read(argc,argv,
      "$VALOR_EXPT/config/apps/bievent.arg"
      ":$VALOR_EXPT/config/apps/seeds.arg");
   LOG("bievent", pNOTICE) << arg;
   Bievent_build(arg);
   LOG("bievent", pNOTICE) << "Done!";
  
   return 0;
}


void Bievent_build(CmdLnArgList& arg)
{
   Messenger::Instance()->Configure(util::CompleteConfigPath(arg.GetString("mesg-thresholds")));
   // read command line parameters 
   int ntoys = arg.GetInt("ntoys");
   string tag = arg.GetString("tag");

   std::string config_path = util::CompleteConfigPath(arg.GetString("config"));
  
   string hypothesis       = arg.GetString("physics-hypothesis");
   string fplopt           = arg.GetString("fit-param-opt");
   bool disable_banff_tuning = false; //UWAGA: !!!!!!!! this is not meant to be this way. It is just to make sure systematics is actually applied 
   //if(arg.Exists("force-untuned"))
   //{
   //   disable_banff_tuning = true;
   //   LOG("bievent", pNOTICE) << "Forcing paramters to their untuned values";
   //};

   double deltaCP_first = -TMath::Pi();
   double deltaCP_last = +TMath::Pi();
   const int dcp_npt = arg.GetInt("dcp-npt"); 
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt-1);
   const int N23 = arg.GetInt("N23");
   const int NMH = arg.GetInt("NMH");

   // Best fit values NH
   double mixing12 = arg.GetDouble("input-mixing12");
   double mixing13 = arg.GetDouble("input-mixing13");
   double mixing23 = arg.GetDouble("input-mixing23");
   double deltacp  = arg.GetDouble("input-deltacp");
   double dm12sq   = arg.GetDouble("input-dm12sq");
   double dm23sq   = arg.GetDouble("input-dm23sq");
   double betabar  = arg.GetDouble("input-betabar");


 // PDG 2016 best fits
   vector<double> ss12 =arg.GetDoubleVec("ss12");
   vector<double> dm12 =arg.GetDoubleVec("dm12");
   // 2018 best fit values from 2D fits (including IH shift for dmsq23)
   vector<double> ss13 = arg.GetDoubleVec("ss13");
   vector<double> dm23 = arg.GetDoubleVec("dm23");
   //dm23[1] = -0.00246;
   cout<<"dm23 for NH is "<< dm23[0]<< " and for IH "<< dm23[1];
      //dm23[2] = dm23[2] - dm12[2]; //IH shift for dmsq23 // not shure if we stil need this. Do we? 

   // Alternative values to use for sinsqtheta23
   vector<double> ss23 = arg.GetDoubleVec("ss23");


   //Setting global options 
   //Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);
   Opt::Instance()->gOscMixingParamChoice = ( (arg.GetBool("sinsq2theta"))  ? conv::kNuMixSinsq2Theta : conv::kNuMixSinsqTheta);
   Opt::Instance()->gOscOctantChoice      = ( (arg.GetBool("second-octant"))? conv::kNuOctantTwo      : conv::kNuOctantOne);
   conv::NuOscDmsqScheme_t dmsq_scheme = util::DmsqSchemeFromString(arg.GetString("dmsq-scheme"));
   
   if(dmsq_scheme == kNuDmsqUndefined) 
   {
    LOG("bievent", pFATAL) << "Could not parse the input Dmsq scheme: [" << arg.GetString("dmsq-scheme") << "]";
    util::Die(-1);
   }
    Opt::Instance()->gOscDmsqScheme = dmsq_scheme;
  if(Opt::Instance()->Undefined())
   {
      LOG("bievent", pFATAL) << "Global options not defined correctly. Exiting...";
      util::Die(-1);
   }

   // Print job options
   cout << "\n[Job inputs]" << endl;
   cout << "- Reading analysis configuration from: " << arg.GetString("config") << endl;
   cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
   cout << "- Physics hypothesis: " << hypothesis << endl;
   cout << *Opt::Instance() << endl;

  

   CompoundFitElements * elements = new CompoundFitElements();
   elements->Read(config_path.c_str());
   FitModelABC * model = elements->FitModel();

   // Set and update fpl 
   //const CorrelatedSystematics * syst = cm->GetSystematics(); //  probably we need this line 
   FitParamList * fp = util::GetFitParameterList(hypothesis.c_str(), fplopt.c_str());
   LOG("bievent", pDEBUG) << "is phys before updtae " << util::AllWithinLimits(*fp);
   ////DBUG CODE -----------
   //for (int iparam = 0 ; iparam < fp->Size(); iparam++) {
   // PhysParamId_t p = fp->WhichParam(iparam);
   // if(!ParameterPool::Instance()->IsParamInPool(p)) {
   //   LOG("bievent", pFATAL) << "The parameter \'" << util::NameRegistry::Instance()->ParamAsString(p)
   //                            << "\' was in the FPL but not the ParameterPool";
   //   util::Die(-1);
   // }
   // const Parameter & param = ParameterPool::Instance()->GetParam(p);
//
   // double value = fp->CurrValueByIndex(iparam);
   // if (not param.IsWithinLimits(value)) {
   //   LOG("bievent", pDEBUG) << "prameter  " << util::NameRegistry::Instance()->ParamAsString(p) << " is not within PL with value " << value;
   //       }
   //}

   //DBUG CODE -----------

   //Update fpl
   if (fp->Exists(conv::kPhypMixing12)) { fp->Update(conv::kPhypMixing12,mixing12); }// set  of the parameters for 3flav osc 
   if (fp->Exists(conv::kPhypMixing13)) { fp->Update(conv::kPhypMixing13,mixing13); }
   if (fp->Exists(conv::kPhypMixing23)) { fp->Update(conv::kPhypMixing23,mixing23); }
   if (fp->Exists(conv::kPhypDm12sq  )) { fp->Update(conv::kPhypDm12sq  ,dm12sq  ); }
   if (fp->Exists(conv::kPhypDm23sq  )) { fp->Update(conv::kPhypDm23sq  ,dm23sq  ); }
   if (fp->Exists(conv::kPhypCPV     )) { fp->Update(conv::kPhypCPV     ,deltacp ); }

   if (fp->Exists(conv::kPhypMixing12Bar)) { fp->Update(conv::kPhypMixing12Bar,mixing12); } // additional params for neutrino oscillations 
   if (fp->Exists(conv::kPhypMixing13Bar)) { fp->Update(conv::kPhypMixing13Bar,mixing13); } //(3-flavour framework) - with separate mixings 
   if (fp->Exists(conv::kPhypMixing23Bar)) { fp->Update(conv::kPhypMixing23Bar,mixing23); } // for neutrino and antineutrino
   if (fp->Exists(conv::kPhypDm12sqBar  )) { fp->Update(conv::kPhypDm12sqBar  ,dm12sq  ); }
   if (fp->Exists(conv::kPhypDm23sqBar  )) { fp->Update(conv::kPhypDm23sqBar  ,dm23sq  ); }
   if (fp->Exists(conv::kPhypCPVBar     )) { fp->Update(conv::kPhypCPVBar     ,deltacp ); }

   if (fp->Exists(conv::kPhypBetaBar    )) { fp->Update(conv::kPhypBetaBar    ,betabar); }  
   
   //DBUG CODE -----------
   LOG("bievent", pDEBUG) << "is phys after updtae " << util::AllWithinLimits(*fp);

   for (int iparam = 0 ; iparam < fp->Size(); iparam++) {
    PhysParamId_t p = fp->WhichParam(iparam);
    if(!ParameterPool::Instance()->IsParamInPool(p)) {
      LOG("bievent", pFATAL) << "The parameter \'" << util::NameRegistry::Instance()->ParamAsString(p)
                               << "\' was in the FPL but not the ParameterPool";
      util::Die(-1);
    }
    const Parameter & param = ParameterPool::Instance()->GetParam(p);

    double value = fp->CurrValueByIndex(iparam);
    if (not param.IsWithinLimits(value)) {
      LOG("bievent", pDEBUG) << "prameter  " << util::NameRegistry::Instance()->ParamAsString(p) << " is not within PL with value " << value;
          }
  }
   //DBUG CODE -----------
   // set priors so we can read mean and sigma for osc param randomization 
   util::oscillation::SetParameterPriors(arg, *fp); 

   //
   // - Open loop over values of deltaCP (each octant and mass hierarchy)
   //   and fill TGraph with a point of N_fhc_elike Vs N_rhc_elike
   // - Set lines of the TGraph
   // - Calculate N_fhc_elike Vs N_rhc_elike point from data and plot it on top of the MC one
   // 


   TH2D* hEventsSyst[dcp_npt][N23][NMH];
   TH2D* hEvents[dcp_npt][N23][NMH];
   TObjArray* nominals = new TObjArray(0);
   TH1D* hFHCEvents[dcp_npt][N23][NMH];
   TH1D* hRHCEvents[dcp_npt][N23][NMH];
   for(int idcp=0; idcp < dcp_npt; idcp++)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            const char* name_syst = Form("hEventsSyst_%d_%d_%d", idcp, iss23, imh);
            hEventsSyst[idcp][iss23][imh] = new TH2D(name_syst, name_syst, 600, 0, 150, 100, 0, 25);
            const char* name = Form("hEvents_%d_%d_%d", idcp, iss23, imh);
            hEvents[idcp][iss23][imh] = new TH2D(name, name, 150, 0, 150, 25, 0, 25);
            const char* fhcname = Form("hEvents_fhc_%d_%d_%d", idcp, iss23, imh);
            hFHCEvents[idcp][iss23][imh] = new TH1D(fhcname, fhcname, 150, 0, 150);
            const char* rhcname = Form("hEvents_rhc_%d_%d_%d", idcp, iss23, imh);
            hRHCEvents[idcp][iss23][imh] = new TH1D(rhcname, rhcname, 25, 0, 25);
         }
      }
   }
   for(int e = 0; e < ntoys; ++e)
   {

      if(e % 100 == 0)
      {
         cout << "Throwing toy " << e << endl;
      } 
      LOG("bievent", pDEBUG) << "rnd seeds ";
      Rnd::Instance()->SetSeeds(arg);
      LOG("bievent", pDEBUG) << "rset is isPhysical ";
      bool isPhysical = false;
      while(!isPhysical)
      {
         LOG("bievent", pDEBUG) << "enter loop and reset syst";
         LOG("bievent", pDEBUG) << "is phys " << util::AllWithinLimits(*fp);
        if(disable_banff_tuning) // force untuned from sigma var. Maybe rename ? 
         {
            cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
            model->ForceSystematicsToUntuned(*fp);
         }
         else
         {
            model->ResetSystematics(*fp, "all");  // this causes parameters be outside the limits 
         }        
         if(e > 0)
         {  // Generate one nominal set and use random throws for the rest
             LOG("bievent", pDEBUG) << "Generate one nominal set and use random throws for the rest";
            model->RandomizeSystematics(*fp); 
            LOG("bievent", pDEBUG) << "RandomizeOsc";
            RandomizeOsc(fp);
         }
         LOG("bievent", pDEBUG) << "is phys " << util::AllWithinLimits(*fp);
         LOG("bievent", pDEBUG) << "tweak fpl";
         fp->FromPhysParamList(model->Tweaks(), "current silent");
         LOG("bievent", pDEBUG) << "is phys " << util::AllWithinLimits(*fp);
        
         isPhysical = true;//util::AllWithinLimits(*fp); // think on how to fix this check 
      }

      for(int idcp=0; idcp < dcp_npt; idcp++)
      {
         for(int iss23 = 0; iss23 < N23; ++iss23)
         {
            for(int imh = 0; imh < NMH; ++imh)
            {
               double Nfhc = -999;
               double Nrhc = -999;
               double Ncc1pi = -999;
    
               double dcp = deltaCP_first + idcp * deltaCP_width;
       
               if(imh == 0)
               {
                  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
                  Set_MH(model); 
               }
               else
               {
                  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
                  Set_MH(model);
               }

               fp->SetCurrentValue(conv::kPhypDm23sq, dm23[imh]);
               fp->SetCurrentValue(conv::kPhypCPV, dcp);
               fp->SetCurrentValue(conv::kPhypMixing23, ss23[iss23]);
   
               LOG("bievent", pDEBUG) << "about to getnerate PDF";
               // Generate predictions 
               model->Generate(*fp,"+SB");

               // Array to store results
               TObjArray *array = new TObjArray();
               array->SetOwner(true);
               LOG("bievent", pDEBUG) << "writting";
               // Write stuff out
                model->Write(array,"curr","pdft_rf +SB");
    
                       
               Get_events(array,model,Nfhc,Nrhc);

               hEventsSyst[idcp][iss23][imh]->Fill(Nfhc, Nrhc); //UWAGA: not quite sure why we are creating syst contour without applyig stat fluctuations. 

               LOG("bievent", pDEBUG) << "chose dmsq value "<< dm23[imh]<< "cpv = "<< dcp<<" theta23 = "<<ss23[iss23];
              
               if(e == 0) 
               {  
                 nominals->Add(new TMarker(Nfhc, Nrhc, 20)); // NOTE: saving nominal values for points at the plot
               }

                 // Apply statistical fluctuations 
               Nfhc = Rnd::Instance()->RndToyMCstat().Poisson(Nfhc);
               Nrhc = Rnd::Instance()->RndToyMCstat().Poisson(Nrhc);
               hEvents[idcp][iss23][imh]->Fill(Nfhc, Nrhc);

               hFHCEvents[idcp][iss23][imh]->Fill(Nfhc);
               hRHCEvents[idcp][iss23][imh]->Fill(Nrhc);
            }
         }
      }
   }

   TFile file(Form("%s.root", tag.c_str()), "RECREATE");
   for(int idcp=0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            hEvents[idcp][iss23][imh]->Write();
            hEventsSyst[idcp][iss23][imh]->Write();
            hFHCEvents[idcp][iss23][imh]->Write();
            hRHCEvents[idcp][iss23][imh]->Write();
         }
      }
   }
   TFile filenom(Form("%snominal.root", tag.c_str()), "RECREATE");
   nominals->Write("Nominals", TObject::kSingleKey);
   filenom.Close();
   file.Close();

   return;
}
//________________________________________________________________________________________
void Set_MH(valor::FitModelABC* model)
{  // NOTE: accesses MH @NuOsc3Flavour through an instance in GSSFM for each submodel 
   //       this is a work around as this instance is private 
    CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
   assert(cm);
   for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) 
   {
      GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
      this_model->UpdateMassHierarchyFromOpt();
   }
}
//________________________________________________________________________________________
void Get_events( TObjArray* array, valor::FitModelABC* model, double& Nfhc, double& Nrhc)
{
   CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
   assert(cm); 
   for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) 
   {
      GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
      BeamConf_t b = this_model->GetBeam();
      Sample_t   s = this_model->GetSample();
      const SingleSampleKineBinning & binning = *this_model->GetBinning();
      const VALORBinning & reco_smp = *(binning.RecoSmpBinning());
      const std::vector<KineVar_t>  & reco_vars = reco_smp.KineVars();
      conv::KineSpace_t ks = util::KineSpaceFromKineVars(reco_vars); 
      if ((b == util::BeamConfFromString("JPARC_FHC")) && (s == util::SampleFromString("1Re") || s == util::SampleFromString("MultiRe")))
      {
         TString name = TString::Format("curr_pdft_sb_rf_%s_%s", this_model->Name().c_str(), util::AsString(ks).c_str());
         TH1 * curr_pdf = (TH1*) array->FindObject(name.Data());
         if (!curr_pdf) 
         {
            LOG("bievent", pERROR) << "Couldn't find pdf histogram: " << name;
               util::Die(-1);
         }
         if (Nfhc != -999)
         {
            Nfhc = Nfhc + curr_pdf->Integral();
         }
         else 
         {
           Nfhc = curr_pdf->Integral(); 
         }
         LOG("bievent", pDEBUG) << "Nfhc  " << Nfhc;
      }

      if ((b == util::BeamConfFromString("JPARC_RHC")) && (s == util::SampleFromString("1Re")))
      {
         TString name = TString::Format("curr_pdft_sb_rf_%s_%s", this_model->Name().c_str(), util::AsString(ks).c_str());
         TH1 * curr_pdf = (TH1*) array->FindObject(name.Data());
         if (!curr_pdf) 
         {
            LOG("bievent", pERROR) << "Couldn't find pdf histogram: " << name;
               util::Die(-1);
         }
           Nrhc = curr_pdf->Integral();
           LOG("bievent", pDEBUG) << "Nrhc  " << Nrhc;
      }
   }// loop models 
   // TODO: add checks  if n events is not -999
  
}
//________________________________________________________________________________________
bool OscParamsPhysical(double mix12val, double mix13val, double dm12val)
{
   if(mix12val < 0 || mix12val > 1)
   {
      return false;
   }
   if(mix13val < 0 || mix13val > 1)
   {
      return false;
   }
   if(dm12val < kMinFitDm12sq || dm12val > kMaxFitDm12sq)
   {
      return false;
   }

   return true;
}
//________________________________________________________________________________________
void RandomizeOsc(valor::FitParamList *fp)

{
     const util::NameRegistry * nr = util::NameRegistry::Instance();
   string prior_scheme_name_mixing13  = fp->GetPriorSchemeName(nr->PhysParamFromString("mixing13"));
   string prior_scheme_name_mixing12  = fp->GetPriorSchemeName(nr->PhysParamFromString("mixing12"));
   string target_scheme_name_mixing13 = ""; // Leave this as an empty string to allow util::ConvertOscParamToUsedScheme to set the target scheme according to the values set in Opt::Instance
   string target_scheme_name_mixing12 = ""; // Leave this as an empty string to allow util::ConvertOscParamToUsedScheme to set the target scheme according to the values set in Opt::Instance

     
   for( ; ; )
   { 
    
      double mix12val = Rnd::Instance()->RndToyMCosc().Gaus(fp->GetMeanGaus(nr->PhysParamFromString("mixing12")), 
         fp->GetSigmaGaus(nr->PhysParamFromString("mixing12")));
      
      double mix13val = Rnd::Instance()->RndToyMCosc().Gaus(fp->GetMeanGaus(nr->PhysParamFromString("mixing13")), 
         fp->GetSigmaGaus(nr->PhysParamFromString("mixing13")));
      
      double dm12val = Rnd::Instance()->RndToyMCosc().Gaus(fp->GetMeanGaus(nr->PhysParamFromString("dm12sq")), 
         fp->GetSigmaGaus(nr->PhysParamFromString("dm12sq")));

      LOG("bievent", pDEBUG) << "mixing 12 before correcting to the mixing scheme "<< mix12val;
      LOG("bievent", pDEBUG) << "mixing 13 before correcting to the mixing scheme "<< mix13val;
      LOG("bievent", pDEBUG) << "dmsq 12   before correcting to the mixing scheme "<< dm12val;  

      mix12val = util::ConvertOscParamToUsedScheme(nr->PhysParamFromString("mixing12"), mix12val, prior_scheme_name_mixing12, target_scheme_name_mixing12); //util::FromSinsq2thToSinsqth(mix12val);
      mix13val = util::ConvertOscParamToUsedScheme(nr->PhysParamFromString("mixing13"), mix13val, prior_scheme_name_mixing13, target_scheme_name_mixing13);//util::FromSinsq2thToSinsqth(mix13val);


      LOG("bievent", pDEBUG) << "mixing 12 after correcting to the mixing scheme "<< mix12val;
      LOG("bievent", pDEBUG) << "mixing 13 after correcting to the mixing scheme "<< mix13val;
      LOG("bievent", pDEBUG) << "dmsq 12 after correcting to the mixing scheme "<< dm12val;     
      //Update FPL
      fp->Update(nr->PhysParamFromString("mixing12"),mix12val);
      fp->Update(nr->PhysParamFromString("mixing13"),mix13val);
      fp->Update(nr->PhysParamFromString("dm12sq"),dm12val);

      LOG("bievent", pDEBUG) << "randomizing osc";

      // check if all parameters are within the  physics limits 
      if(OscParamsPhysical(mix12val, mix13val, dm12val))
      {
         
        break;
      }
      
   }
}


