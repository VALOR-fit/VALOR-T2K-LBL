//
// ------------------------------------------------------------------------------------------------
// Generate input for the plot of  contributions from various event modes in the predicted SuperK 
// 1) e-like neutrino and antineutrino spectra  (--compare-by-nu-antinu)
// 2) e-like neutrino +antineutrino  events with Erec below or above 500MeV (--compare-by-ereco)
//
// April 2021 
// Maria Antonova  <maria.antonova \at ific.uv.es>
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
// TODO:check all includes

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
void Get_events( TObjArray *array, valor::FitModelABC* model, double& Nfhc, double& Nrhc, bool compare_by_nu_antinu, bool compare_by_ereco);
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
   bool disable_banff_tuning = false; 
   /*if(arg.Exists("force-untuned") && arg.GetString("force_untuned") != "-") )
   {
      disable_banff_tuning = true;
      LOG("bievent", pNOTICE) << "Forcing paramters to their untuned values";
   };
*/
   // choose the mode we are running in i.e. fhc vs rhc or Erec<500 MeV vs Erec>MeV
   bool compare_by_nu_antinu = false;//
   bool compare_by_ereco  = false;//
   if (arg.GetBool("nu-antinu"))
   {
      LOG("bievent", pINFO) << "Will generate bievent plot  in numu vs antinumu comparison mode";
      compare_by_nu_antinu = true;
   }
   else 
   {
      LOG("bievent", pINFO) << "Will generate bievent plot  in Ereco cut comparison mode";
      compare_by_ereco = true;
   }

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


   //Data BF point to generate systematic contour 
   double bf_dm23sq   = arg.GetDouble("bf-dm23sq");
   double bf_deltacp  = arg.GetDouble("bf-deltacp");
   double bf_mixing23 = arg.GetDouble("bf-mixing23");
   string bf_hier     = arg.GetString("bf-hier");

 // PDG 2016 best fits
   vector<double> ss12 =arg.GetDoubleVec("ss12");
   vector<double> dm12 =arg.GetDoubleVec("dm12");
   vector<double> dm23 = arg.GetDoubleVec("dm23");
   vector<double> ss23 = arg.GetDoubleVec("ss23");


   //Setting global options 
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

   FitParamList * fp = util::GetFitParameterList(hypothesis.c_str(), fplopt.c_str());


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
   

   // set priors so we can read mean and sigma for osc param randomization 
   util::oscillation::SetParameterPriors(arg, *fp); 

   //
   // - Open loop over values of deltaCP (each octant and mass hierarchy)
   //   and fill TGraph with a point of N_fhc_elike Vs N_rhc_elike
   // - Set lines of the TGraph
   // - Calculate N_fhc_elike Vs N_rhc_elike point from data and plot it on top of the MC one
   // 


   TH2D* hEventsSyst = new TH2D("hEventsSyst", "hEventsSyst", 600, 0, 150, 100, 0, 25);
   TH2D* hEvents = new TH2D("hEvents", "hEvents", 600, 0, 150, 100, 0, 25);
   TObjArray* nominals = new TObjArray(0);
   TH1D* hFHCEvents = new TH1D("hEvents_fhc", "hEvents_fhc", 150, 0, 150);
   TH1D* hRHCEvents = new TH1D("hEvents_rhc", "hEvents_rhc", 150, 0, 150);
   


   Rnd::Instance()->SetSeeds(arg);

   for(int e = 0; e < ntoys; ++e)
   {

      if(e % 100 == 0)
      {
         LOG("bievent", pINFO) <<"Throwing toy " << e;
      } 

      model->ResetSystematics(*fp, "all"); 
      

      //NOTE: write nominal for fixed points 
      
      if(e == 0) 
      { 
         for(int idcp=0; idcp < dcp_npt; idcp++)
         {
           for(int iss23 = 0; iss23 < N23; ++iss23)
           {
              for(int imh = 0; imh < NMH; ++imh)
              {
                 double Nfhc = -999;
                 double Nrhc = -999;
      
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
   

                 model->Generate(*fp,"+SB");

                 // Array to store results
                 TObjArray *array = new TObjArray();
                 array->SetOwner(true);
                 LOG("bievent", pDEBUG) << "writting";
                 // Write stuff out
                  model->Write(array,"curr","pdft_rf +SB");
    
                         
                 Get_events(array,model,Nfhc,Nrhc,compare_by_nu_antinu,compare_by_ereco);

                 LOG("bievent", pDEBUG) << "chose dmsq value "<< dm23[imh]<< "cpv = "<< dcp<<" theta23 = "<<ss23[iss23];
                
 
                 nominals->Add(new TMarker(Nfhc, Nrhc, 20)); 
                 
               }
            }
         }
      }//end of the loop over all params to fill nominal vals 

      bool isPhysical = false;
      while(!isPhysical)
      { 
        model->ResetSystematics(*fp, "all"); 
        
        if(e!=0) {
          model->RandomizeSystematics(*fp);
        }

        fp->FromPhysParamList(model->Tweaks(), "current silent");
        if (util::AllWithinLimits(*fp)) isPhysical = true;
      }

      // throwing toys at data bf to generate syst err contour 

      //set osc params to bf
      fp->SetCurrentValue(conv::kPhypDm23sq, bf_dm23sq);
      fp->SetCurrentValue(conv::kPhypCPV, bf_deltacp);
      fp->SetCurrentValue(conv::kPhypMixing23, bf_mixing23);
      Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(bf_hier);
      Set_MH(model); 

      //Generate model
      model->Generate(*fp,"+SB");
       // Array to store results
      TObjArray *array = new TObjArray();
      array->SetOwner(true);
      LOG("bievent", pDEBUG) << "writting";
      // Write stuff out
      double Nfhc = -999;// UWAGA: so far naming stays the same as for standard bi-event plot, just for the sake of time saving. I#ll rename when it will finally work 
      double Nrhc = -999;
      model->Write(array,"curr","pdft_rf +SB");                   
      Get_events(array,model,Nfhc,Nrhc,compare_by_nu_antinu,compare_by_ereco);
      if (e==0){
         nominals->Add(new TMarker(Nfhc, Nrhc, 20)); // NOTE: saving nominal for bf point
      }
      hEventsSyst->Fill(Nfhc, Nrhc); //NOTE: fill histo with data at every toy 
      //

      // Apply statistical fluctuations 
      Nfhc = Rnd::Instance()->RndToyMCstat().Poisson(Nfhc);
      Nrhc = Rnd::Instance()->RndToyMCstat().Poisson(Nrhc);
      hEvents->Fill(Nfhc, Nrhc);
      hFHCEvents->Fill(Nfhc);
      hRHCEvents->Fill(Nrhc);
   }

   TFile file(Form("%s.root", tag.c_str()), "RECREATE");
   hEvents->Write();
   hEventsSyst->Write();
   hFHCEvents->Write();
   hRHCEvents->Write();

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
void Get_events( TObjArray* array, valor::FitModelABC* model, double& Nfhc, double& Nrhc, bool compare_by_nu_antinu,bool compare_by_ereco)
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
      if (compare_by_nu_antinu)
      {
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
         }
      }
      else if (compare_by_ereco){
         if ( s == util::SampleFromString("1Re") || s == util::SampleFromString("MultiRe"))
         {
            // Nfhc -> ereco>500, Nrhc Ereco<500
            TString name = TString::Format("curr_pdft_sb_rf_%s_%s", this_model->Name().c_str(), util::AsString(ks).c_str());
            TH1 * curr_pdf = (TH1*) array->FindObject(name.Data());
            if (!curr_pdf) 
            {
               LOG("bievent", pERROR) << "Couldn't find pdf histogram: " << name;
                  util::Die(-1);
            }
            int cut_bin  = curr_pdf->GetXaxis()->FindBin(0.540);
            LOG("bievent", pDEBUG) <<"x min "<<curr_pdf->GetXaxis()->GetBinCenter(1)<< " x max "<<curr_pdf->GetXaxis()->GetBinCenter(24);

            if (Nfhc != -999 && Nrhc!=-999)
            {
               Nfhc = Nfhc + curr_pdf->Integral(cut_bin+1, curr_pdf->GetNbinsX());
               Nrhc = Nrhc + curr_pdf->Integral(0, cut_bin);
               LOG("bievent", pDEBUG) << " Nevts tot " << curr_pdf->Integral()<< " and summ of both "<< (curr_pdf->Integral(cut_bin, curr_pdf->GetNbinsX()+1)) + (curr_pdf->Integral(1, cut_bin-1)) ;
            }

            else 
            {
              Nfhc = curr_pdf->Integral(cut_bin+1, curr_pdf->GetNbinsX()); 
              Nrhc = curr_pdf->Integral(0, cut_bin); 
              LOG("bievent", pDEBUG) << " Nevts tot " << curr_pdf->Integral()<< " and summ of both "<< (curr_pdf->Integral(cut_bin, curr_pdf->GetNbinsX()+1)) + (curr_pdf->Integral(1, cut_bin-1)) ;
            }
         }
      }
      else{
           LOG("bievent", pERROR) << "No plot generation mode was selected. Nothing to be done here.";
            util::Die(-1);
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



