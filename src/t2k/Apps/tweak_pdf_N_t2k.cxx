//____________________________________________________________________________
/*!

\app       valor-tweak_pdf_t2k

\brief     Vary the MC spectra by throwing sets of systematics

\options     --output              : Where to save tweaked histograms to. [def: tweaks_[PARAMS].root]
             --ntoys               : Number of toy-MC experiments [def: 200]
             --physics-hypothesis  : physics hypothesis [def: "3f_osc"]
             --fit-param-opt       : Options used when setting up the FitParamList
                                     [def: "NONE"]
             --config              : Analysis configuration file [required]
             --mesg-thresholds     : Message thresholds [def: "verbosity.xml"]
             --enu-reco-max        : Maximum Ereco to show on the plot [def: 10]
             --params              : Which parameters to tweak [def: allsyst]

             --input-mixing12,etc  : Option to give input values to osc params (in .arg)

             --load-prior-from-file : Option for file to give constrained prior matrix [ use default prior ]
             --prior-fit-number     : Option to specify the entry in the tree for above file [ 0 ]

             [--nthreads ]         : Number of CPU threads to use [ std::thread::hardware_concurrency ]


\author    Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
           University of Liverpool & STFC Rutherford Appleton Laboratory

           Steve Dennis <s.r.dennis \at liverpool.ac.uk>
           Jaiden Parlone <j.parlone \at liverpool.ac.uk>
           University of Liverpool

\created   October 2023

\copyright Copyright (c) 2010-2020, The VALOR group. All rights reserved

*/
//____________________________________________________________________________

// System
#include <sstream>
#include <fstream>
#include <iomanip>
// ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#include "TPRegexp.h"
// VALOR/SDK
#include "valor/Conventions/Opt.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/FitCore/CorrelatedSystematics.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/FitCore/Likelihood.h"
#include "valor/FitCore/FitterABC.h"
#include "valor/FitCore/UtlFit.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/UtlOscillation.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/UtlForMacros.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/VThreadPool.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/SingleSampleKineBinning.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"
#include "valor/KineBinning/VALORBinning.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

const double enu_reco_max = 10;
const KineVar_t kv = conv::kvEreco;


int main(int argc, char ** argv)
{
  valor::util::PrintBanner(Env::Instance()->Banner(),500);

  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_FMWK/config/apps/make_table_systematic_effects.arg:$VALOR_FMWK/config/apps/seeds.arg:$VALOR_FMWK/config/apps/oscparams.arg");
  LOG("valor", pNOTICE) << arg;

  // Read the command line
  int ntoys                 = arg.GetInt("ntoys");
  int n_threads             = arg.GetInt("nthreads");
  string hypname            = arg.GetString("physics-hypothesis");
  string fplopt             = arg.GetString("fit-param-opt");
  string output_rootfile    = arg.GetString("output");
  string params_to_tweak    = arg.GetString("params");
  string params_to_exclude  = arg.GetString("params-exclude");
  bool rethrow_unphysical = arg.GetBool("rethrow-unphysical");

  std::string mesg_thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(mesg_thresholds);

  LOG("valor",pINFO) << "Will rethrow unphysical toys";

  // Set up the thread pool
  if(n_threads>0) {//without n-cores specified, use as many as the host has
    util::VThreadPool::Instance()->SetNThreads(n_threads);
  }

  // Set the random seeds
  Rnd::Instance()->SetSeeds(arg);

  //change if something different to default (i.e. normal, sinsqtheta, first octant and dmsq12abs23 (NH) / dmsq12abs31 (IH))
  Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);
  Opt::Instance()->gOscMixingParamChoice = ( (arg.GetBool("sinsq2theta"))  ? conv::kNuMixSinsq2Theta : conv::kNuMixSinsqTheta);
  Opt::Instance()->gOscOctantChoice      = ( (arg.GetBool("second-octant"))? conv::kNuOctantTwo      : conv::kNuOctantOne);

  conv::NuOscDmsqScheme_t dmsq_scheme = util::DmsqSchemeFromString(arg.GetString("dmsq-scheme"));
  if(dmsq_scheme == kNuDmsqUndefined) {
    LOG("valor", pFATAL) << "Could not parse the input Dmsq scheme: [" << arg.GetString("dmsq-scheme") << "]";
    util::Die(-1);
  }
  Opt::Instance()->gOscDmsqScheme = dmsq_scheme;

  // Config
  std::string config_path = util::CompleteConfigPath(arg.GetString("config"));

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config_path << endl;
  cout << "- Physics hypothesis used : " << hypname << endl;

  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());
  FitModelABC * model = elements->FitModel();
  CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
  assert(cm);
  const util::NameRegistry * nr = util::NameRegistry::Instance();

  // Arrays to store the spectra for the default, +1 sigma and -1 sigma variations
  vector<TH1D *> nom_pdf_werr;
  vector<TH1D *> nom_pdf     ;
  vector<TH1D *> err_p1      ;
  vector<TH1D *> err_m1      ;

  //conv::PhysParamId_t debug_param = nr->PhysParamFromString("2p2h_Edep_lowEnu"); //DEBUG

  vector<double> max(cm->GetNSubModels());
  vector<vector<TH1D *> > twk_pdf(cm->GetNSubModels());

  // Array to store the number of events in each sample for each toy
  double events [cm->GetNSubModels()][ntoys];

  std::vector<PhysParamId_t> param_set = util::GetParameterSet(params_to_tweak, true);   // The 2nd arg makes it only return exact matches. This prevents e.g. askign for 'mixing12' and getting 'mixing12_bar' as well as 'mixing12'
  std::vector<PhysParamId_t> params_to_remove;
  if(params_to_exclude != "NONE") {
    params_to_remove = util::GetParameterSet(params_to_exclude, true);
  }
  if(params_to_remove.size() > 0) {
    util::RemoveItemsFromVector<PhysParamId_t>(param_set, params_to_remove);
  }
  assert(param_set.size());

  std::map<PhysParamId_t,bool> syst_param_map;
  std::vector<PhysParamId_t> osc_params;
  for (unsigned int ip = 0 ; ip < param_set.size() ; ip++) {
    PhysParamId_t p = param_set[ip];
    if(!util::IsSystematic(p)) { osc_params.push_back(p); }
    else                    { syst_param_map[p] = true; }
  }

  cout << "Randomising osc params:" << endl;
  for(PhysParamId_t p : osc_params) {
    cout << "\t" << util::NameRegistry::Instance()->ParamAsStringVerySimple(p) << endl;
  }
  cout << "Randomising syst params:" << endl;
  for(auto const& pair : syst_param_map) {
    cout << "\t" << util::NameRegistry::Instance()->ParamAsStringVerySimple(pair.first) << endl;
  }

  FitParamList * fp = util::GetFitParameterList(hypname.c_str(), fplopt.c_str());
  int i = 0;

  util::oscillation::SetParameterPriors(arg, *fp);

  cout << "Randomising syst params:" << endl;
  for(auto const& pair : syst_param_map) {
      cout << "\t" << util::NameRegistry::Instance()->ParamAsStringVerySimple(pair.first) << endl;
  }


  // Generate pdfs
  unsigned int n_rethrew = 0;
  for(int itoy = -1; itoy < ntoys; itoy++) {
      // Randomise the systematics in the model, and copy them to the FPL, toy -1 should be nominal.

      model->ResetSystematics(*fp, "all");
      if (itoy == -1) {
          fp->Print();
      }
      if (itoy > -1) {

          model->RandomizeSystematics(syst_param_map, *fp);

      }


    //std::cout << "At DEBUG " << i++ <<  " Value of " << util::NameRegistry::Instance()->ParamAsString(debug_param) << " is " << fp->CurrValue(debug_param) << std::endl; //DEBUG
    fp->FromPhysParamList(model->Tweaks(),"current-silent");
    //std::cout << "At DEBUG " << i++ <<  " Value of " << util::NameRegistry::Instance()->ParamAsString(debug_param) << " is " << fp->CurrValue(debug_param) << std::endl; //DEBUG
    // Set input parameters
    RgKeyList input_args = arg.FindKeys("input-");
    TPRegexp input_match("input-(\\w*)");
    for (RgKey argname : input_args) { // RgKey is just a typename for std::string right now
      TObjArray * matches = input_match.MatchS(argname.c_str());

      if (matches->GetEntries() != 2) {
        LOG("valor",pFATAL) << "Cannot parse argument "<<argname<<" to get parameter name";
        util::Die(-1);
      }
      TObjString * tos_param_name = dynamic_cast<TObjString*>(matches->At(1));
      std::string param_name(tos_param_name->String().Data());
      delete matches;

      conv::PhysParamId_t param = nr->PhysParamFromString(param_name);
      if (param == conv::kPhypUndefined) {
        LOG("valor",pFATAL) << "Cannot find parameter to fix from argument string "<<argname;
        util::Die(-1);
      }

      //if (debug_param == conv::kPhypUndefined) {
      //    LOG("valor",pFATAL) << "Cannot find parameter to fix from argument string "<< util::NameRegistry::Instance()->ParamAsString(debug_param);
      //    util::Die(-1);
      //}

      double param_val = arg.GetDouble(argname);
      if ( std::fabs(param_val + 999. ) < conv::kASmallNum ) {
        param_val = Opt::Instance()->GetGlobalBestFit(param);
      }
      if(fp->Exists(param)) {
        fp->Update      (param, param_val);
        fp->SetInitValue(param, param_val);

      }
    }


    for(PhysParamId_t p : osc_params) {

      double rnd_value = util::oscillation::GetRandomizedParam(p,*fp);
      fp->Update(p, rnd_value);
      if(p == conv::kPhypMixing14    || p == conv::kPhypMixing15    || p == conv::kPhypMixing24    ||
         p == conv::kPhypMixing25    || p == conv::kPhypMixing34    || p == conv::kPhypMixing35    ||
         p == conv::kPhypDm41sq      || p == conv::kPhypDm51sq      || p == conv::kPhypDelta12     ||
         p == conv::kPhypDelta13     || p == conv::kPhypDelta15     || p == conv::kPhypDelta34     ||
         p == conv::kPhypDelta35     || p == conv::kPhypNSIepsMuTau || p == conv::kPhypMixingEE    ||
         p == conv::kPhypMixingMuMu  || p == conv::kPhypMixingMuE) {
        LOG("valor", pWARN) << "Randomising sterile params may not work as intended";
      }
    }

    fp->FromPhysParamList(model->Tweaks(),"current-silent");

    // Check that the parameter values are physical. Rethrow if they aren't
    if (!util::AllWithinLimits(*fp) && rethrow_unphysical) {
      itoy--;
      n_rethrew++;
      continue;
    }

    if(itoy % 100 == 0) { cout << "** Toy-MC: #" << itoy << endl; }
    //if(itoy % 1 == 0) { cout << "** Toy-MC: #" << itoy << endl; }

    // Generate predictions
    model->Generate(*fp,"+SB");
    //std::cout << "Print tag 2" << std::endl;
    // fp->Print(); //DEBUG

    // Array to store results
    TObjArray * array = new TObjArray();
    array->SetOwner(true);

    // Write stuff out
    model->Write(array,"curr","pdfts_rf +SB");

      // Pick results wanted
      for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
          SingleSampleFitModelABC * this_model = dynamic_cast<SingleSampleFitModelABC*>(cm->GetSubModel(imodel));
          BeamConf_t b = this_model->GetBeam();
          Sample_t   s = this_model->GetSample();
          TString name = TString::Format("curr_pdfts_sb_rf_%s_%s", this_model->Name().c_str(), util::AsString(kv).c_str());
          TH1D * curr_pdf = (TH1D*) array->FindObject(name.Data());
          if (!curr_pdf) {
              cerr<<"Could not find histogram "<<name<<endl;
              exit(-1);
          }
          if(itoy < 0) {
              nom_pdf_werr.push_back(new TH1D(*curr_pdf));
              nom_pdf     .push_back(new TH1D(*curr_pdf));
          } else {
              twk_pdf[imodel].push_back(new TH1D(*curr_pdf));
          }
      } //imodel
      array->Clear("C");
      delete array;
  }//itoy

    for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
        err_m1.push_back(new TH1D(*nom_pdf[imodel]));
        err_p1.push_back(new TH1D(*nom_pdf[imodel]));
        err_m1.back()->Reset();
        err_p1.back()->Reset();
        for(int ebin = 1; ebin <= nom_pdf[imodel]->GetNbinsX(); ebin++) {
            double nom_val = nom_pdf[imodel]->GetBinContent(ebin);
            double sumsq   = 0;
            for(int itoy = 0; itoy < ntoys; itoy++) {
                double curr_val = twk_pdf[imodel][itoy]->GetBinContent(ebin);
                sumsq += TMath::Power(curr_val-nom_val,2.);
            }//itoy
            double one_sigma = TMath::Sqrt(sumsq/ntoys);
            err_p1[imodel]->SetBinContent(ebin, nom_val + one_sigma);
            err_m1[imodel]->SetBinContent(ebin, nom_val - one_sigma);
            nom_pdf_werr [imodel] -> SetBinError(ebin,one_sigma);
            if( nom_pdf[imodel]->GetBinCenter(ebin) < enu_reco_max ) {
                max [imodel] = TMath::Max(max[imodel], nom_val + one_sigma);
            }
        }//ebin
    }//imodel


    // Draw and (if necessary), write to file
    TFile * outfile = NULL;
    if (output_rootfile!="-" && output_rootfile!="none") {
        output_rootfile = util::ReplaceAll(output_rootfile,"[PARAMS]",params_to_tweak);
        outfile = new TFile(output_rootfile.c_str(),"RECREATE");
        LOG("valor",pNOTICE)<<"Writing to file: "<<output_rootfile;
    }

    for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
        SingleSampleFitModelABC * this_model = dynamic_cast<SingleSampleFitModelABC*>(cm->GetSubModel(imodel));
        BeamConf_t b = this_model->GetBeam();
        Sample_t   s = this_model->GetSample();
        //~ ic++;
        //~ c[ibeam]->cd(ic);
        //~ TH1F * hframe = c[ibeam]->GetPad(ic)->DrawFrame(0.,0.,enu_reco_max, 1.24*max[imodel]);
        //~ hframe -> GetXaxis() -> SetTitle("Reco. E_{#nu} (GeV)");
        //~ hframe -> GetYaxis() -> SetTitle("");
        //~ //hframe -> GetYaxis() -> SetTitle("Events per bin");
        //~ hframe -> Draw();
        //~ nom_pdf_werr [imodel] -> SetFillColor(kRed);
        //~ nom_pdf_werr [imodel] -> SetLineColor(kRed);
        //~ nom_pdf_werr [imodel] -> SetLineWidth(0);
        //~ nom_pdf_werr [imodel] -> SetMarkerSize(0);
        //~ nom_pdf_werr [imodel] -> Draw("E2SAME");
        //~ nom_pdf      [imodel] -> SetLineColor(kBlack);
        //~ nom_pdf      [imodel] -> SetMarkerSize(0);
        //~ nom_pdf      [imodel] -> SetLineWidth(1);
        //~ nom_pdf      [imodel] -> Draw("SAME");
        //~ TLatex * title = new TLatex(1,1.03*max[imodel],Form("%s / %s", util::AsString(b).c_str(), util::AsStringTeX(s,true).c_str()));
        //~ title->SetTextSize(0.08);
        //~ title->Draw();
        //~ TLegend * leg = new TLegend(0.50, 0.50, 0.92, 0.80);
        //~ leg->SetFillColor(0);
        //~ //leg->SetHeader(Form("%s / %s", util::AsString(b).c_str(), util::AsStringTeX(s,true).c_str()));
        //~ leg->AddEntry(nom_pdf     [imodel], "nominal",        "L");
        //~ leg->AddEntry(nom_pdf_werr[imodel], "a-priori error", "F");
        //~ leg->Draw();

        if (outfile && outfile->IsOpen()) {
            TString name = TString::Format("env_sb_%s_%s",
                                           util::AsString(b).c_str(), util::AsString(s).c_str() );
            nom_pdf[imodel]->Write((name+"_nom").Data());
            err_p1 [imodel]->Write((name+"_p1" ).Data());
            err_m1 [imodel]->Write((name+"_m1" ).Data());
        }
    }//imodel
    //~ c[ibeam]->Update();

    if (outfile) {
        outfile->Close();
        delete outfile;
    }

  std::cout <<"Finished generating " << ntoys << " toys. Rethrew " << n_rethrew <<std::endl;
  std::cout <<"Done!"<<std::endl;
}
