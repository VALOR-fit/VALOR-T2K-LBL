//____________________________________________________________________________
/*!

\app       valor_dune-run_fit_marg

\brief     Fit a set of marginalisation toys.

\options     --tag                 : Tag for input marginalisation toy file [def: "margtoy"]
             --run                 : Run # [def: 0]
             --nmargtoy            : Number of marginalisation toys [def: 1]
             --physics-hypothesis  : Physics hypothesis (see UtlParam.cxx)
                                     [def: "nuosc3f"]
             --fit-param-opt       : Options used when setting up the FitParamList
                                     [def: "NONE"]
             --param-values-option : Decide how to set the values of physics 
                                     params used for toy-MC generation 
                                    [def:"randomize-all;prefit-errors"]
             --config              : Analysis configuration file [required]
             --mesg-thresholds     : Message thresholds [def: "verbosity.xml"]
            [--with-stat-fluct |   : Include statistical fluctuation in toy 
                                     MC generation? [def: yes]
             --without-stat-fluct] 
            [--store-spectra   |   : Store toy-MC/best-fit spectra? [def: no]
             --discard-spectra   ] 
             --n-sigma-limit       : The limit to apply to systematic parameter
                                     throws
             --use-simple-binning  : Update simple binning templates too [def: no-simple-binning]

\author    Andy Chappell <A.Chappell \at warwick.ac.uk>
           University of Warwick

\created   January 14, 2019

\copyright (c) 2010-2019, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________

// root
#include "TFile.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TStopwatch.h"

#include "cmath"
#include "fstream"      // std::ofstream

// fmwk
#include "valor/CmdLine/CmdLnArgList.h"

#include "valor/Conventions/KineSpace.h"
#include "valor/Conventions/Opt.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/FitCore/FitInputs.h"
#include "valor/FitCore/FitterABC.h"
#include "valor/FitCore/UtlFit.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/Marginalisation.h"
#include "valor/FitCore/UtlOscillation.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/ToyStream.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"

using namespace valor;
using namespace valor::conv;
//using namespace valor::marginalisation;

void Fit(CmdLnArgList& arg);
void Init(CmdLnArgList& arg);

int main(int argc, char** argv)
{

   CmdLnArgList arg;
   arg.Read(argc,argv,
      "$VALOR_EXPT/config/apps/run_fit_marg.arg"
      ":$VALOR_EXPT/config/apps/seeds.arg"
      ":$VALOR_EXPT/config/apps/oscparams.arg");
   LOG("valor_dune", pNOTICE) << arg;
   Fit(arg);
   LOG("valor_dune", pNOTICE) << "Finished";
  
   return 0;
}

void Fit(CmdLnArgList& arg)
{
   Init(arg);

   // Get arguments
   string hypname         = arg.GetString("physics-hypothesis");
   string fplopt           = arg.GetString("fit-param-opt");
   string param_val_opts   = arg.GetString("param-values-option");
   string fitter_name      = arg.GetString("fitter");
   int Ntoys_fake_data     = arg.GetInt("ntoys-fakedata");
   int start_toy_fake_data = arg.GetInt("start-toy-fakedata");
   bool with_stat_fluct    = arg.GetBool("with-stat-fluct");
   bool rethrow_unphysical = arg.GetBool("rethrow-unphysical");
   bool use_simple_binning = arg.GetBool("use-simple-binning");

   string fp1_name         = arg.GetString("fit-param1");
   const double fp1_val    = arg.GetDouble("fit-param1-value");
   const int fp1_bin       = arg.GetInt("fit-param1-bin");
   const double fp1_minval = arg.GetDouble("fit-param1-min-value");
   const double fp1_maxval = arg.GetDouble("fit-param1-max-value");
   const int fp1_nbins     = arg.GetInt("fit-param1-nbins");
   const int fp1_start_bin = arg.GetInt("fit-param1-start-bin");
   string penaltyparam  = arg.GetString("use-penalty-term");

   bool fit_is_2D = false;
   string fp2_name;
   double fp2_val = 0.0;
   int fp2_bin = 0;
   double fp2_minval = 0.0;
   double fp2_maxval = 0.0;
   int fp2_nbins = 0;
   int fp2_start_bin = 0;

   LOG("valor_dune", pNOTICE) << "Fit parameter " << fp1_name << "(" <<
      fp1_bin << ") = " << fp1_val << " [" << fp1_nbins << ": " << fp1_minval <<
      ", " << fp1_maxval << "]";

   if(arg.Exists("fit-param2") and arg.GetString("fit-param2") != "-")
   {
      fit_is_2D     = true;
      fp2_name      = arg.GetString("fit-param2");
      fp2_val       = arg.GetDouble("fit-param2-value");
      fp2_bin       = arg.GetInt("fit-param2-bin");
      fp2_minval    = arg.GetDouble("fit-param2-min-value");
      fp2_maxval    = arg.GetDouble("fit-param2-max-value");
      fp2_nbins     = arg.GetInt("fit-param2-nbins");
      fp2_start_bin = arg.GetInt("fit-param2-start-bin");
      LOG("valor_dune", pNOTICE) << "Fit parameter " << fp2_name << "(" <<
         fp2_bin << ") = " << fp2_val << " [" << fp2_nbins << ": " <<
         fp2_minval << ", " << fp2_maxval << "]";
   }

   bool read_data_from_binary = (arg.GetString("data-file") != "-" && !arg.GetString("data-file").empty());
   if(read_data_from_binary && arg.GetInt("gen-fakedata")) {
      LOG("valor_dune", pFATAL) << "Can't both read a fake data file and generate fake data "
                                << "(don't use --data-file and --gen-fakedata together)";
      util::Die(-1);
   }

    bool fusepenalty = false;

    if(arg.Exists("use-penalty-term") and arg.GetString("use-penalty-term") != "-" ){

        if (penaltyparam == fp1_name || penaltyparam == fp2_name ){fusepenalty = true;}
        else{
            LOG("valor_dune", pFATAL) << "Penalty term specified for a parameter that is different from fitting one(s)" <<
                                      " penalty set for " << penaltyparam;
            util::Die(-1);
        }
    }


   CompoundFitElements* elements = new CompoundFitElements();
   elements->Read(util::CompleteConfigPath(arg.GetString("config")));
   FitModelABC* model = elements->FitModel();
   FitDataABC* data = elements->FitData();

   const char* read_opt_data = "";
   if (use_simple_binning) {
      read_opt_data = "+SB";
   }

   util::ToyStream* toy_stream = nullptr;
   if(arg.GetInt("gen-fakedata")) {
      // Generate a fake data set to fit to

      if(Ntoys_fake_data > 1) {
         LOG("valor_dune", pFATAL) << "Currently only fits with 1 fakedata toy supported when generating on the fly. "
                                      "It should be relatively easy to add support for > 1 toy if you need it!";
         util::Die(-1);
      }

      FitParamList* fpl = util::GetFitParameterList(hypname.c_str(), fplopt.c_str());

      // Set the priors
      util::oscillation::SetOscillationParameterPriors(arg, *fpl);
      // Fix or set marginalisation flag for oscillation parameters
      util::oscillation::SetOscillationParameters(arg, *fpl);

      // Initialize non delta M squared oscillation parameters to zero
      // util::oscillation::SetUnoscillated(*fpl);
      LOG("valor_dune", pNOTICE) << *fpl;

      // Call Generate - forces a call to Book as well as providing the fake data
      data->Generate(model, *fpl, with_stat_fluct);
   }
   else if(read_data_from_binary) {
      // First generate a dummy fake data toy to be overwritten (so that the data arrays can be initialised)
      FitParamList* fpl_void = util::GetFitParameterList("void", "");
      data->Generate(model, *fpl_void, with_stat_fluct);

      // Fit to data read from an input binary toy file
      std::string data_file_name = arg.GetString("data-file");
      LOG("valor_dune", pNOTICE) << "Will read fake/real data from file [" << data_file_name << "]";
      toy_stream = data->InitToyStream(data_file_name, "read", read_opt_data);
   }
   else {
      // Fit to a data file specified in the analysis configuration XML file.
      // This has already been read in CompoundFitElements::AdoptFitData() (when reading the config file),
      // so there's no need to do anything here.
      if(!data->DataExists(read_opt_data)) {
         LOG("valor_dune", pFATAL) << "Tried to fit to fake/real data specified in the analysis configuration file, "
                                   << "but the data arrays for one or more samples were empty";
         util::Die(-1);
      }
      LOG("valor_dune", pNOTICE) << "Will fit to fake/real data specified in the analysis configuration file";
   }

   // Get the requested batch number and number of marginalisation toys for the
   // batch
   int Nbatches = arg.GetInt("nbatches");
   int firstBatch = arg.GetInt("first-batch");
   int Ntoys = arg.GetInt("nmargtoy");
   util::oscillation::PenaltyTerm penaltyset;

   TStopwatch clock;
   clock.Start();
   Marginalisation marginalisation(*elements, arg.GetString("margtoy-tag"), use_simple_binning);
   if(!fit_is_2D)
   {  // Initialise for 1D fit
      marginalisation.InitFit(fp1_nbins, fp1_start_bin, fp1_minval, fp1_maxval);
   }
   else
   {  // Initalise for 2D fit
      marginalisation.InitFit(fp1_nbins, fp1_start_bin, fp1_minval, fp1_maxval,
                              fp2_nbins, fp2_start_bin, fp2_minval, fp2_maxval);
   }
   FitResult* result = nullptr;

   std::vector<double>* chi2_best_fit_toy = new std::vector<double>(Ntoys_fake_data);

   // Skip forward in the input fakedata binary file to the position of the first toy to read
   if(read_data_from_binary && arg.GetInt("gen-fakedata") == 0 && start_toy_fake_data != 0)
   {
     toy_stream->SetPositionToToy(start_toy_fake_data);
   }

   for(int ifake = start_toy_fake_data; ifake < start_toy_fake_data + Ntoys_fake_data; ifake++)
   {
      // Read a fakedata toy, if requested
      if(read_data_from_binary && arg.GetInt("gen-fakedata") == 0)
      {
        data->Read(toy_stream);

#ifdef VALOR_EXPT_LOW_LEVEL_MESG_ENABLED
        LOG("valor_dune", pNOTICE) << "Fake data toy #" << ifake << " event rates:\n" << data->GetPDFIntegralAsString();
#endif
      }

      // Perform a fit and populate the entire chi2 and likelihood surfaces
      // This version of the fit works for both 1D and 2D
      if(fusepenalty)
      {
        penaltyset = util::oscillation::SetPenaltyTermParameters(arg, penaltyparam);
      }
      else { penaltyset = {0,0,"",""};}
      result = marginalisation.Fit(data, Ntoys, Nbatches, penaltyset, firstBatch);

      // Perform an optimised fit starting at the specified value
   /*   if(!fit_is_2D)
      {
         result = marginalisation.Fit(fp1_val, data, Ntoys, Nbatches, firstBatch);
      }
      else
      {
         result = marginalisation.Fit(fp1_val, fp2_val, data, Ntoys, Nbatches,
               firstBatch);
      }*/

      if(ifake % 100 == 0) {
         if(!fit_is_2D)
         {
            LOG("valor_dune", pNOTICE) << "For fakedata toy #" << ifake <<
               ", Best fit: " << fp1_name << " = " <<
               result->GetParam1BestFit() << " Chi^2 = " << result->GetChi2BestFit() <<
               " " << " L = " << result->GetLikelihoodBestFit();
         }
         else
         {
            LOG("valor_dune", pNOTICE) << "For fakedata toy #" << ifake <<
               ", Best fit param 1: " << fp1_name << " = " << result->GetParam1BestFit() <<
               ", Best fit param 2: " << fp2_name << " = " << result->GetParam2BestFit() <<
               result->GetChi2BestFit() << " " << " L = "  << result->GetLikelihoodBestFit();
         }
      }

      chi2_best_fit_toy->at(ifake - start_toy_fake_data) = result->GetChi2BestFit();
   }
   delete toy_stream;

   std::ostringstream oss;
   oss << arg.GetString("tag");
   oss << "_" << start_toy_fake_data << "_" << start_toy_fake_data + Ntoys_fake_data;
   oss << "_result.root";
   string filename = oss.str();
   TFile* fitFile = new TFile(filename.c_str(), "RECREATE");

   // Store surfaces
   LOG("valor_dune", pNOTICE) << "Only storing the chi2 and likelihood surfaces of the last fakedata toy (to avoid bloat)";
   if(!fit_is_2D)
   {
      TH1D* chi2Surface = result->GetChi2Surface1D();
      TH1D* likelihoodSurface = result->GetLikelihoodSurface1D();
      chi2Surface->Write("chi2");
      likelihoodSurface->Write("likelihood");
   }
   else
   {
      TH2D* chi2Surface = result->GetChi2Surface2D();
      TH2D* likelihoodSurface = result->GetLikelihoodSurface2D();
      chi2Surface->Write("chi2");
      likelihoodSurface->Write("likelihood");
   }

   // Store vector of marginal chi2 values
   fitFile->WriteObject(chi2_best_fit_toy, "chi2_best_fit_toy");

   fitFile->Close();
   delete fitFile;

/*   marginalisation.ComputeStatistics(fp1_bin, data, Ntoys, Nbatches, firstBatch);
   double likelihood = marginalisation.GetLikelihood();
   double chi2 = marginalisation.GetChi2();
   LOG("valor_dune", pNOTICE) << "chi2(marg): " << chi2 << " L(marg): " <<
      likelihood;*/
   clock.Stop();

   LOG("valor_dune", pNOTICE) << "Read: " << Ntoys * Ntoys_fake_data <<
      " spectra. Real time: " << clock.RealTime() << "   CPU time: " <<
      clock.CpuTime();
}

void Init(CmdLnArgList& arg)
{
   Messenger::Instance()->Configure(util::CompleteConfigPath(arg.GetString("mesg-thresholds")));

   Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);
   Opt::Instance()->gOscMixingParamChoice = ( (arg.GetBool("sinsq2theta"))  ? conv::kNuMixSinsq2Theta : conv::kNuMixSinsqTheta);
   Opt::Instance()->gOscOctantChoice      = ( (arg.GetBool("second-octant"))? conv::kNuOctantTwo      : conv::kNuOctantOne);

  conv::NuOscDmsqScheme_t dmsq_scheme = util::DmsqSchemeFromString(arg.GetString("dmsq-scheme"));
  if(dmsq_scheme == kNuDmsqUndefined) {
    LOG("valor_dune", pFATAL) << "Could not parse the input Dmsq scheme: [" << arg.GetString("dmsq-scheme") << "]";
    util::Die(-1);
  }
  Opt::Instance()->gOscDmsqScheme = dmsq_scheme;

   Rnd::Instance()->SetSeeds(arg);
}

