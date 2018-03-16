#include <iostream>
#include <fstream>
#include <cstring>
#include <limits>

#include "SuperkMC.h"
#include "SuperkObservation.h"
#include "SuperkSamples.h"
#include "valor/VArray/VALORArrayF.h"
#include "valor/Utilities/Env.h"

#include <TMath.h>
#include <TSystem.h>

using namespace valor;
using namespace valor::conv;

int min_bin_sub_param0 = 0;
int max_bin_sub_param0 = 0;
int min_bin_sub_param1 = 0;
int max_bin_sub_param1 = 0;

double compute_chisq(VALORArrayF* varray, SuperkObservation* fakedata,
      SuperkExpectation* skexp, SuperkSample_t s, int iexp)
{
   double normPOT = skexp->GetSuperkMC()->AnalyzedSamplePOT(s) / 1.0E20;
   double chisq = 0.0;
   for(int ie = 0; ie < util::NumErecoBins(s); ++ie)
   {
      for(int it = 0; it < util::NumQ2Bins(s); ++it)
      {
         // For 1 point per file the index to the parameter is always 0
         double nexp = varray->GetEntry(iexp, 0, 0, ie,it);
         nexp *= normPOT;
         double nobs = fakedata->TwoDSpectrumVArray(s)->GetEntry(ie, it);
         chisq += valor::util::ChisqPoisson(nobs, nexp);
      }
   }
   return chisq;
}

void read_header(const char* header_filename)
{
   std::ifstream toys_stream(header_filename, std::ios::binary | std::ios::in);
   assert(toys_stream);
   toys_stream.read( (char*)&min_bin_sub_param0, sizeof(int));
   toys_stream.read( (char*)&max_bin_sub_param0, sizeof(int));
   toys_stream.read( (char*)&min_bin_sub_param1, sizeof(int));
   toys_stream.read( (char*)&max_bin_sub_param1, sizeof(int));
   toys_stream.close();
}

VALORArrayF* read_sample(std::ifstream& toys_stream)
{
    int ndim = 0;
    toys_stream.read((char*)&ndim, sizeof(int));
    int dims[6];
    toys_stream.read((char*)dims, 6 * sizeof(int));
    VALORArrayF* varray = new VALORArrayF(ndim, dims);
    int blocksize = varray->GetSize(-1);
    toys_stream.read((char*)varray->GetArray(), blocksize * sizeof(float));
    
    return varray;
}

double* read_spectra(const char* main_filename, SuperkObservation* fakedata,
      SuperkExpectation* skexp, int num_toys)
{
   std::ifstream toys_stream(main_filename, std::ios::binary | std::ios::in);
   assert(toys_stream);
   VALORArrayF* varray_1re = read_sample(toys_stream);
   VALORArrayF* varray_1rmu = read_sample(toys_stream);
   VALORArrayF* varray_rhc1re = read_sample(toys_stream);
   VALORArrayF* varray_rhc1rmu = read_sample(toys_stream);
   VALORArrayF* varray_multire = read_sample(toys_stream);

   double* chisq = new double[num_toys];
   for(int i = 0; i < num_toys; ++i)
   {
      chisq[i] = compute_chisq(varray_1re, fakedata, skexp, kSKsFCFV1RingELike, i);
      chisq[i] += compute_chisq(varray_1rmu, fakedata, skexp, kSKsFCFV1RingMuLike, i);
      chisq[i] += compute_chisq(varray_rhc1re, fakedata, skexp, kSKsRHCFCFV1RingELike, i);
      chisq[i] += compute_chisq(varray_rhc1rmu, fakedata, skexp, kSKsRHCFCFV1RingMuLike, i);
      chisq[i] += compute_chisq(varray_multire, fakedata, skexp, kSKsFCFVMultiRingELike, i);
   }

   toys_stream.close();
   delete varray_1re;
   delete varray_1rmu;
   delete varray_rhc1re;
   delete varray_rhc1rmu;
   delete varray_multire;

   return chisq;
}

double** study(const char* directory, int num_toys, int num_points,
      const char* long_hierarchy, const char* short_hierarchy,
      SuperkExpectation* skexp, SuperkObservation* fakedata)
{
   double** chisq = new double*[num_points];
   for(int i = 0; i < num_points; ++i)
   {
      std::cout << "Processing point " << (i + 1) << " of " << num_points << std::endl;
/*      char full_path_header[256];
      sprintf(full_path_header, "%s/oa2018_%s_reactor_deltacp_file%d_toy%d_x_oscpar%d_oscpar1_MH%s_subrange.bin",
               directory, long_hierarchy, i, num_toys, num_points, short_hierarchy);
      read_header(full_path_header);
      std::cout << min_bin_sub_param0 << " " << max_bin_sub_param0 << " " <<
         min_bin_sub_param1 << " " << max_bin_sub_param1 << std::endl;*/

      char full_path_main[256];
      sprintf(full_path_main, "%s/oa2018_%s_reactor_deltacp_file%d_toy%d_x_oscpar%d_oscpar1_MH%s.bin",
               directory, long_hierarchy, i, num_toys, num_points, short_hierarchy);
      chisq[i] = read_spectra(full_path_main, fakedata, skexp, num_toys);
   }

   return chisq;
}

void AsimovA(PhysHypothesis& ph)
{
   ph.SetCurrValue(conv::kPhypMixing12, 0.304);
   ph.SetCurrValue(conv::kPhypMixing13, 0.0219);
   ph.SetCurrValue(conv::kPhypMixing23, 0.528);
   ph.SetCurrValue(conv::kPhypDm12sq, 7.53E-5);
   ph.SetCurrValue(conv::kPhypDm23sq, 2.509E-3);
   ph.SetCurrValue(conv::kPhypCPV, -1.601);
}

SuperkObservation* generate_fake_data(SuperkExpectation& skexp)
{
  
   // Global options
   Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
   Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString("auto");
   Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString("sinsqtheta");
   Opt::Instance()->gOscOctantChoice      = util::OctantFromString("two");
   Opt::Instance()->SetDefaultParameterValues();
  

   SuperkObservation* skobs_fakedata = new SuperkObservation();
   PhysHypothesis ph("3f_osc");
   AsimovA(ph);
    
   // No need to throw systematics
   skobs_fakedata->CreateMockData(skexp, ph, false);

   return skobs_fakedata;
}

void marginalise(double** chisq, const int N)
{
   const double SCALE = 80;
   int num_toys = 40000 / N;
   double** average_chisq = new double*[N];
   double** dchisq = new double*[N];
   for(int i = 0; i < N; ++i)
   {
      average_chisq[i] = new double[101];
      dchisq[i] = new double[101];
   }
   for(int i = 0; i < N; ++i)
   {
      int first = i * num_toys;
      int last = (i + 1) * num_toys;
/*      std::cout << "i, first, last, #toys " << i << " " << first <<
         " " << last << " " << num_toys << std::endl;*/
      for(int j = 0; j < 101; ++j)
      {
         // Find the minima and maximum chisq values in the range
         double min = std::numeric_limits<double>::max();
         double max = -std::numeric_limits<double>::max();
         for(int x = first; x < last; ++x)
         {
            if(chisq[j][x] < min)
            {
               min = chisq[j][x];
            }
            if(chisq[j][x] > max)
            {
               max = chisq[j][x];
            }
         }
         // Yes, two factors of 0.5 is correct
         double diff_chisq = 0.5 * 0.5 * (max - min);
         double scale = TMath::Min(diff_chisq, SCALE);
         double likelihood = 0.0;
         for(int x = first; x < last; ++x)
         {
            likelihood += exp(-chisq[j][x] * 0.5 + scale) / num_toys;
         }
         average_chisq[i][j] = 2 * scale - 2 * log(likelihood); 
/*         std::cout << "j, min, max, diff, scale, likelihood, mean " << j <<
            " " << min << " " << max << " " << diff_chisq << " " << scale <<
            " " << likelihood << " " << average_chisq[i][j] << std::endl;*/
      }
      double min = std::numeric_limits<double>::max();
      for(int j = 0; j < 101; ++j)
      {
         if(average_chisq[i][j] < min)
         {
            min = average_chisq[i][j];
         }
      }
      for(int j = 0; j < 101; ++j)
      {
         dchisq[i][j] = average_chisq[i][j] - min;
      }
   }

   for(int i = 0; i < N; ++i)
   {
      for(int j = 0; j < 101; ++j)
      {
         std::cout << dchisq[i][j] << ",";
      }
      std::cout << std::endl;
   }
   for(int i = 0; i < N; ++i)
   {
      delete[] average_chisq[i];
      delete[] dchisq[i];
   }
   delete[] average_chisq;
   delete[] dchisq;
}

void marginalisation_study(const char* directory, const char* long_hierarchy, const char* short_hierarchy)
{
   Env::Instance()->AutoSet("t2xk");

   // Analysis configuration
   std::string config(Form("%s/config/analysis_configuration_oa2018_run1-9c_v1_postfit.xml",
            gSystem->Getenv("T2KOSCANA")));

   // Get superk spectrum predictor & configure
   SuperkExpectation* skexp = new SuperkExpectation;
   skexp->Read(config);

   SuperkObservation* fakedata = generate_fake_data(*skexp);
   for(int is = 0; is < kns; ++is)
   {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined)
      {
         continue;
      }
      std::cout << util::AsStringCompact(s) << " " <<
         fakedata->TwoDSpectrumVArray(s)->Integral() << endl;
   }

   double** chisq = study(directory, 40000, 101, long_hierarchy, short_hierarchy, skexp, fakedata);
   std::cout << std::endl << "400" << std::endl << std::endl;
   marginalise(chisq, 400);
   std::cout << std::endl << "200" << std::endl << std::endl;
   marginalise(chisq, 200);
   std::cout << std::endl << "100" << std::endl << std::endl;
   marginalise(chisq, 100);

/*   for(int i = 0; i < 40000; ++i)
   {
      for(int j = 0; j < 101; ++j)
      {
         std::cout << chisq[j][i] << ",";
      }
      std::cout << std::endl;
   }*/
   for(int j = 0; j < 101; ++j)
   {
      delete[] chisq[j];
   }
   delete[] chisq;

   delete skexp;
   delete fakedata;
}

