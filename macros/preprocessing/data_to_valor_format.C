#include <cmath>
#include <vector>
#include <utility>

#include "Build.h"
#include "AnalysisOptions.h"

const double M_n = 939.565346;
const double M_p = 938.272013;
const double M_m = 105.65836;
const double Eb = 27.0;

bool elike_cut(const double fqevis, const double erec)
{
   return fqevis > 100 && erec < 1250;
}

bool mulike_cut(const double fqmomm)
{
   return fqmomm > 200;
}

double erec_ccqe(double p_lep, double M_lep, double costheta, double numode)
{
   double M_1 = M_p;
   double M_2 = M_n;
   if(numode)
   {
      M_1 = M_p;
      M_2 = M_n;
   }
   else
   {
      M_1 = M_n;
      M_2 = M_p;
   }
   double E_lep = sqrt(p_lep * p_lep + M_lep * M_lep);
   double E_rec = (M_1 * M_1 - (M_2 - Eb) * (M_2 - Eb) - M_lep * M_lep + 2 * (M_2 - Eb) * E_lep) / (2 * (M_2 - Eb - E_lep + p_lep * costheta));
   
   return E_rec;
}

double erec_cc1pi(double p_lep, double M_lep, double costheta)
{
   double E_lep = sqrt(p_lep * p_lep + M_lep * M_lep);
   double E_rec = M_p * E_lep - 0.5 * M_lep * M_lep;
   E_rec += 0.5 * (M_d * M_d - M_p * M_p);
   E_rec /= (M_p - E_lep + p_lep * costheta);

   return E_rec;
}

double coulomb_correction(double fqmomm, bool numode)
{
    // Apply Coulomb correction to muon momentum (FHC: +4.3 MeV; RHC: -3.3 MeV - inverting the MC simulation)
    if(numode)
    {
        return fqmomm + 4.3;
    }
    else
    {
        return fqmomm - 3.3;
    }    
}

void data_to_valor_format(const char * tag = "SKdata_Run18",
      const char * input_file = "$T2KOSCANA/inputs/data_run1to8/t2ksk.14c.fqv4r0.run1-8.mr74.data.minituple.root",
      char *treename = "minituple")
{
//   gROOT->ProcessLine(".L $T2KOSCANA/macros/preprocessing/data_to_valor_format.h");
   using namespace std;
   using namespace valor;
   using namespace valor::conv;
  
   const int good_horn_current_fhc_a = 200;
   const int good_horn_current_fhc_b = 250;
   const int good_horn_current_rhc = -250;
  
   const int nereco_1rmu = util::NumErecoBins(conv::kSKsFCFV1RingMuLike);
   const int nq2_1rmu = util::NumQ2Bins(conv::kSKsFCFV1RingMuLike);
   TH2D *mock_data_1Rmu_0 = new TH2D("mock_data_1Rmu_0","mock_data_1Rmu_0",nq2_1rmu,util::Q2BinEdges(kSKsFCFV1RingMuLike),nereco_1rmu,util::ErecoBinEdges(kSKsFCFV1RingMuLike));
  
   const int nereco_rhc1rmu = util::NumErecoBins(conv::kSKsRHCFCFV1RingMuLike);
   const int nq2_rhc1rmu = util::NumQ2Bins(conv::kSKsRHCFCFV1RingMuLike);
   TH2D *mock_data_RHC1Rmu_0 = new TH2D("mock_data_RHC1Rmu_0","mock_data_RHC1Rmu_0",nq2_rhc1rmu,util::Q2BinEdges(kSKsRHCFCFV1RingMuLike),nereco_rhc1rmu,util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike));
  
   const int nereco_1re = util::NumErecoBins(conv::kSKsFCFV1RingELike);
   const int nq2_1re = util::NumQ2Bins(conv::kSKsFCFV1RingELike);
   TH2D *mock_data_1Re_0 = new TH2D("mock_data_1Re_0","mock_data_1Re_0",nq2_1re,util::Q2BinEdges(kSKsFCFV1RingELike),nereco_1re,util::ErecoBinEdges(kSKsFCFV1RingELike));
  
   const int nereco_rhc1re = util::NumErecoBins(conv::kSKsRHCFCFV1RingELike);
   const int nq2_rhc1re = util::NumQ2Bins(conv::kSKsRHCFCFV1RingELike);
   TH2D *mock_data_RHC1Re_0 = new TH2D("mock_data_RHC1Re_0","mock_data_RHC1Re_0",nq2_rhc1re,util::Q2BinEdges(kSKsRHCFCFV1RingELike),nereco_rhc1re,util::ErecoBinEdges(kSKsRHCFCFV1RingELike));

   const int nereco_multire = util::NumErecoBins(conv::kSKsFCFVMultiRingELike);
   const int nq2_multire = util::NumQ2Bins(conv::kSKsFCFVMultiRingELike);
   TH2D *mock_data_MultiRe_0 = new TH2D("mock_data_MultiRe_0","mock_data_MultiRe_0",nq2_multire,util::Q2BinEdges(kSKsFCFVMultiRingELike),nereco_multire,util::ErecoBinEdges(kSKsFCFVMultiRingELike));

   TFile * fin = new TFile(input_file,"OPEN");
   assert(fin);
   TTree * tree = (TTree*) fin->Get(treename);
 
   float fqmomm;
   float fqmome;
   float fqnse;
   double erec;
   float fqevis;
   int iclass;
   float dir[2][7][3];
   double horn_current;
  
   tree->SetBranchAddress("erec", &erec);
   tree->SetBranchAddress("iclass", &iclass);
//  tree->SetBranchAddress("fqnse", &fqnse);
   tree->SetBranchAddress("fqmomm", &fqmomm);
   tree->SetBranchAddress("fqmome", &fqmome);
   tree->SetBranchAddress("horn_current", &horn_current);
   tree->SetBranchAddress("fq1rdir", &dir);
   tree->SetBranchAddress("fqevis", &fqevis);
   int nentries = tree->GetEntriesFast();
  
  
   // Dataset info
   int iexperiment = 0;
  
   int n_fhc_mu=0,n_rhc_mu=0,n_fhc_e=0,n_rhc_e=0,n_multie=0;
  
   SuperkMC skmc;

   std::vector<std::pair<float,float>> fhc1RePts;
   std::vector<std::pair<float,float>> rhc1RePts;
   std::vector<std::pair<float,float>> multiRePts;
  
   for(int i=0;i<nentries;i++)
   {
      tree->GetEntry(i);
    
      if(erec < 0.)
      {
         continue;
      }

      double dir_nu[3] = { 0.669764, -0.742179, 0.024223 }; // T2K beam direction
      bool is_electron = iclass == 11 || iclass == 14;
      double costheta = 0.0;

      if(is_electron)
      {
         costheta = dir_nu[0] * dir[0][1][0] + dir_nu[1] * dir[0][1][1] + dir_nu[2] * dir[0][1][2];
      }
      else
      {
         costheta = dir_nu[0] * dir[0][2][0] + dir_nu[1] * dir[0][2][1] + dir_nu[2] * dir[0][2][2];
      }

      double LepTheta = TMath::ACos(costheta); // radiant
      double theta = LepTheta / TMath::Pi() * 180.;
    
      if(theta < 0.)
      {
         continue;
      }
    
      if(iclass == 11)
      {  // e-like
         std::cout << (horn_current > 0 ? "e-:: " : "e+:: ") << "evis = " << fqevis << " erec = " << erec;
         if(!elike_cut(fqevis, erec))
         {
            std::cout << " : Not selected" << std::endl;
            continue;
         }
         std::cout << " : Selected" << std::endl;
         erec /= 1000.0;
         if(horn_current > 0)
         {  // fhc
            mock_data_1Re_0->Fill(theta, erec);
            fhc1RePts.push_back(std::pair<float,float>(erec, theta));
            n_fhc_e++;
            cout << "Erec = " << erec << " - Theta = " << theta << endl;
         }
         else if(horn_current < 0)
         {  // rhc
            mock_data_RHC1Re_0->Fill(theta, erec);
            rhc1RePts.push_back(std::pair<float,float>(erec, theta));
            n_rhc_e++;
         }
      }
      else if(iclass == 14)
      {  //MultiRe
         std::cout << (horn_current > 0 ? "e- cc1pi:: " : "e+ cc1pi:: ") << "evis = " << fqevis << " erec = " << erec;
         if(!elike_cut(fqevis, erec))
         {
            std::cout << " : Not selected" << std::endl;
            continue;
         }
         std::cout << " : Selected" << std::endl;
         erec /= 1000.0;
         if(horn_current > 0)
         {  // fhc
            mock_data_MultiRe_0->Fill(theta,erec);
            multiRePts.push_back(std::pair<float,float>(erec, theta));
            n_multie++;
         }
      }
      else if(iclass == 13)
      {  // mu-like
         std::cout << (horn_current > 0 ? "mu-:: " : "mu+:: ") << "p = " << fqmomm << " erec= " << erec;
         fqmomm = coulomb_correction(fqmomm, horn_current > 0);
         erec = erec_ccqe(fqmomm, M_m, costheta, horn_current > 0);
         std::cout << " => p = " << fqmomm << " erec = " << erec;
         if(!mulike_cut(fqmomm))
         {
            std::cout << " : Not selected" << std::endl;
            continue;
         }
         std::cout << " : Selected" << std::endl;

         erec /= 1000.0;
         if(horn_current > 0)
         {  // fhc
            mock_data_1Rmu_0->Fill(theta, erec);
            n_fhc_mu++;
         }
         else if(horn_current < 0)
         {  // rhc
            mock_data_RHC1Rmu_0->Fill(theta, erec);
            n_rhc_mu++;
         }
      }
   }
  
   TGraph gfhc1Re(fhc1RePts.size());
   gfhc1Re.SetTitle("g1Re");
   gfhc1Re.SetName("g1Re");
   for(int i = 0; i < fhc1RePts.size(); ++i)
   {
      std::pair<float, float> val = fhc1RePts[i];
      gfhc1Re.SetPoint(i, val.first, val.second);
   }

   TGraph grhc1Re(rhc1RePts.size());
   grhc1Re.SetTitle("gRHC1Re");
   grhc1Re.SetName("gRHC1Re");
   for(int i = 0; i < rhc1RePts.size(); ++i)
   {
      std::pair<float, float> val = rhc1RePts[i];
      grhc1Re.SetPoint(i, val.first, val.second);
   }
 
   TGraph gmultiRe(multiRePts.size());
   gmultiRe.SetTitle("gMultiRe");
   gmultiRe.SetName("gMultiRe");
   for(int i = 0; i < multiRePts.size(); ++i)
   {
      std::pair<float, float> val = multiRePts[i];
      gmultiRe.SetPoint(i, val.first, val.second);
   }
  
   // Print # of events
   cout <<"\n\nEvent rates:"<<endl;
   cout <<" FHC 1Rmu: "<<n_fhc_mu<<endl;
   cout <<" FHC 1Re:  "<<n_fhc_e <<endl;
   cout <<" RHC 1Rmu: "<<n_rhc_mu<<endl;
   cout <<" RHC 1Re:  "<<n_rhc_e <<endl;
   cout <<" MultiRe:  "<<n_multie <<endl;
  
   fin->Close();
  
   TCanvas *c_1Rmu = new TCanvas("c_1Rmu","c_1Rmu");
   mock_data_1Rmu_0->Draw("colz");
  
   TCanvas *c_RHC1Rmu = new TCanvas("c_RHC1Rmu","c_RHC1Rmu");
   mock_data_RHC1Rmu_0->Draw("colz");
  
   TCanvas *c_1Re = new TCanvas("c_1Re","c_1Re");
   mock_data_1Re_0->Draw("colz");
  
   TCanvas *c_RHC1Re = new TCanvas("c_RHC1Re","c_RHC1Re");
   mock_data_RHC1Re_0->Draw("colz");
  
   TCanvas *c_MultiRe = new TCanvas("c_MultiRe","c_MultiRe");
   mock_data_MultiRe_0->Draw("colz");
  
   // Store histograms
   char name[200];
   sprintf(name,"%s.root",tag);
   TFile *fout = new TFile(name,"RECREATE");
   mock_data_1Rmu_0->Write();
   mock_data_RHC1Rmu_0->Write();
   mock_data_1Re_0->Write();
   mock_data_RHC1Re_0->Write();
   mock_data_MultiRe_0->Write();
   gfhc1Re.Write();
   grhc1Re.Write();
   gmultiRe.Write();
   fout->Close();
  
   return;
}

