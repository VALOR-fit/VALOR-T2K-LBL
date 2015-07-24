
//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"


//
// It takes expected histogram from test_sk_ereco_spectrum_nominal_mainmode_nuebarapp.C output
//

using namespace std;
using namespace valor;
using namespace valor::conv;

void plot_ereco_on_top_of_prediction
(
 const char * sk_sample            = "RHC1Re",                         /* options: "1Rmu", "1Re", "MultiRmu" */
 const char * tag    = "fitmarg_betabar_DEV",    /* tag for output file */
 
 const char * expfile = "sk_ereco_spectrum_RHC1Re_prova-sinsqtheta[12_0.306,13_0.025,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007500,23_0.002509]eV2,betabar_1.000-RHC1Re.root",
 const char * xmlfile = "../../inputs/data_run1to6/skdata_run1to6.xml",
 )
{



   // //example of macro illustrating how to superimpose two histograms
   // //with different scales in the "same" pad.
   // // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   // //Author: Rene Brun

   // TCanvas *c1 = new TCanvas("c1","hists with different scales",600,400);

   // //create/fill draw h1
   // gStyle->SetOptStat(kFALSE);
   // TH1F *h1 = new TH1F("h1","my histogram",100,-3,3);
   // Int_t i;
   // for (i=0;i<10000;i++) h1->Fill(gRandom->Gaus(0,1));
   // h1->Draw();
   // c1->Update();

   // //create hint1 filled with the bins integral of h1
   // TH1F *hint1 = new TH1F("hint1","h1 bins integral",100,-3,3);
   // Float_t sum = 0;
   // for (i=1;i<=100;i++) {
   //    sum += h1->GetBinContent(i);
   //    hint1->SetBinContent(i,sum);
   // }

   // //scale hint1 to the pad coordinates
   // Float_t rightmax = 1.1*hint1->GetMaximum();
   // Float_t scale = gPad->GetUymax()/rightmax;
   // hint1->SetLineColor(kRed);
   // hint1->Scale(scale);
   // hint1->Draw("same");

   // //draw an axis on the right side
   // TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
   //       gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   // axis->SetLineColor(kRed);
   // axis->SetLabelColor(kRed);
   // axis->Draw();

   // return;










  // Get Expectation
  
  TFile *fexp = new TFile(expfile,"READ");
  TH1D *hErecoExp = (TH1D*) fexp->Get("sk_RHC1Retotal");

  TH1D *hErecoMode_NumuBarToNueBar = (TH1D*) fexp->Get("hErecoMode_NumuBarToNueBar");
  TH1D *hErecoMode_NC = (TH1D*) fexp->Get("hErecoMode_NC");
  TH1D *hErecoMode_BkgOther = (TH1D*) fexp->Get("hErecoMode_BkgOther");
  TH1D *hErecoMode_NumuToNue = (TH1D*) fexp->Get("hErecoMode_NumuToNue");
 
  TH1D *hEreco_AllBkg = (TH1D*) hErecoMode_NC->Clone("hEreco_AllBkg");
  hEreco_AllBkg->Add(hErecoMode_BkgOther);
  hEreco_AllBkg->Add(hErecoMode_NumuToNue);
 
  TH1D *hErecoMode_NumuBarToNueBar_norm = (TH1D*) fexp->Get("hErecoMode_NumuBarToNueBar_norm");
  TH1D *hErecoMode_NC_norm = (TH1D*) fexp->Get("hErecoMode_NC_norm");
  TH1D *hErecoMode_BkgOther_norm = (TH1D*) fexp->Get("hErecoMode_BkgOther_norm");
  TH1D *hErecoMode_NumuToNue_norm = (TH1D*) fexp->Get("hErecoMode_NumuToNue_norm");
  
  TH1D *hEreco_AllBkg_norm = (TH1D*) hEreco_AllBkg->Clone("hEreco_AllBkg_norm");
  hEreco_AllBkg_norm->Scale(1./hEreco_AllBkg_norm->Integral());

  cout << "Tot = " << hErecoExp->Integral() << endl;
  cout << "Tot 2 = " 
       << hErecoMode_NumuBarToNueBar->Integral() + 
    hErecoMode_NC->Integral() + 
    hErecoMode_BkgOther->Integral() + 
    hErecoMode_NumuToNue->Integral()
       << endl; 
  cout << "hErecoMode_NumuBarToNueBar = " << hErecoMode_NumuBarToNueBar->Integral() << endl; 
  cout << "hErecoMode_NC = " << hErecoMode_NC->Integral() << endl; 
  cout << "hErecoMode_BkgOther = " << hErecoMode_BkgOther->Integral() << endl; 
  cout << "hErecoMode_NumuToNue = " << hErecoMode_NumuToNue->Integral() << endl;
  // return;


  // Get Data

  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));
  
  FileStat_t dummy;
  if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
    {
      xmlfile2 = xmlfile;
      if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
	cerr<<"Data input file not found: "<<xmlfile<<endl;
	exit(-1);
      }
    }
 
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  bool using_xml_file = true;
  if(using_xml_file){
    skobs_fakedata->Read(xmlfile2.Data(),0,using_xml_file);
  }
  else{
    skobs_fakedata->Read(rootfile2.Data(),0,using_xml_file);
  }
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  TH1D *hErecoData;   
  hErecoData = (TH1D*)skobs_fakedata->ErecoSpectrum(s)->Clone("hErecoData");
  
  // Draw 
  
  TCanvas *Ereco_ExpVsData = new TCanvas("Ereco_ExpVsData","Ereco_ExpVsData");
  hErecoData->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoData->GetYaxis()->SetTitle("Number of events");  
  hErecoExp->SetLineColor(kBlue);
  hErecoData->DrawCopy("pe");
  hErecoExp->DrawCopy("histo same");  
  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hErecoExp,"Expected","L");
  leg->AddEntry(hErecoData,"Data","P");
  leg->Draw();
  
  // TCanvas *cEreco_ExpVsData_Component = new TCanvas("cEreco_ExpVsData_Component","cEreco_ExpVsData_Component");
  // hErecoMode_NC->SetFillColor(0);
  // hErecoMode_NC->SetLineColor(kRed);
  // hErecoMode_BkgOther->SetFillColor(0);
  // hErecoMode_BkgOther->SetLineColor(kBlue);
  // hErecoMode_NumuToNue->SetFillColor(0);
  // hErecoMode_NumuToNue->SetLineColor(kGray+3);
  // hErecoMode_NumuBarToNueBar->SetFillColor(0);
  // hErecoMode_NumuBarToNueBar->SetLineColor(kGreen);
  // hErecoData->DrawCopy("pe");
  // //scale hint1 to the pad coordinates
  // float myscale = 3.;
  // Float_t rightmax = 1.1*hErecoData->GetMaximum();
  // Float_t scale = gPad->GetUymax()/rightmax * myscale;
  // hErecoMode_NumuBarToNueBar->Scale(scale);
  // hErecoMode_NumuBarToNueBar->Draw("same");
  // hErecoMode_BkgOther->Scale(scale);
  // hErecoMode_BkgOther->DrawCopy("same");
  // hErecoMode_NC->Scale(scale);
  // hErecoMode_NC->DrawCopy("same");
  // hErecoMode_NumuToNue->Scale(scale);
  // hErecoMode_NumuToNue->DrawCopy("same");
  // cout << "rightmax = " << rightmax << endl;
  // cout << "gPad->GetUxmax() = " << gPad->GetUxmax() << endl;  
  // cout << "scale = " << scale << endl;
  // //draw an axis on the right side
  // //TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
  // TGaxis *axis = new TGaxis(1.25,gPad->GetUymin(),1.25, gPad->GetUymax()*myscale,0,rightmax,510,"+L");
  // //TGaxis *axis = new TGaxis(1.25,gPad->GetUymin(),1.25, gPad->GetUymax()*myscale,0,rightmax,510,"+L");
  // axis->SetLineColor(kBlack);
  // axis->SetLabelColor(kBlack);
  // axis->SetTitle("# of observed events");
  // axis->ImportAxisAttributes(hErecoData->GetYaxis());
  // axis->Draw();  
  // gPad->SetTicks(1,0); // remove ticks on the left
  // return;







  TCanvas *cEreco_ExpVsData_Component = new TCanvas("cEreco_ExpVsData_Component","cEreco_ExpVsData_Component");
  hErecoMode_NC->SetFillColor(0);
  hErecoMode_NC->SetLineColor(kRed);
  hErecoMode_BkgOther->SetFillColor(0);
  hErecoMode_BkgOther->SetLineColor(kBlue);
  hErecoMode_NumuToNue->SetFillColor(0);
  hErecoMode_NumuToNue->SetLineColor(kGray+3);
  hErecoMode_NumuBarToNueBar->SetFillColor(0);
  hErecoMode_NumuBarToNueBar->SetLineColor(kGreen);
  hEreco_AllBkg->SetFillColor(0);
  hEreco_AllBkg->SetLineColor(kBlack);
  hEreco_AllBkg->SetLineStyle(2);
  hEreco_AllBkg->SetLineWidth(3);
  // rescale MC
  hErecoMode_NumuBarToNueBar->Scale(7.);
  hErecoMode_BkgOther->Scale(7.);
  hErecoMode_NC->Scale(7.);
  hErecoMode_NumuToNue->Scale(7.);
  hEreco_AllBkg->Scale(7.);
  //
  hErecoData->DrawCopy("pe");
  hErecoMode_NumuBarToNueBar->Draw("same");
  hErecoMode_BkgOther->DrawCopy("same");
  hErecoMode_NC->DrawCopy("same");
  hErecoMode_NumuToNue->DrawCopy("same");
  hEreco_AllBkg->DrawCopy("same");
  TLegend *leg_comp = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg_comp->SetFillColor(0);
  leg_comp->SetBorderSize(0);
  leg_comp->AddEntry(hErecoData,"Data","P");
  leg_comp->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","L");
  leg_comp->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","L");
  leg_comp->AddEntry(hErecoMode_NC,"NC","L");
  leg_comp->AddEntry(hErecoMode_BkgOther,"Bkg other","L");
  leg_comp->AddEntry(hEreco_AllBkg,"All Bkg","L");
  leg_comp->Draw();
  TLatex *tex_comp = new TLatex(0.7,0.7,"MC #times 7");
  tex_comp->SetTextSize(0.07);
  tex_comp->Draw();
  
   

  TCanvas *Ereco_ExpVsData_Component_Norm = new TCanvas("Ereco_ExpVsData_Component_Norm","Ereco_ExpVsData_Component_Norm");
  hErecoMode_NC_norm->SetFillColor(0);
  hErecoMode_NC_norm->SetLineColor(kRed);
  hErecoMode_BkgOther_norm->SetFillColor(0);
  hErecoMode_BkgOther_norm->SetLineColor(kBlue);
  hErecoMode_NumuToNue_norm->SetFillColor(0);
  hErecoMode_NumuToNue_norm->SetLineColor(kGray+3);
  hErecoMode_NumuBarToNueBar_norm->SetFillColor(0);
  hErecoMode_NumuBarToNueBar_norm->SetLineColor(kGreen);
  hEreco_AllBkg_norm->SetFillColor(0);
  hEreco_AllBkg_norm->SetLineColor(kBlack);
  hEreco_AllBkg_norm->SetLineStyle(2);
  hEreco_AllBkg_norm->SetLineWidth(3);
  // rescale MC
  hErecoMode_NumuBarToNueBar_norm->Scale(7.);
  hErecoMode_BkgOther_norm->Scale(7.);
  hErecoMode_NC_norm->Scale(7.);
  hErecoMode_NumuToNue_norm->Scale(7.);
  hEreco_AllBkg_norm->Scale(7.);
  //
  hErecoData->DrawCopy("pe");
  hErecoMode_NumuBarToNueBar_norm->DrawCopy("same");
  hErecoMode_BkgOther_norm->DrawCopy("same");
  hErecoMode_NC_norm->DrawCopy("same");
  hErecoMode_NumuToNue_norm->DrawCopy("same");
  hEreco_AllBkg_norm->DrawCopy("same");
  TLegend *leg_comp_norm = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg_comp_norm->SetFillColor(0);
  leg_comp_norm->SetBorderSize(0);
  leg_comp_norm->AddEntry(hErecoData,"Data","P");
  leg_comp_norm->AddEntry(hErecoMode_NumuBarToNueBar_norm,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","L");
  leg_comp_norm->AddEntry(hErecoMode_NumuToNue_norm,"#nu_{#mu} #rightarrow #nu_{e}","L");
  leg_comp_norm->AddEntry(hErecoMode_NC_norm,"NC","L");
  leg_comp_norm->AddEntry(hErecoMode_BkgOther_norm,"Bkg other","L");
  leg_comp_norm->AddEntry(hEreco_AllBkg_norm,"All Bkg","L");
  leg_comp_norm->Draw();
 








}
