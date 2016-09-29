
//
// ------------------------------------------------------------------------------------------------
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareContour1D
(
 const char* file_1 = "",
 const char* file_2 = "",
 const char* tag    = "",
 const char* tag1   = "",
 const char* tag2   = "",
 
 const bool doGlobalMin  = false,
 const bool same_mh_color = true, // use same color for NH and IH
 const bool same_mh_style = false, // use same style for NH and IH
 const bool draw_preliminary = false,
 
 const bool draw_1sigmaCL = false,
 const bool draw_90CL     = false,
 const bool draw_2sigmaCL = false,
 const bool draw_3sigmaCL = false,
 const bool draw_4sigmaCL = false,
 const bool draw_5sigmaCL = true,
 
 const bool draw_reactors = false
 )
{
  TGaxis::SetMaxDigits(3);

  const int Nfiles = 2;
  TString latexname[Nfiles] = {
    "CC#pi",
    "1Ring"
  };
  
  if(tag1 != "") latexname[0] = tag1;
  if(tag2 != "") latexname[1] = tag2;
  
  TCanvas *ccontours = new TCanvas("ccontours","ccontours",1200,900);

  TFile f1(file_1,"read");
  TFile f2(file_2,"read");
  
  TH1D* h1 = f1.Get("h_deltachisq");
  TH1D* h2 = f2.Get("h_deltachisq");
  h1->SetLineColor(kRed);
  
  h1->GetYaxis()->SetRangeUser(0,7);
  h1->Draw();
  h2->SetLineStyle(2);
  h2->Draw("same");
  
  TLine *l1sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),1.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),1.);
  TLine *l90 = new TLine(h1->GetXaxis()->GetBinLowEdge(1),2.71,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),4.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),4.);
  TLine *l3sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),9.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),9.);
  TLine *l4sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),16.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),16.);
  TLine *l5sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),25.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),25.);
  l1sigma->SetLineWidth(2);
  l90->SetLineWidth(2);
  l2sigma->SetLineWidth(2);
  l3sigma->SetLineWidth(2);
  l4sigma->SetLineWidth(2);
  l5sigma->SetLineWidth(2);
  l1sigma->SetLineColor(kBlue);
  l90->SetLineColor(kBlue);
  l2sigma->SetLineColor(kBlue);
  l3sigma->SetLineColor(kBlue);
  l4sigma->SetLineColor(kBlue);
  l5sigma->SetLineColor(kBlue);
  if(draw_1sigmaCL) l1sigma->Draw();
  if(draw_90CL)     l90->Draw();
  if(draw_2sigmaCL) l2sigma->Draw();
  if(draw_3sigmaCL) l3sigma->Draw();
  if(draw_4sigmaCL) l4sigma->Draw();
  if(draw_5sigmaCL) l5sigma->Draw();
  
  TLatex *latLine1Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),4.1,"2#sigma");
  TLatex *latLine3Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),9.1,"3#sigma");
  TLatex *latLine4Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),16.1,"4#sigma");
  TLatex *latLine5Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),25.1,"5#sigma");
  latLine1Sigma->SetTextColor(kBlue);
  latLine90CL->SetTextColor(kBlue);
  latLine2Sigma->SetTextColor(kBlue);
  latLine3Sigma->SetTextColor(kBlue);
  latLine4Sigma->SetTextColor(kBlue);
  latLine5Sigma->SetTextColor(kBlue);
  latLine1Sigma->SetTextSize(0.07);
  latLine90CL->SetTextSize(0.07);
  latLine2Sigma->SetTextSize(0.07);
  latLine3Sigma->SetTextSize(0.07);
  latLine4Sigma->SetTextSize(0.07);
  latLine5Sigma->SetTextSize(0.07);
  if(draw_1sigmaCL) latLine1Sigma->Draw("same");
  if(draw_90CL)     latLine90CL->Draw("same");
  if(draw_2sigmaCL) latLine2Sigma->Draw("same");
  if(draw_3sigmaCL) latLine3Sigma->Draw("same");
  if(draw_4sigmaCL) latLine4Sigma->Draw("same");
  if(draw_5sigmaCL) latLine5Sigma->Draw("same");
  TLegend *lcontours = new TLegend(0.65,0.73,0.85,0.9);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  lcontours->AddEntry(h1,latexname[0]);
  lcontours->AddEntry(h2,latexname[1]);
  lcontours->Draw();
  
  TLatex *latPreliminary = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),1.1,"PRELIMINARY");
  if(draw_preliminary) latPreliminary->Draw("same");
  ccontours->SaveAs(Form("%s.pdf",tag));
  // // Store in root file
  // if(SaveRootFile){
  //   TFile* f2 = new TFile(Form("%s_contour.root",tag),"recreate");
  //   for(int ifile=0;ifile<Nfiles;ifile++){
  //     for(int imh=0;imh<Nmh;imh++){
  // 	//sprintf(name,"hdeltachisq_%i_%i",ifile);
  // 	hdeltachisq1D_GlobMin[ifile]->Write();//(name);
  //     }
  //   }
  //   f2->Close();
  // }
  
  return;
}


