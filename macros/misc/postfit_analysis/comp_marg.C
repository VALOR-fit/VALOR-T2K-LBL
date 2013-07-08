//
//  Compare marginalised chisq from different analyses
//  Runs on the output of plot_marginalised_chisq.C which writes out marginalised dchisq distributions,
//  for fits to various datasets from a given analysis, like:
//	  KEY: TH1D	marg_sinsq2theta23_0;1	
//	  KEY: TH1D	marg_dm32sq_0;1	
//	  KEY: TH1D	marg_sinsq2theta23_1;1	
//	  KEY: TH1D	marg_dm32sq_1;1	
//	  KEY: TH1D	marg_sinsq2theta23_2;1	
//	  KEY: TH1D	marg_dm32sq_2;1	
//  The number is a dataset id.
//

void comp_marg
(
	const char * file1,
	const char * label1 = "analysis a",
	const char * file2,
	const char * label2 = "analysis b",
        const char * type, 
        int id,
	double ymax
)

{
  gStyle->SetOptStat(0);

  const char * hname = Form("marg_%s_%d",type,id);
  cout << hname << endl;

  TFile * f1 = new TFile(file1,"read");
  f1->cd();
  TH1D * h1 = (TH1D *) f1->Get(hname);
  assert(h1);

  TFile * f2 = new TFile(file2,"read");
  f2->cd();
  TH1D * h2 = (TH1D *) f2->Get(hname);
  assert(h2);

  TCanvas * c1 = new TCanvas("c1","",10,10,800,500);

  h1->SetLineWidth(2);
  h1->SetLineColor(kBlack);
  h1->SetLineStyle(kDashed);

  h2->SetLineWidth(2);
  h2->SetLineColor(kRed);
  h2->SetLineStyle(kSolid);

  h1->GetYaxis()->SetTitle("#Delta#chi^{2}");
  if(type == "dm32sq") {
     h1->GetXaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
     h1->GetXaxis()->SetNdivisions(505);
  } 
  else
  if(type == "sinsq2theta23") {
     h1->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  }

  h1->GetYaxis()->SetRangeUser(0,ymax);

  h1->Draw("L");
  h2->Draw("Lsame");


  TLegend * legend = new TLegend(0.5, 0.6, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->AddEntry(h1,label1,"L");
  legend->AddEntry(h2,label2,"L");
  legend->Draw();

  c1->SaveAs(Form("marg_dchisq_comp_%s.eps",type));
}
