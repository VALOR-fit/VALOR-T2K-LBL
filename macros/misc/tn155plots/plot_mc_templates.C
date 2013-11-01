//run test_nsk.C with:
//  #define __save_templates_debug__ 
//  #define PHYSHYPWEIGHTSAREONE
//  //#define PHYSHYPWEIGHTAREONEANDNUEISZERO
//
// Template_mcnorm_fluxtuned_banfftuned_syst_phys.root is what is used here

void plot_mc_templates() 
{
  using namespace t2kana;

  TFile f("Template_mcnorm_fluxtuned_banfftuned_syst_phys.root");
  double templemax = 10.0;
  bool debug = false; // check the templates are the correct ones

  //get the templates
  const int ntempl = 32;
  TH2D * h_template[ntempl] = {0};
  for(int i = 0; i < ntempl; i++) {
    f.GetObject(Form("nsk_fin_0_%d",i),h_template[i]);
    //h_template[i]->SetTitle("");
    h_template[i]->GetXaxis()->SetRangeUser(0,templemax-conv::kOneMeV);
    h_template[i]->GetYaxis()->SetRangeUser(0,templemax-conv::kOneMeV);
    h_template[i]->GetYaxis()->SetTitleOffset(0.65);
    if(debug)
      h_template[i]->GetXaxis()->SetTitle(h_template[i]->GetTitle());
  }

  TCanvas * ctempl = new TCanvas("ctempl","ctempl",0,0,700,750);
  ctempl->Divide(2, 3);
  TCanvas * ctemplosc = new TCanvas();
  ctemplosc->Divide(2,2);
  int ic = 0; //number on the 6x6 canvas
  int icosc = 0; //number on the 4x4 canvas
  int iplot = 0; // number of the plot
  for(int it = 0; it < ntempl; it++) {
    //oscillated (CC only)
    if((it >= 6  && it <= 9) || (it >= 28 && it <= 31)) {
      ctemplosc->cd((icosc%4) +1)->SetRightMargin(0.15);
      ctemplosc->cd((icosc%4) +1)->SetLeftMargin(0.09);
      h_template[it]->Draw("colz");
      if(icosc%4 == 3) {
	ctemplosc->SaveAs(Form("mc_template_%d.pdf",iplot));
	ctemplosc->SaveAs(Form("mc_template_%d.eps",iplot));
	iplot++;
      }
      icosc++;
    }//osc
    //not oscillated (CC+NC)
    else {
      ctempl->cd((ic%6) +1)->SetRightMargin(0.12);
      ctempl->cd((ic%6) +1)->SetLeftMargin(0.08);
      ctempl->cd((ic%6) +1)->SetBottomMargin(0.13);
      ctempl->cd((ic%6) +1)->SetTopMargin(0.06);
      h_template[it]->Draw("colz");
      if(ic%6 == 5) {
	ctempl->SaveAs(Form("mc_template_%d.pdf",iplot));
	ctempl->SaveAs(Form("mc_template_%d.eps",iplot));
	iplot++;
      }
      ic++;
    }//not osc
  }//loop over templates
}

