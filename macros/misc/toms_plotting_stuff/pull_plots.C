{
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.15);
  TChain ch("nfitnt");
  ch.Add("28/pull28all-*root");
  cout << ch.GetEntries() << endl;
  TH2D * h2 = new TH2D("h2","h2",20,-5,5,40,-5,5);
  TH1D * hsame1 = new TH1D("hsame1","hsame1",250,-5,5);
  TH1D * hsame2 = new TH1D("hsame2","hsame2",250,-5,5);
  TH1D * hsf1 = new TH1D("hsf1","hsf1",200,-1,1);
  TH1D * hsf2 = new TH1D("hsf2","hsf2",200,-1.5,1);
  hsame2->SetLineColor(kRed);
  int x,y,z;

  gStyle->SetOptStat(0);
						
  //compare pull1 with pull2 (for ccnqe skdet)
  hsame1->GetXaxis()->SetTitle("pull");
  ch.Draw("(bestfit_skdet_ccnqe-input_skdet_ccnqe)/errmigrad_skdet_ccnqe>>hsame1","migradfitstatus!=4&&fitstatus==3");
  ch.Draw("(bestfit_skdet_ccnqe-input_skdet_ccnqe)/(sqrt(1-pow(errmigrad_skdet_ccnqe,2)))>>hsame2","(migradfitstatus!=4)&&fitstatus==3&&(errmigrad_skdet_ccnqe<1)","SAME");
  TLegend * leg = new TLegend(0.2,0.8,0.4,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hsame1,"pull_{1}","l");
  leg->AddEntry(hsame2,"pull_{2}","l");
  leg->Draw();
  c->SaveAs("28_all_0.900000_skdet_ccnqe_twopulls_rev.eps");

  //bestfit of SF
  hsf1->GetXaxis()->SetTitle("bestfit");
  ch.Draw("bestfit_xsec_sf>>hsf1","(migradfitstatus!=4)&&(fitstatus==3)");
  TLegend * leg2 = new TLegend(-1,0,0,hsf1->GetMaximum()*1.05,"","br");
  leg2->SetFillStyle(3004);
  leg2->SetFillColor(kRed);
  leg2->AddEntry((TObject*)0,"Not","");
  leg2->AddEntry((TObject*)0,"Allowed","");
  leg2->Draw();
  c->SaveAs("28_all_0.900000_sf_bestfit_rev.eps");
  //pull1 of SF
  hsf2->GetXaxis()->SetTitle("pull_{1}");
  ch.Draw("(bestfit_xsec_sf-input_xsec_sf)/errmigrad_xsec_sf>>hsf2","(migradfitstatus!=4)&&(fitstatus==3)");
  TLine l(0,0,0,hsf2->GetMaximum()*1.05);
  l->SetLineStyle(2);
  l->Draw();
  c->SaveAs("28_all_0.900000_sf_pull_rev.eps");
  
  gStyle->SetOptStat(1100);
  c->SetTopMargin(0.1);

  //example of 2D mixing vs syst
  h2->GetXaxis()->SetTitle("f_{CCnQE; r}^{SK} residual");
  h2->GetYaxis()->SetTitle("sin^{2}#theta2_{23} pull_{1}");
  ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_skdet_ccnqe-input_skdet_ccnqe>>h2","migradfitstatus!=4&&fitstatus==3","colz");
  c->SaveAs("28_all_0.900000_mixpull_skdetccnqeres_rev.eps");


  // unrequired plots

  /*
    ch.Draw("errmigrad_skdet_energyscale>>hnew:(120,0,1.2),"migradfitstatus!=4&&fitstatus==3"");
    c->SaveAs("28_all_0.900000_escale_err.eps");
    ch.Draw("errmigrad_skdet_nc>>hnew:(120,0,1.2),"migradfitstatus!=4&&fitstatus==3"");
    c->SaveAs("28_all_0.900000_nc_err.eps");

    //mix pull vs syst residual
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_xsec_sf-input_xsec_sf>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_sfres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_xsec_pilessdelta-input_xsec_pilessdelta>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_pddres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_xsec_binde-input_xsec_binde>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_binderes.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_skdet_nc-input_skdet_nc>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_ncres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_skdet_ccnqe-input_skdet_ccnqe>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_ccnqeres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_skdet_energyscale-input_skdet_energyscale>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_escaleres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_banff_maqe-input_banff_maqe>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_maqeres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:bestfit_banff_ccqeE0-input_banff_ccqeE0>>h2","migradfitstatus!=4&&fitstatus==3","colz");
    h2->GetBinXYZ(h2->GetMaximumBin(),x,y,z);
    cout << h2->GetYaxis()->GetBinCenter(y) <<  "\t" << h2->GetXaxis()->GetBinCenter(x) << endl;
    c->SaveAs("mixpull_ccqee0res.pdf");

    //mix residual vs syst residual
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_xsec_sf-input_xsec_sf>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_sfres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_xsec_pilessdelta-input_xsec_pilessdelta>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_pddres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_xsec_binde-input_xsec_binde>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_binderes.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_skdet_nc-input_skdet_nc>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_ncres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_skdet_ccnqe-input_skdet_ccnqe>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_ccnqeres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_skdet_energyscale-input_skdet_energyscale>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_escaleres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_banff_maqe-input_banff_maqe>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_maqeres.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):bestfit_banff_ccqeE0-input_banff_ccqeE0>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_ccqee0res.pdf");

    //mix pull vs syst pull_osc
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_xsec_sf-input_xsec_sf)/errmigrad_xsec_sf>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_sfpull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_xsec_pilessdelta-input_xsec_pilessdelta)/errmigrad_xsec_pilessdelta>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_pddpull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_xsec_binde-input_xsec_binde)/errmigrad_xsec_binde>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_bindepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_skdet_nc-input_skdet_nc)/errmigrad_skdet_nc>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_ncpull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_skdet_ccnqe-input_skdet_ccnqe)/errmigrad_skdet_ccnqe>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_ccnqepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_skdet_energyscale-input_skdet_energyscale)/errmigrad_skdet_energyscale>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_escalepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_banff_maqe-input_banff_maqe)/errmigrad_banff_maqe>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_maqepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23:(bestfit_banff_ccqeE0-input_banff_ccqeE0)/errmigrad_banff_ccqeE0>>hnew(20,-5,5,40,-5,5)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixpull_ccqee0pull.pdf");

    //mix residual vs syst pull_osc
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_xsec_sf-input_xsec_sf)/errmigrad_xsec_sf>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_sfpull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_xsec_pilessdelta-input_xsec_pilessdelta)/errmigrad_xsec_pilessdelta>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_pddpull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_xsec_binde-input_xsec_binde)/errmigrad_xsec_binde>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_bindepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_skdet_nc-input_skdet_nc)/errmigrad_skdet_nc>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_ncpull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_skdet_ccnqe-input_skdet_ccnqe)/errmigrad_skdet_ccnqe>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_ccnqepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_skdet_energyscale-input_skdet_energyscale)/errmigrad_skdet_energyscale>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_escalepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_banff_maqe-input_banff_maqe)/errmigrad_banff_maqe>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_maqepull.pdf");
    ch.Draw("(bestfit_mixing23-input_mixing23):(bestfit_banff_ccqeE0-input_banff_ccqeE0)/errmigrad_banff_ccqeE0>>hnew(20,-5,5,40,-0.25,0.25)","migradfitstatus!=4&&fitstatus==3","colz");
    c->SaveAs("mixres_ccqee0pull.pdf");


    ch.Draw("(bestfit_mixing23-input_mixing23)/errmigrad_mixing23>>hnew(100,-5,5)","migradfitstatus!=4&&fitstatus==3","");
  */

}
