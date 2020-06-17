void make_plots_total_inputs(){

  TFile *_file0 = TFile::Open("data/dune/ndfit/2016b/flux/2016b_flux_DUNE_optimised.root");
  TH2D * h = new TH2D("h","Fractional covariance matrix", 264,0.,264,264,0,264);

  for(int i=1; i<209; i++){
    for(int j=1; j<209; j++){
      double value = total_covariance->GetBinContent(i,j);
      std::cout << "value = " << value << std::endl;
      h->SetBinContent(i+1,j+1,value);
    }
  }
  TFile *_file1 = TFile::Open("data/dune/ndfit/2016b/xsec/new_xsec_covariance_2016b_LBNEND.root");
  _file1.cd();

  for(int i=0; i<43; i++){
    for(int j=0; j<43; j++){
      double value = frac_cov_max(i,j);
      h->SetBinContent(i+209,j+209,value);
    }
  }

  for(int i=251; i<264; i++)
    h->SetBinContent(i+1,i+1,0.01);
  
  //  for(int i=1; i<157; i++)
  //{
      double sigma = 100*sqrt(h->GetBinContent(i,i));
      //    std::cout << "One sigma error for: " << i << " = " << sigma << std::endl;
      //}



  //  TCanvas *c0 = new TCanvas("c0","c0",500,500);
  //gStyle->SetOptStat(0);
  //lartpc_cov->Draw("colz");

  TCanvas *c1 = new TCanvas("c1","c1",750,750);
  c1->SetLogz();
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 

  //  h->GetXaxis()->SetLabelSize("0.025");
  h->GetXaxis()->SetBinLabel(1,"FHC #nu_{#mu}");
  h->GetXaxis()->SetBinLabel(20,"FHC #bar{#nu}_{#mu}");
  h->GetXaxis()->SetBinLabel(39,"FHC #nu_{e}");
  h->GetXaxis()->SetBinLabel(46,"FHC #bar{#nu}_{e}");
  h->GetXaxis()->SetBinLabel(53,"RHC #nu_{#mu}");
  h->GetXaxis()->SetBinLabel(72,"RHC #bar{#nu}_{#mu}");
  h->GetXaxis()->SetBinLabel(91,"RHC #nu_{e}");
  h->GetXaxis()->SetBinLabel(98,"RHC #bar{#nu}_{e}");

  h->GetXaxis()->SetBinLabel(105,"FHC #nu_{#mu}");
  h->GetXaxis()->SetBinLabel(124,"FHC #bar{#nu}_{#mu}");
  h->GetXaxis()->SetBinLabel(143,"FHC #nu_{e}");
  h->GetXaxis()->SetBinLabel(150,"FHC #bar{#nu}_{e}");
  h->GetXaxis()->SetBinLabel(157,"RHC #nu_{#mu}");
  h->GetXaxis()->SetBinLabel(176,"RHC #bar{#nu}_{#mu}");
  h->GetXaxis()->SetBinLabel(195,"RHC #nu_{e}");
  h->GetXaxis()->SetBinLabel(202,"RHC #bar{#nu}_{e}");

  /*  h->GetXaxis()->SetBinLabel(1+104,"CCQE");
  //  h->GetXaxis()->SetBinLabel(4+104,"#bar{#nu}_{#mu} CCQE");
  h->GetXaxis()->SetBinLabel(7+104,"MEC");
  //  h->GetXaxis()->SetBinLabel(8+104,"#bar{#nu}_{#mu} MEC");
  h->GetXaxis()->SetBinLabel(9+104,"CC1#pi^{#pm}");
  //  h->GetXaxis()->SetBinLabel(12+104,"#bar{#nu}_{#mu} CC1#pi^{#pm}");
  h->GetXaxis()->SetBinLabel(15+104,"CC1#pi^{0}");
  //  h->GetXaxis()->SetBinLabel(18+104,"#bar{#nu}_{#mu} CC1#pi^{0}");
  h->GetXaxis()->SetBinLabel(21+104,"CC2#pi");
  h->GetXaxis()->SetBinLabel(23+104,"DIS");
  //  h->GetXaxis()->SetBinLabel(26+104,"#bar{#nu}_{#mu} DIS");
  h->GetXaxis()->SetBinLabel(29+104,"CC coh");
  h->GetXaxis()->SetBinLabel(31+104,"NC");
  //  h->GetXaxis()->SetBinLabel(33+104,"CC #nu_{e}/#nu_{#mu}");
  //  h->GetXaxis()->SetBinLabel(34+104,"MFP #pi");
  // h->GetXaxis()->SetBinLabel(35+104,"MFP N");
  h->GetXaxis()->SetBinLabel(36+104,"#pi FSI");
  //  h->GetXaxis()->SetBinLabel(37+104,"FrInel #pi");
  //h->GetXaxis()->SetBinLabel(38+104,"FrAbs #pi");
  //h->GetXaxis()->SetBinLabel(39+104,"FrPiProd #pi");
  h->GetXaxis()->SetBinLabel(40+104,"N FSI");
  //  h->GetXaxis()->SetBinLabel(41+104,"FrInel N");
  //h->GetXaxis()->SetBinLabel(42+104,"FrAbs N");
  //h->GetXaxis()->SetBinLabel(43+104,"FrPiProd N");
  //h->GetXaxis()->LabelsOption("v");
  */
  h->GetXaxis()->SetBinLabel(252,"ND eff");//"Eff #nu_{#mu} CCQE");
  h->GetXaxis()->SetBinLabel(260,"FD eff");//"Eff #nu_{#mu} CCQE");
  /*  h->GetXaxis()->SetBinLabel(150,"Eff #nu_{#mu} CC#pi");
  //  h->GetXaxis()->SetBinLabel(152,"Eff #nu_{#mu} CC2#pi");
  h->GetXaxis()->SetBinLabel(153.5,"Eff #nu_{e}");
  h->GetXaxis()->SetBinLabel(155.5,"Eff wrong sign");
  //h->GetXaxis()->SetBinLabel(156.5,"Eff NC");
  */
  h->GetXaxis()->SetTickLength(0.);
  h->GetXaxis()->LabelsOption("v");


  //Y AXIS
  //  h->GetYaxis()->SetLabelSize("0.025");
  h->GetYaxis()->SetBinLabel(1,"FHC #nu_{#mu}");
  h->GetYaxis()->SetBinLabel(20,"FHC #bar{#nu}_{#mu}");
  h->GetYaxis()->SetBinLabel(39,"FHC #nu_{e}");
  h->GetYaxis()->SetBinLabel(46,"FHC #bar{#nu}_{e}");
  h->GetYaxis()->SetBinLabel(53,"RHC #nu_{#mu}");
  h->GetYaxis()->SetBinLabel(72,"RHC #bar{#nu}_{#mu}");
  h->GetYaxis()->SetBinLabel(91,"RHC #nu_{e}");
  h->GetYaxis()->SetBinLabel(98,"RHC #bar{#nu}_{e}");

  h->GetYaxis()->SetBinLabel(105,"FHC #nu_{#mu}");
  h->GetYaxis()->SetBinLabel(124,"FHC #bar{#nu}_{#mu}");
  h->GetYaxis()->SetBinLabel(143,"FHC #nu_{e}");
  h->GetYaxis()->SetBinLabel(150,"FHC #bar{#nu}_{e}");
  h->GetYaxis()->SetBinLabel(157,"RHC #nu_{#mu}");
  h->GetYaxis()->SetBinLabel(176,"RHC #bar{#nu}_{#mu}");
  h->GetYaxis()->SetBinLabel(195,"RHC #nu_{e}");
  h->GetYaxis()->SetBinLabel(202,"RHC #bar{#nu}_{e}");

  h->GetYaxis()->SetBinLabel(252,"Near detector eff");//"Eff #nu_{#mu} CCQE"); 
  h->GetYaxis()->SetBinLabel(260,"Far detector eff");//"Eff #nu_{#mu} CCQE");  

  /*  h->GetYaxis()->SetBinLabel(10,"FHC #nu_{#mu}");
  h->GetYaxis()->SetBinLabel(30,"FHC #bar{#nu}_{#mu}");
  h->GetYaxis()->SetBinLabel(42,"FHC #nu_{e}");
  h->GetYaxis()->SetBinLabel(49,"FHC #bar{#nu}_{e}");
  h->GetYaxis()->SetBinLabel(63,"RHC #nu_{#mu}");
  h->GetYaxis()->SetBinLabel(82,"RHC #bar{#nu}_{#mu}");
  h->GetYaxis()->SetBinLabel(95,"RHC #nu_{e}");
  h->GetYaxis()->SetBinLabel(101,"RHC #bar{#nu}_{e}");

  h->GetYaxis()->SetBinLabel(4+104,"CCQE");
  //  h->GetYaxis()->SetBinLabel(4+104,"#bar{#nu}_{#mu} CCQE");
  h->GetYaxis()->SetBinLabel(8+104,"MEC");
  //  h->GetYaxis()->SetBinLabel(8+104,"#bar{#nu}_{#mu} MEC");
  h->GetYaxis()->SetBinLabel(12+104,"CC1#pi^{#pm}");
  //  h->GetYaxis()->SetBinLabel(12+104,"#bar{#nu}_{#mu} CC1#pi^{#pm}");
  h->GetYaxis()->SetBinLabel(18+104,"CC1#pi^{0}");
  //  h->GetYaxis()->SetBinLabel(18+104,"#bar{#nu}_{#mu} CC1#pi^{0}");
  h->GetYaxis()->SetBinLabel(21+104,"CC2#pi");
  h->GetYaxis()->SetBinLabel(26+104,"DIS");
  //  h->GetYaxis()->SetBinLabel(26+104,"#bar{#nu}_{#mu} DIS");
  h->GetYaxis()->SetBinLabel(30+104,"CC coh");
  h->GetYaxis()->SetBinLabel(32+104,"NC");
  //  h->GetYaxis()->SetBinLabel(33+104,"CC #nu_{e}/#nu_{#mu}");
  //  h->GetYaxis()->SetBinLabel(34+104,"MFP #pi");
  // h->GetYaxis()->SetBinLabel(35+104,"MFP N");
  h->GetYaxis()->SetBinLabel(38+104,"#pi FSI");
  //  h->GetYaxis()->SetBinLabel(37+104,"FrInel #pi");
  //h->GetYaxis()->SetBinLabel(38+104,"FrAbs #pi");
  //h->GetYaxis()->SetBinLabel(39+104,"FrPiProd #pi");
  h->GetYaxis()->SetBinLabel(43+104,"N FSI");
  //  h->GetYaxis()->SetBinLabel(41+104,"FrInel N");
  //h->GetYaxis()->SetBinLabel(42+104,"FrAbs N");
  //h->GetYaxis()->SetBinLabel(43+104,"FrPiProd N");
  //h->GetYaxis()->LabelsOption("v");

  h->GetYaxis()->SetBinLabel(148,"Eff #nu_{#mu} CCQE");
  h->GetYaxis()->SetBinLabel(150,"Eff #nu_{#mu} CC#pi");
  //  h->GetYaxis()->SetBinLabel(152,"Eff #nu_{#mu} CC2#pi");
  h->GetYaxis()->SetBinLabel(153.5,"Eff #nu_{e}     ");
  h->GetYaxis()->SetBinLabel(156,"Eff wrong sign");
  //  h->GetYaxis()->SetBinLabel(156.5,"Eff NC");

  */
  h->GetYaxis()->SetTickLength(0.);
  h->GetYaxis()->LabelsOption("v");




  h->Draw("colz");
  //  h->GetZaxis()->SetRangeUser(10E-10,10E-6);
  //  h->Draw("colz");
  /*
  TLine * l0 = new TLine(104,0,104,156);
  TLine * l1 = new TLine(0,104,156,104);

  TLine * l2 = new TLine(6+104,0,6+104,156);
  TLine * l3 = new TLine(0,6+104,156,6+104);

  TLine * l4 = new TLine(8+104,0,8+104,156);
  TLine * l5 = new TLine(0,8+104,156,8+104);

  TLine * l6 = new TLine(19,0,19,156);
  TLine * l7 = new TLine(0,19,156,19);

  TLine * l34 = new TLine(14+104,0,14+104,156);
  TLine * l35 = new TLine(0,14+104,156,14+104);

  TLine * l8 = new TLine(20+104,0,20+104,156);
  TLine * l9 = new TLine(0,20+104,156,20+104);

  TLine * l10 = new TLine(22+104,0,22+104,156);
  TLine * l11 = new TLine(0,22+104,156,22+104);

  TLine * l12 = new TLine(38,0,38,156);
  TLine * l13 = new TLine(0,38,156,38);

  TLine * l14 = new TLine(28+104,0,28+104,156);
  TLine * l15 = new TLine(0,28+104,156,28+104);

  TLine * l16 = new TLine(30+104,0,30+104,156);
  TLine * l17 = new TLine(0,30+104,156,30+104);

  TLine * l18 = new TLine(32+104,0,32+104,156);
  TLine * l19 = new TLine(0,32+104,156,32+104);

  TLine * l20 = new TLine(45,0,45,156);
  TLine * l21 = new TLine(0,45,156,45);

  TLine * l22 = new TLine(35+104,0,35+104,156);
  TLine * l23 = new TLine(0,35+104,156,35+104);

  TLine * l24 = new TLine(52,0,52,156);
  TLine * l25 = new TLine(0,52,156,52);

  TLine * l26 = new TLine(71,0,71,156);
  TLine * l27 = new TLine(0,71,156,71);

  TLine * l28 = new TLine(90,0,90,156);
  TLine * l29 = new TLine(0,90,156,90);

  TLine * l30 = new TLine(39+104,0,39+104,156);
  TLine * l31 = new TLine(0,39+104,156,39+104);

    TLine * l32 = new TLine(97,0,97,156);
  TLine * l33 = new TLine(0,97,156,97);

  
  TLine * l36 = new TLine(148,0,148,156);
  TLine * l37 = new TLine(0,148,156,148);

  TLine * l38 = new TLine(150,0,150,156);
  TLine * l39 = new TLine(0,150,156,150);

  TLine * l40 = new TLine(152,0,152,156);
  TLine * l41 = new TLine(0,152,156,152);

  TLine * l44 = new TLine(154,0,154,156);
  TLine * l45 = new TLine(0,154,156,154);

  TLine * l46 = new TLine(155,0,155,156);
  TLine * l47 = new TLine(0,155,156,155);

  TLine * l48 = new TLine(156,0,156,156);
  TLine * l49 = new TLine(0,156,156,156);
  
  TLine * l42 = new TLine(33+104,0,33+104,156);
  TLine * l43 = new TLine(0,33+104,156,33+104);


  l0->SetLineWidth(2);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l3->SetLineWidth(2);
  l4->SetLineWidth(2);
  l5->SetLineWidth(2);
  l6->SetLineWidth(2);
  l7->SetLineWidth(2);
  l8->SetLineWidth(2);
  l9->SetLineWidth(2);
  l10->SetLineWidth(2);
  l11->SetLineWidth(2);
  l12->SetLineWidth(2);
  l13->SetLineWidth(2);
  l14->SetLineWidth(2);
  l15->SetLineWidth(2);
  l16->SetLineWidth(2);
  l17->SetLineWidth(2);
  l18->SetLineWidth(2);
  l19->SetLineWidth(2);
  l20->SetLineWidth(2);
  l21->SetLineWidth(2);
  l22->SetLineWidth(2);
  l23->SetLineWidth(2);
  l24->SetLineWidth(2);
  l25->SetLineWidth(2);
  l26->SetLineWidth(2);
  l27->SetLineWidth(2);
  l28->SetLineWidth(2);
  l29->SetLineWidth(2);
  l30->SetLineWidth(2);
  l31->SetLineWidth(2);
    l32->SetLineWidth(2);
  l33->SetLineWidth(2);
    l34->SetLineWidth(2);
  l35->SetLineWidth(2);
  l36->SetLineWidth(2);
  l37->SetLineWidth(2);
  l38->SetLineWidth(2);
  l39->SetLineWidth(2);
  l40->SetLineWidth(2);
  l41->SetLineWidth(2);
  l42->SetLineWidth(2);
  l43->SetLineWidth(2);
  l44->SetLineWidth(2);
  l45->SetLineWidth(2);
  l46->SetLineWidth(2);
  l47->SetLineWidth(2);
  l48->SetLineWidth(2);
  l49->SetLineWidth(2);


  l0->Draw("same");
  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");
  l4->Draw("same");
  l5->Draw("same");
  l6->Draw("same");
  l7->Draw("same");
  l8->Draw("same");
  l9->Draw("same");
  l10->Draw("same");
  l11->Draw("same");
  l12->Draw("same");
  l13->Draw("same");
  l14->Draw("same");
  l15->Draw("same");
  l16->Draw("same");
  l17->Draw("same");
  l18->Draw("same");
  l19->Draw("same");
  l20->Draw("same");
  l21->Draw("same");
  l22->Draw("same");
  l23->Draw("same");
  l24->Draw("same");
  l25->Draw("same");
  l26->Draw("same");
  l27->Draw("same");
  l28->Draw("same");
  l29->Draw("same");
  l30->Draw("same");
  l31->Draw("same");
  l32->Draw("same");
  l33->Draw("same");
  l34->Draw("same");
  l35->Draw("same");
  l36->Draw("same");
  l37->Draw("same");
  l38->Draw("same");
  l39->Draw("same");
  l40->Draw("same");
  l41->Draw("same");
  l42->Draw("same");
  l43->Draw("same");
  l44->Draw("same");
  l45->Draw("same");
  l46->Draw("same");
  l47->Draw("same");
  l48->Draw("same");
  l49->Draw("same");
  */

}
