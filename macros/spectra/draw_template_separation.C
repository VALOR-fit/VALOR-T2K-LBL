
void draw_template_separation
(
 double energy_up = 2.,
 double energy_down = 0.,
 const bool NormArea = true,
 
 const bool draw_NumuBarToNueBar = true, 
 const bool draw_NumuToNue = true,
 const bool draw_NueBeam = false,
 const bool draw_NumuBeam = false,
 const bool draw_NumuBarBeam = false,
 const bool draw_NC = false,
)
{
  char name[200];

  const int Nfiles = 3;
  char *filename[Nfiles][2] = {
    {"Template_2D/template_ErecoVsTheta_2D_1Re.root", "Theta"},
    {"Template_2D/template_ErecoVsN_2D_1Re.root", "N"    },
    {"Template_2D/template_ErecoVsQ2_2D_1Re.root", "Q2"    },

    //{"Template_2D/template_ErecoVsTheta_2D_RHC1Re.root", "Theta"},
    //{"Template_2D/template_ErecoVsN_2D_RHC1Re.root", "N"    },
    //{"Template_2D/template_ErecoVsQ2_2D_RHC1Re.root", "Q2"    },
  };  
  
  double separation[Nfiles]={0.};
  double separation_weighted[Nfiles]={0.};
  double separation_2D[Nfiles]={0.};
  double separation_weighted_2D[Nfiles]={0.};

  TFile *file[Nfiles];
  TH1D *hQ2Mode_init[Nfiles]; 
  TH1D *hQ2Mode_NumuBarToNueBar[Nfiles]; 
  TH1D *hQ2Mode_NumuToNue[Nfiles];
  TH1D *hQ2Mode_NueBeam[Nfiles]; 
  TH1D *hQ2Mode_NumuBeam[Nfiles];
  TH1D *hQ2Mode_NumuBarBeam[Nfiles];
  TH1D *hQ2Mode_NC[Nfiles]; 

  TH2D *hErecoVsQ2Mode_init[Nfiles]; 
  TH2D *hErecoVsQ2Mode_NumuBarToNueBar[Nfiles]; 
  TH2D *hErecoVsQ2Mode_NumuToNue[Nfiles];
  TH2D *hErecoVsQ2Mode_NueBeam[Nfiles]; 
  TH2D *hErecoVsQ2Mode_NumuBeam[Nfiles];
  TH2D *hErecoVsQ2Mode_NumuBarBeam[Nfiles];
  TH2D *hErecoVsQ2Mode_NC[Nfiles]; 

  TH1D *hSeparation_Q2[Nfiles];
  TH2D *hSeparation_ErecoVsQ2[Nfiles];
  
  TCanvas *cSeparation_Q2[Nfiles];
  TCanvas *cNuNubarQ2[Nfiles]; 
  
  TCanvas *cErecoVsQ2_NumuBarToNueBar[Nfiles]; 
  TCanvas *cErecoVsQ2_NumuToNue[Nfiles]; 

  TLegend *legX[Nfiles];

  for(int ifile=0;ifile<Nfiles;ifile++){
    
    file[ifile] = new TFile(filename[ifile][0],"READ");
    if(!file[ifile]->IsOpen()){
      cerr << endl;
      cerr << "File " << ifile << " --> " << filename[ifile][1] << " is not open!" << endl;
      cerr << endl;
      exit(1);
    }

    // Compare the components 
    
    hQ2Mode_NumuBarToNueBar[ifile] = (TH1D*) file[ifile]->Get("hQ2Mode_NumuBarToNueBar");
    hQ2Mode_NumuToNue[ifile]       = (TH1D*) file[ifile]->Get("hQ2Mode_NumuToNue");
    hQ2Mode_NueBeam[ifile]         = (TH1D*) file[ifile]->Get("hQ2Mode_NueBeam");
    hQ2Mode_NumuBeam[ifile]        = (TH1D*) file[ifile]->Get("hQ2Mode_NumuBeam");
    hQ2Mode_NumuBarBeam[ifile]     = (TH1D*) file[ifile]->Get("hQ2Mode_NumuBarBeam");
    hQ2Mode_NC[ifile]              = (TH1D*) file[ifile]->Get("hQ2Mode_NC");
    hQ2Mode_init[ifile] = (TH1D*) hQ2Mode_NC[ifile]->Clone("hQ2Mode_init");
    hQ2Mode_init[ifile]->Reset();

    
    hErecoVsQ2Mode_NumuBarToNueBar[ifile] = (TH2D*) file[ifile]->Get("hErecoVsQ2Mode_NumuBarToNueBar");
    hErecoVsQ2Mode_NumuToNue[ifile]       = (TH2D*) file[ifile]->Get("hErecoVsQ2Mode_NumuToNue");
    hErecoVsQ2Mode_NueBeam[ifile]         = (TH2D*) file[ifile]->Get("hErecoVsQ2Mode_NueBeam");
    hErecoVsQ2Mode_NumuBeam[ifile]        = (TH2D*) file[ifile]->Get("hErecoVsQ2Mode_NumuBeam");
    hErecoVsQ2Mode_NumuBarBeam[ifile]     = (TH2D*) file[ifile]->Get("hErecoVsQ2Mode_NumuBarBeam");
    hErecoVsQ2Mode_NC[ifile]              = (TH2D*) file[ifile]->Get("hErecoVsQ2Mode_NC");
    hErecoVsQ2Mode_init[ifile] = (TH2D*) hErecoVsQ2Mode_NC[ifile]->Clone("hErecoVsQ2Mode_init");
    hErecoVsQ2Mode_init[ifile]->Reset();

    
    if(NormArea){
      if(draw_NumuBarToNueBar){
	hQ2Mode_NumuBarToNueBar[ifile]->Scale(1. / hQ2Mode_NumuBarToNueBar[ifile]->Integral() );
	hErecoVsQ2Mode_NumuBarToNueBar[ifile]->Scale(1. / hErecoVsQ2Mode_NumuBarToNueBar[ifile]->Integral() );
      }
      if(draw_NumuToNue){
	hQ2Mode_NumuToNue[ifile]->Scale(1. / hQ2Mode_NumuToNue[ifile]->Integral() );
	hErecoVsQ2Mode_NumuToNue[ifile]->Scale(1. / hErecoVsQ2Mode_NumuToNue[ifile]->Integral() );
      }
      if(draw_NueBeam){
    	hQ2Mode_NueBeam[ifile]->Scale(1. / hQ2Mode_NueBeam[ifile]->Integral() );
    	hErecoVsQ2Mode_NueBeam[ifile]->Scale(1. / hErecoVsQ2Mode_NueBeam[ifile]->Integral() );
      }
      if(draw_NumuBeam){
  	hQ2Mode_NumuBeam[ifile]->Scale(1. / hQ2Mode_NumuBeam[ifile]->Integral() );
  	hErecoVsQ2Mode_NumuBeam[ifile]->Scale(1. / hErecoVsQ2Mode_NumuBeam[ifile]->Integral() );
      }
      if(draw_NumuBarBeam){
	hQ2Mode_NumuBarBeam[ifile]->Scale(1. / hQ2Mode_NumuBarBeam[ifile]->Integral() );
	hErecoVsQ2Mode_NumuBarBeam[ifile]->Scale(1. / hErecoVsQ2Mode_NumuBarBeam[ifile]->Integral() );
      }
      if(draw_NC){
	hQ2Mode_NC[ifile]->Scale(1. / hQ2Mode_NC[ifile]->Integral() );
	hErecoVsQ2Mode_NC[ifile]->Scale(1. / hErecoVsQ2Mode_NC[ifile]->Integral() );
      }
    }
    
    double max = -999;
    if(draw_NumuBarToNueBar)
      if(hQ2Mode_NumuBarToNueBar[ifile]->GetMaximum() > max) max = hQ2Mode_NumuBarToNueBar[ifile]->GetMaximum(); 
    if(draw_NumuToNue)    
      if(hQ2Mode_NumuToNue[ifile]->GetMaximum() > max)       max = hQ2Mode_NumuToNue[ifile]->GetMaximum(); 
    if(draw_NueBeam)
      if(hQ2Mode_NueBeam[ifile]->GetMaximum() > max)         max = hQ2Mode_NueBeam[ifile]->GetMaximum(); 
    if(draw_NumuBeam)
      if(hQ2Mode_NumuBeam[ifile]->GetMaximum() > max)        max = hQ2Mode_NumuBeam[ifile]->GetMaximum(); 
    if(draw_NumuBarBeam)
      if(hQ2Mode_NumuBarBeam[ifile]->GetMaximum() > max)     max = hQ2Mode_NumuBarBeam[ifile]->GetMaximum(); 
    if(draw_NC)
      if(hQ2Mode_NC[ifile]->GetMaximum() > max)              max = hQ2Mode_NC[ifile]->GetMaximum(); 
    
    hQ2Mode_init[ifile]->SetMaximum(max*1.2);

    hQ2Mode_NumuBarToNueBar[ifile]->SetFillColor(0);
    hQ2Mode_NumuToNue[ifile]->SetFillColor(0);
    hQ2Mode_NueBeam[ifile]->SetFillColor(0);
    hQ2Mode_NumuBeam[ifile]->SetFillColor(0);
    hQ2Mode_NumuBarBeam[ifile]->SetFillColor(0);
    hQ2Mode_NC[ifile]->SetFillColor(0);
    
    hQ2Mode_NumuBarBeam[ifile]->SetLineColor(kBlue);
    hQ2Mode_NumuBeam[ifile]->SetLineColor(kCyan);
    hQ2Mode_NueBeam[ifile]->SetLineColor(kYellow);
    hQ2Mode_NumuBarToNueBar[ifile]->SetLineColor(kGreen);
    hQ2Mode_NumuToNue[ifile]->SetLineColor(kGray+3);
    hQ2Mode_NC[ifile]->SetLineColor(kRed);
    
    hQ2Mode_NumuBarBeam[ifile]->SetLineWidth(2);
    hQ2Mode_NumuBeam[ifile]->SetLineWidth(2);
    hQ2Mode_NueBeam[ifile]->SetLineWidth(2);
    hQ2Mode_NumuBarToNueBar[ifile]->SetLineWidth(2);
    hQ2Mode_NumuToNue[ifile]->SetLineWidth(2);
    hQ2Mode_NC[ifile]->SetLineWidth(2);
    
    //sprintf(name,"cNuNubarQ2_file%i",ifile);
    sprintf(name,"cNuNubarQ2_%s",filename[ifile][1]);
    cNuNubarQ2[ifile] = new TCanvas(name,name);
    hQ2Mode_init[ifile]->DrawCopy("");
    if(draw_NumuToNue)
      hQ2Mode_NumuToNue[ifile]->DrawCopy("same");
    if(draw_NumuBarToNueBar)
      hQ2Mode_NumuBarToNueBar[ifile]->DrawCopy("same");
    if(draw_NueBeam)
      hQ2Mode_NueBeam[ifile]->DrawCopy("same");
    if(draw_NumuBeam)
      hQ2Mode_NumuBeam[ifile]->DrawCopy("same");
    if(draw_NumuBarBeam)
      hQ2Mode_NumuBarBeam[ifile]->DrawCopy("same");
    if(draw_NC)
      hQ2Mode_NC[ifile]->DrawCopy("same");
   
    legX[ifile] = new TLegend(0.6, 0.4, 0.9, 0.9); 
    legX[ifile]->SetFillColor(0);
    legX[ifile]->SetBorderSize(0);
    if(draw_NumuToNue)
      legX[ifile]->AddEntry(hQ2Mode_NumuToNue[ifile],"#nu_{#mu} #rightarrow #nu_{e}","l");
    if(draw_NumuBarToNueBar)
      legX[ifile]->AddEntry(hQ2Mode_NumuBarToNueBar[ifile],"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","l");
    if(draw_NC)
      legX[ifile]->AddEntry(hQ2Mode_NC[ifile],"NC","l");
    if(draw_NueBeam)
      legX[ifile]->AddEntry(hQ2Mode_NueBeam[ifile],"#nu_{e}/#bar{#nu}_{e} intrinsic","l");
    if(draw_NumuBeam)
      legX[ifile]->AddEntry(hQ2Mode_NumuBeam[ifile],"#nu_{#mu} intrinsic","l");
    if(draw_NumuBarBeam)
      legX[ifile]->AddEntry(hQ2Mode_NumuBarBeam[ifile],"#bar{#nu}_{#mu} intrinsic","l");
    legX[ifile]->Draw();


    // Calculate separability

    sprintf(name,"hSeparation_%s",filename[ifile][1]);
    hSeparation_Q2[ifile] = (TH1D*)hQ2Mode_NumuToNue[ifile]->Clone(name);
    hSeparation_Q2[ifile]->Reset();
    for(int ibin=0;ibin<hSeparation_Q2[ifile]->GetNbinsX();ibin++){
      double separ = 
    	hQ2Mode_NumuToNue[ifile]->GetBinContent(ibin+1) - hQ2Mode_NumuBarToNueBar[ifile]->GetBinContent(ibin+1);      
      hSeparation_Q2[ifile]->SetBinContent(ibin+1,separ);

      separation[ifile] += TMath::Abs(separ);
    }

    sprintf(name,"cSeparation_1D_%s",filename[ifile][1]);
    cSeparation_Q2[ifile] = new TCanvas(name,name);
    hSeparation_Q2[ifile]->Draw();
    hSeparation_Q2[ifile]->SetLineColor(ifile+1);
    

    // 2D templates

    sprintf(name,"cErecoVsQ2_NumuBarToNueBar_%s",filename[ifile][1]);
    cErecoVsQ2_NumuBarToNueBar[ifile] = new TCanvas(name,name); 
    hErecoVsQ2Mode_NumuBarToNueBar[ifile]->Draw("colz");

    sprintf(name,"cErecoVsQ2_NumuToNue_%s",filename[ifile][1]);
    cErecoVsQ2_NumuToNue[ifile] = new TCanvas(name,name); 
    hErecoVsQ2Mode_NumuToNue[ifile]->Draw("colz");

    
    // Calculate separability
    
    sprintf(name,"hSeparation_2D_%s",filename[ifile][1]);
    hSeparation_ErecoVsQ2[ifile] = (TH2D*)hErecoVsQ2Mode_NumuToNue[ifile]->Clone(name);
    hSeparation_ErecoVsQ2[ifile]->Reset();
    
    for(int ibinX=0;ibinX<hSeparation_ErecoVsQ2[ifile]->GetNbinsX();ibinX++){
      for(int ibinY=0;ibinY<hSeparation_ErecoVsQ2[ifile]->GetNbinsY();ibinY++){

	double separ = 
	  hErecoVsQ2Mode_NumuToNue[ifile]->GetBinContent(ibinX+1,ibinY+1) - hErecoVsQ2Mode_NumuBarToNueBar[ifile]->GetBinContent(ibinX+1,ibinY+1);      

	if( hSeparation_ErecoVsQ2[ifile]->GetYaxis()->GetBinCenter(ibinY+1) > energy_up ) separ = 0.;

	hSeparation_ErecoVsQ2[ifile]->SetBinContent(ibinX+1,ibinY+1,separ);
	separation_2D[ifile] += TMath::Abs(separ);

	double ntot_bin = hErecoVsQ2Mode_NumuToNue[ifile]->GetBinContent(ibinX+1,ibinY+1) + hErecoVsQ2Mode_NumuBarToNueBar[ifile]->GetBinContent(ibinX+1,ibinY+1);  

	separation_weighted_2D[ifile] += TMath::Abs(separ) * ntot_bin;
      }
    }
    double ntot = hErecoVsQ2Mode_NumuToNue[ifile]->Integral() + hErecoVsQ2Mode_NumuBarToNueBar[ifile]->Integral();  
    separation_weighted_2D[ifile] /= ntot;
    
    sprintf(name,"cSeparation_2D_%s",filename[ifile][1]);
    cSeparation_Q2[ifile] = new TCanvas(name,name);
    hSeparation_ErecoVsQ2[ifile]->Draw("colz");
    
  } // end loop over files
  
  cout << endl;
  cout << "Check energy range: " << energy_down << " - " << energy_up << endl;
  cout << endl;
  cout << "Total separation from 1D templates:" << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << " - " << filename[ifile][1] << " --> " << separation[ifile] << endl;
  }
  cout << endl;
  cout << "Total separation from 2D templates:" << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << " - " << filename[ifile][1] << " --> " << separation_2D[ifile] << endl;
  }
  cout << "Total weighted separation from 2D templates:" << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << " - " << filename[ifile][1] << " --> " << separation_weighted_2D[ifile] << endl;
  }

  
}
