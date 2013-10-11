/*
// -----------------------------------------------------------------------------------------------------------
This macro writes to text files the bin numbers of TH2D histograms in bands around the confidence contours. 
These bin numbers are written in pairs of (sinsq2theta, dmsq), and each pair is on a new line. 

N.Grant
// -----------------------------------------------------------------------------------------------------------
*/

void write_contour_bins_to_textfile
(
  /* Input ROOT file with confidence regions in TH2D format */
  const char * infile,
  
  /* The number of bins included on both sides of the confidence contours. 
     Bin numbers are written to the edges of the input histograms but not beyond them. */
  const int bandWidth = 1,

  /* If write_single_list is true, a single list of bins is written without duplication for all the "write_cl" that are true.
     However this list is written only if the numbers of bins in the different TH2Ds are equal.
     If write_single_list is false, a separate list is written for each cl for which "write_cl" is true; there is no
     duplication in any single list, but there may be duplication between different lists. 
  */
  bool write_single_list = true,

  bool write_contours_with_systematics=true,
  bool write_contours_without_systematics=true,
  bool write_68cl=true,
  bool write_90cl=true,
  bool write_95cl=false,
  bool write_99cl=false,

  /* Names of TH2D histograms (confidence regions) to allow retrieval from file 
     The default names are the ones produced by the plot_feldman_cousins_contour.C macro (and its variatons).
     But the option to specify names is provided here since the contours might be produced from other code too
     (eg non-FC contours, etc). */
  const char * hstname_cl68_statonly = "hCL68StatOnly", 
  const char * hstname_cl90_statonly = "hCL90StatOnly", 
  const char * hstname_cl95_statonly = "hCL95StatOnly", 
  const char * hstname_cl99_statonly = "hCL99StatOnly",
  const char * hstname_cl68_statsyst = "hCL68StatSyst", 
  const char * hstname_cl90_statsyst = "hCL90StatSyst", 
  const char * hstname_cl95_statsyst = "hCL95StatSyst",
  const char * hstname_cl99_statsyst = "hCL99StatSyst",

  //set this if you are inputting a single dchisq surface, that has not had dchisqcrit subtracted
  //(it will subtract the conventional dchisqcrit, e.g. 2.30 for 68%)
  bool constdchisqsurface = true,

  //convert output of a limited range histogram into one that sits on the nominal fc surface
  bool convert_to_fc_surface_range = true
)
{
  cout << "Writing bin numbers in a band " << bandWidth
       << ((bandWidth == 1) ? " bin" : " bins")
       << " wide around confidence contours from " << infile << endl;
 
  // Read input contours
  TFile *contourfile = TFile::Open(infile, "READ");
  contourfile->cd();

  TH2D *hContStatOnlyF[4] = {0};
  hContStatOnlyF[0] = (TH2D*) contourfile->Get(hstname_cl68_statonly);
  hContStatOnlyF[1] = (TH2D*) contourfile->Get(hstname_cl90_statonly);
  hContStatOnlyF[2] = (TH2D*) contourfile->Get(hstname_cl95_statonly);
  hContStatOnlyF[3] = (TH2D*) contourfile->Get(hstname_cl99_statonly);
  //clone histograms (for case where the inputs are the same histogram)
  TH2D *hContStatOnly[4] = {0};
  hContStatOnly[0] = (TH2D*)hContStatOnlyF[0]->Clone("hContStatOnly0");
  hContStatOnly[1] = (TH2D*)hContStatOnlyF[1]->Clone("hContStatOnly1");
  hContStatOnly[2] = (TH2D*)hContStatOnlyF[2]->Clone("hContStatOnly2");
  hContStatOnly[3] = (TH2D*)hContStatOnlyF[3]->Clone("hContStatOnly3");

  TH2D *hContStatSystF[4] = {0};
  hContStatSystF[0] = (TH2D*) contourfile->Get(hstname_cl68_statsyst);
  hContStatSystF[1] = (TH2D*) contourfile->Get(hstname_cl90_statsyst);
  hContStatSystF[2] = (TH2D*) contourfile->Get(hstname_cl95_statsyst);
  hContStatSystF[3] = (TH2D*) contourfile->Get(hstname_cl99_statsyst);
  //clone histograms (for case where the inputs are the same histogram)
  TH2D *hContStatSyst[4] = {0};
  hContStatSyst[0] = (TH2D*)hContStatSystF[0]->Clone("hContStatSyst0");
  hContStatSyst[1] = (TH2D*)hContStatSystF[1]->Clone("hContStatSyst1");
  hContStatSyst[2] = (TH2D*)hContStatSystF[2]->Clone("hContStatSyst2");
  hContStatSyst[3] = (TH2D*)hContStatSystF[3]->Clone("hContStatSyst3");

  // The 'default' FC surface binning
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.5;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int)(((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)(((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);

  TString outfilename;

  bool firstXBinFound;
  int yZeroBin;
  int yMaxBin;

  bool equalBinsStatOnly = true;
  for(int i=0; i<4; i++)
    {
    if(!hContStatOnly[i]) continue;
    for(int j=0; j<4; j++)
      {
      if(hContStatOnly[i]->GetNbinsX() != hContStatOnly[j]->GetNbinsX() || hContStatOnly[i]->GetNbinsY() != hContStatOnly[j]->GetNbinsY())
	equalBinsStatOnly = false;
      }
    }  

  bool equalBinsStatSyst = true;
  for(int i=0; i<4; i++)
    {
    if(!hContStatSyst[i]) continue;
    for(int j=0; j<4; j++)
      {
      if(hContStatSyst[i]->GetNbinsX() != hContStatSyst[j]->GetNbinsX() || hContStatSyst[i]->GetNbinsY() != hContStatSyst[j]->GetNbinsY())
	equalBinsStatSyst = false;
      }
    }

  if(constdchisqsurface) {
    double constdchisq[4] = {2.30, 4.61, 5.99, 9.21};
    cout << "Converting constant dchisq surface into const dchisqcrit-subtracted surface" << endl;

    //stat only
    for(int is=0; is<4; is++) {
      if(!hContStatOnly[is]) continue;
      for(int ix = 1; ix <= hContStatOnly[is]->GetNbinsX(); ix++)
	for(int iy = 1; iy <= hContStatOnly[is]->GetNbinsY(); iy++)
	  hContStatOnly[is]->SetBinContent(ix,iy, hContStatOnly[is]->GetBinContent(ix,iy) - constdchisq[is]);
    }

    //stat+syst
    for(int is=0; is<4; is++) {
      if(!hContStatSyst[is]) continue;
      for(int ix = 1; ix <= hContStatSyst[is]->GetNbinsX(); ix++)
	for(int iy = 1; iy <= hContStatSyst[is]->GetNbinsY(); iy++)
	  hContStatSyst[is]->SetBinContent(ix,iy, hContStatSyst[is]->GetBinContent(ix,iy) - constdchisq[is]);
    }
  }//constdchisqsurface


  //if requested, convert the input dchisq surface histogram range to the 'FC' binning scheme
  int nbinx_addon_statonly = 0, nbiny_addon_statonly = 0;
  if(equalBinsStatOnly) {
    cout << hContStatOnly[0]->GetNbinsX() << " bins in x, from " << hContStatOnly[0]->GetXaxis()->GetXmin() << " to " << hContStatOnly[0]->GetXaxis()->GetXmax() << endl;
    cout << hContStatOnly[0]->GetNbinsY() << " bins in y, from " << hContStatOnly[0]->GetYaxis()->GetXmin() << " to " << hContStatOnly[0]->GetYaxis()->GetXmax() << endl;
    if(convert_to_fc_surface_range) {
      cout << mixing23_nbins << " bins in x(fc), from " << mixing23_min << " to " << mixing23_max << endl
	   << dm23sq_nbins   << " bins in y(fc), from " << dm23sq_min   << " to " << dm23sq_max   << endl;

      //find conversion for x
      double xstep = (hContStatOnly[0]->GetXaxis()->GetXmax() - hContStatOnly[0]->GetXaxis()->GetXmin()) / (double)hContStatOnly[0]->GetNbinsX();
      double xnbins = hContStatOnly[0]->GetNbinsX();
      if(TMath::Abs(xstep - mixing23_step) > 1E-8)
	cout << "Cannot convert xaxis to fc axis. Different step size" << endl;
      else {
	nbinx_addon_statonly = (int)(((hContStatOnly[0]->GetXaxis()->GetXmin() - mixing23_min) / mixing23_step) + 0.01);
	cout << "Adding " << nbinx_addon_statonly << " to the x axis bins" << endl;
      }

      //find conversion for y
      double ystep = (hContStatOnly[0]->GetYaxis()->GetXmax() - hContStatOnly[0]->GetYaxis()->GetXmin()) / (double)hContStatOnly[0]->GetNbinsY();
      double ynbins = hContStatOnly[0]->GetNbinsY();
      if(TMath::Abs(ystep - dm23sq_step) > 1E-8)
	cout << "Cannot convert yaxis to fc axis. Different step size" << endl;
      else {
	nbiny_addon_statonly = (int)(((hContStatOnly[0]->GetYaxis()->GetXmin() - dm23sq_min) / dm23sq_step) + 0.01);
	cout << "Adding " << nbiny_addon_statonly << " to the y axis bins" << endl;
      }
    }
  }

  //setup the visualisation histograms (with appropriate binning)
  TH2D * visualise_output_statonly = 0;
  TH2D * visualise_output_statsyst = 0;
  TRandom3 rand(0);
  if(convert_to_fc_surface_range) {
    visualise_output_statonly = new TH2D("","",mixing23_nbins+1,mixing23_min,mixing23_max+mixing23_step,dm23sq_nbins+1,dm23sq_min,dm23sq_max+mixing23_step);
    visualise_output_statsyst = new TH2D("","",mixing23_nbins+1,mixing23_min,mixing23_max+mixing23_step,dm23sq_nbins+1,dm23sq_min,dm23sq_max+mixing23_step);
  }
  else {
    visualise_output_statonly = new TH2D("","",hContStatOnly[0]->GetNbinsX(),hContStatOnly[0]->GetXaxis()->GetXmin(),hContStatOnly[0]->GetXaxis()->GetXmax(),hContStatOnly[0]->GetNbinsY(),hContStatOnly[0]->GetYaxis()->GetXmin(),hContStatOnly[0]->GetYaxis()->GetXmax());
    visualise_output_statsyst = new TH2D("","",hContStatSyst[0]->GetNbinsX(),hContStatSyst[0]->GetXaxis()->GetXmin(),hContStatSyst[0]->GetXaxis()->GetXmax(),hContStatSyst[0]->GetNbinsY(),hContStatSyst[0]->GetYaxis()->GetXmin(),hContStatSyst[0]->GetYaxis()->GetXmax());
  }


  //if requested, convert the input dchisq surface histogram range to the 'FC' binning scheme
  int nbinx_addon_statsyst = 0, nbiny_addon_statsyst = 0;
  if(equalBinsStatSyst) {
    cout << hContStatSyst[0]->GetNbinsX() << " bins in x, from " << hContStatSyst[0]->GetXaxis()->GetXmin() << " to " << hContStatSyst[0]->GetXaxis()->GetXmax() << endl;
    cout << hContStatSyst[0]->GetNbinsY() << " bins in y, from " << hContStatSyst[0]->GetYaxis()->GetXmin() << " to " << hContStatSyst[0]->GetYaxis()->GetXmax() << endl;
    if(convert_to_fc_surface_range) {
      cout << mixing23_nbins << " bins in x(fc), from " << mixing23_min << " to " << mixing23_max << endl
	   << dm23sq_nbins   << " bins in y(fc), from " << dm23sq_min   << " to " << dm23sq_max   << endl;

      //find conversion for x
      double xstep = (hContStatSyst[0]->GetXaxis()->GetXmax() - hContStatSyst[0]->GetXaxis()->GetXmin()) / (double)hContStatSyst[0]->GetNbinsX();
      double xnbins = hContStatSyst[0]->GetNbinsX();
      if(TMath::Abs(xstep - mixing23_step) > 1E-8)
	cout << "Cannot convert xaxis to fc axis. Different step size" << endl;
      else {
	nbinx_addon_statsyst = (int)(((hContStatSyst[0]->GetXaxis()->GetXmin() - mixing23_min) / mixing23_step) + 0.01);
	cout << "Adding " << nbinx_addon_statsyst << " to the x axis bins" << endl;
      }

      //find conversion for y
      double ystep = (hContStatSyst[0]->GetYaxis()->GetXmax() - hContStatSyst[0]->GetYaxis()->GetXmin()) / (double)hContStatSyst[0]->GetNbinsY();
      double ynbins = hContStatSyst[0]->GetNbinsY();
      if(TMath::Abs(ystep - dm23sq_step) > 1E-8)
	cout << "Cannot convert yaxis to fc axis. Different step size" << endl;
      else {
	nbiny_addon_statsyst = (int)(((hContStatSyst[0]->GetYaxis()->GetXmin() - dm23sq_min) / dm23sq_step) + 0.01);
	cout << "Adding " << nbiny_addon_statsyst << " to the y axis bins" << endl;
      }
    }
  }

  if(write_contours_without_systematics)
    {
    if(write_single_list && !equalBinsStatOnly)
      cout<<"Not writing stat only confidence contour bin numbers to a single list since binnings in different TH2Ds are not equal."<<endl;
    else 
      {
      cout << "Writing stat only bin numbers for ";
      if(write_68cl) cout<<"68% "; 
      if(write_90cl) cout<<"90% "; 
      if(write_95cl) cout<<"95% "; 
      if(write_99cl) cout<<"99% "; 
      cout << "confidence contours ";
      if(write_single_list) cout << "to a single list";
      else                  cout << "to separate lists";
      cout << endl;
      }
    
    outfilename = "ContourBinsStatOnly";
    if(write_68cl) outfilename += "_68";
    if(write_90cl) outfilename += "_90";
    if(write_95cl) outfilename += "_95";
    if(write_99cl) outfilename += "_99";
    outfilename += ".data";

    if(!(write_single_list && !equalBinsStatOnly)) 
      ofstream StatOnlyBins(outfilename);

    if(write_single_list)
      {
      if(equalBinsStatOnly)
	{
        firstXBinFound = false;
        for(int i=1; i<=hContStatOnly[0]->GetNbinsX(); i++)
          {
          yMaxBin = 0;
          //start at bin 2 since bin j is compared with bin j-1                                                                                                                                   
          for(int j=2; j<=hContStatOnly[0]->GetNbinsY(); j++)
            {
	    for(int h=0; h<4; h++)
	      {
	      if((h==0 && !write_68cl) || (h==1 && !write_90cl) || (h==2 && !write_95cl) || (h==3 && !write_99cl))
	        continue;
              if((hContStatOnly[h]->GetBinContent(i, j) > 0 && hContStatOnly[h]->GetBinContent(i, j-1) < 0)
                 || (hContStatOnly[h]->GetBinContent(i, j) < 0 && hContStatOnly[h]->GetBinContent(i, j-1) > 0))
                {
                if(fabs(hContStatOnly[h]->GetBinContent(i, j)) <= fabs(hContStatOnly[h]->GetBinContent(i, j-1)))
                  yZeroBin = j;
                else
                  yZeroBin = j-1;
                //this is to write the bin numbers that are lower in x than the countours                                                                                                    
                if(!firstXBinFound)
                  {
                  for(int k=i-bandWidth; k<i; k++)
                    {
                    for(int l=yZeroBin-bandWidth; l<=yZeroBin+bandWidth; l++)
                      {
			if(k >=1 && k <= hContStatOnly[h]->GetNbinsX() && l >=1 && l <= hContStatOnly[h]->GetNbinsY()) {
			  StatOnlyBins << k+nbinx_addon_statonly << "\t" << l+nbiny_addon_statonly << "\n";
			  visualise_output_statonly->SetBinContent(k+nbinx_addon_statonly,l+nbiny_addon_statonly,rand.Uniform(10,100));
			}
                      }
                    }
                  firstXBinFound = true;
                  }
                for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
                  {
                  //the if(m > yMaxBin) condition is to avoid duplicating the bin numbers                                                                                                   
		    if(i >= 1 && i <= hContStatOnly[h]->GetNbinsX() && m >= 1 && m <= hContStatOnly[h]->GetNbinsY() && m > yMaxBin) {
		      StatOnlyBins << i+nbinx_addon_statonly << "\t" << m+nbiny_addon_statonly <<"\n";
		      visualise_output_statonly->SetBinContent(i+nbinx_addon_statonly,m+nbiny_addon_statonly,rand.Uniform(10,100));
		    }
                 }
                yMaxBin = yZeroBin+bandWidth;
                }//end of if((hContStatOnly[h]->GetBinContent(i, j) > 0 && hContStatOnly[h]->GetBinContent(i, j-1) < 0).....                                                                     
	      }//end of for(int h=0; h<4; h++)
            }//end of for(int j=1; j<=hContStatOnly[h]->GetNbinsY(); j++)                                                                                                                        
          }//end of for(int i=1; i<=hContStatOnly[h]->GetNbinsX(); i++)  
	}//end of if(equalBinsStatOnly)
      }//end of if(write_single_list)
      else
	{
        for(int h=0; h<4; h++)
	  {
	  if((h==0 && !write_68cl) || (h==1 && !write_90cl) || (h==2 && !write_95cl) || (h==3 && !write_99cl))
	    continue;
	  firstXBinFound = false;
	  for(int i=1; i<=hContStatOnly[h]->GetNbinsX(); i++)
	    {
	    yMaxBin = 0;
            //start at bin 2 since bin j is compared with bin j-1 
	    for(int j=2; j<=hContStatOnly[h]->GetNbinsY(); j++)
	      {
	      if((hContStatOnly[h]->GetBinContent(i, j) > 0 && hContStatOnly[h]->GetBinContent(i, j-1) < 0)
		  || (hContStatOnly[h]->GetBinContent(i, j) < 0 && hContStatOnly[h]->GetBinContent(i, j-1) > 0))
	        {
		if(fabs(hContStatOnly[h]->GetBinContent(i, j)) <= fabs(hContStatOnly[h]->GetBinContent(i, j-1)))
	          yZeroBin = j;
		else
		  yZeroBin = j-1;
		//this is to write the bin numbers that are lower in x than the countours 
	        if(!firstXBinFound)
		  {
		  for(int k=i-bandWidth; k<i; k++)
		    {
		    for(int l=yZeroBin-bandWidth; l<=yZeroBin+bandWidth; l++)
		      {
			if(k >=1 && k <= hContStatOnly[h]->GetNbinsX() && l >=1 && l <= hContStatOnly[h]->GetNbinsY()) {
			  StatOnlyBins << k+nbinx_addon_statonly << "\t" << l+nbiny_addon_statonly << "\n";
			  visualise_output_statonly->SetBinContent(k+nbinx_addon_statonly,l+nbiny_addon_statonly,rand.Uniform(10,100));
			}
		      }
		    }
		  firstXBinFound = true;
		  }
		for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
		  {
		  //the if(m > yMaxBin) condition is to avoid duplicating the bin numbers 
		    if(i >= 1 && i <= hContStatOnly[h]->GetNbinsX() && m >= 1 && m <= hContStatOnly[h]->GetNbinsY() && m > yMaxBin) {
		      StatOnlyBins << i+nbinx_addon_statonly << "\t" << m+nbiny_addon_statonly <<"\n";
		      visualise_output_statonly->SetBinContent(i+nbinx_addon_statonly,m+nbiny_addon_statonly,rand.Uniform(10,100));
		    }
		  }
		  yMaxBin = yZeroBin+bandWidth;
		}//end of if((hContStatOnly[h]->GetBinContent(i, j) > 0 && hContStatOnly[h]->GetBinContent(i, j-1) < 0).....                                                                       
	      }//end of for(int j=1; j<=hContStatOnly[h]->GetNbinsY(); j++)                                                                                                                          
	    }//end of for(int i=1; i<=hContStatOnly[h]->GetNbinsX(); i++)                                                                                                                             

	  StatOnlyBins << "\n";

	  }//end of for(int h=0; h<4; h++)                                                                                                                                                            
        }//end of else to if(write_single_list)
      }//end of if(write_contours_without_systematics) 


  if(write_contours_with_systematics)
    {
    if(write_single_list && !equalBinsStatSyst)
      cout<<"Not writing stat+syst confidence contour bin numbers to a single list since binnings in different TH2Ds are not equal."<<endl;
    else
      { 
      cout << "Writing stat+syst bin numbers for ";
      if(write_68cl) cout << "68% ";
      if(write_90cl) cout << "90% ";
      if(write_95cl) cout << "95% ";
      if(write_99cl) cout << "99% ";
      cout << "confidence contours ";   
      if(write_single_list) cout << "to a single list";
      else                  cout << "to separate lists";
      cout << endl;
      }

    outfilename = "ContourBinsStatSyst";
    if(write_68cl) outfilename += "_68";
    if(write_90cl) outfilename += "_90";
    if(write_95cl) outfilename += "_95";
    if(write_99cl) outfilename += "_99";
    outfilename += ".data";

    if(!(write_single_list &&!equalBinsStatSyst))
      ofstream StatSystBins(outfilename);

    if(write_single_list)
      {
      if(equalBinsStatSyst)
	{
	firstXBinFound = false;
	for(int i=1; i<=hContStatSyst[0]->GetNbinsX(); i++)
	  {
	  yMaxBin = 0;
	  //start at bin 2 since bin j is compared with bin j-1                                                                                                                                   
          for(int j=2; j<=hContStatSyst[0]->GetNbinsY(); j++)
	    {
	    for(int h=0; h<4; h++)
	      {
	      if((h==0 && !write_68cl) || (h==1 && !write_90cl) || (h==2 && !write_95cl) ||(h==3 && !write_99cl))
		continue;
	      if((hContStatSyst[h]->GetBinContent(i, j) > 0 && hContStatSyst[h]->GetBinContent(i, j-1) < 0)
		 || (hContStatSyst[h]->GetBinContent(i, j) < 0 && hContStatSyst[h]->GetBinContent(i, j-1) > 0))
		{
		if(fabs(hContStatSyst[h]->GetBinContent(i, j)) <= fabs(hContStatSyst[h]->GetBinContent(i, j-1)))
		  yZeroBin = j;
		else
		  yZeroBin = j-1;
	      if(!firstXBinFound)
		{
		for(int k=i-bandWidth; k<i; k++)
		  {
		  for(int l=yZeroBin-bandWidth; l<=yZeroBin+bandWidth; l++)
		    {
		    if(k >=1 && k <= hContStatSyst[h]->GetNbinsX() && l >= 1 && l <= hContStatSyst[h]->GetNbinsY())
		      {
			StatSystBins << k+nbinx_addon_statsyst << "\t" << l+nbiny_addon_statsyst << "\n";
			visualise_output_statsyst->SetBinContent(k+nbinx_addon_statsyst,l+nbiny_addon_statsyst,rand.Uniform(10,100));
		      }
		    }
		  }
		  firstXBinFound = true;
		}
                for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
		  {
		  //the if(m > yMaxBin) condition is to avoid duplicating the bin numbers
		  if(i >= 1 && i <= hContStatSyst[h]->GetNbinsX() && m >= 1 && m <= hContStatSyst[h]->GetNbinsY() && m > yMaxBin)
		    {
		      StatSystBins << i+nbinx_addon_statsyst << "\t" << m+nbiny_addon_statsyst <<"\n";
		      visualise_output_statsyst->SetBinContent(i+nbinx_addon_statsyst,m+nbiny_addon_statsyst,rand.Uniform(10,100));
		    }
		  }
		    yMaxBin = yZeroBin+bandWidth;
		}//end of if((hContStatSyst[h]->GetBinContent(i, j) > 0 && hContStatSyst[h]->GetBinContent(i, j-1) < 0).....                                                                        
	      }//end of for(int h=0; h<4; h++)
            }//end of for(int j=1; j<=hContStatSyst[h]->GetNbinsY(); j++)                                                                                                                           
          }//end of for(int i=1; i<=hContStatSyst[h]->GetNbinsX(); i++)                                                                                                                                
	}//end of if(equalBinsStatSyst)
      }//end of if(write_single_list)
    else
      {
      for(int h=0; h<4; h++)
        {
        if((h==0 && !write_68cl) || (h==1 && !write_90cl) || (h==2 && !write_95cl) ||(h==3 && !write_99cl))
	  continue;
        firstXBinFound = false;
        for(int i=1; i<=hContStatSyst[h]->GetNbinsX(); i++)
          {
          yMaxBin = 0;
          //start at bin 2 since bin j is compared with bin j-1 
          for(int j=2; j<=hContStatSyst[h]->GetNbinsY(); j++)
            {
            if((hContStatSyst[h]->GetBinContent(i, j) > 0 && hContStatSyst[h]->GetBinContent(i, j-1) < 0)
               || (hContStatSyst[h]->GetBinContent(i, j) < 0 && hContStatSyst[h]->GetBinContent(i, j-1) > 0))
	      {
	      if(fabs(hContStatSyst[h]->GetBinContent(i, j)) <= fabs(hContStatSyst[h]->GetBinContent(i, j-1)))
                yZeroBin = j;
	      else
	        yZeroBin = j-1;
	      if(!firstXBinFound)
	        {
	        for(int k=i-bandWidth; k<i; k++)
		  {
		  for(int l=yZeroBin-bandWidth; l<=yZeroBin+bandWidth; l++)
		    {
		    if(k >=1 && k <= hContStatSyst[h]->GetNbinsX() && l>=1 && l <= hContStatSyst[h]->GetNbinsY())
		      {
			StatSystBins << k+nbinx_addon_statsyst << "\t" << l+nbiny_addon_statsyst << "\n";  
			visualise_output_statsyst->SetBinContent(k+nbinx_addon_statsyst,l+nbiny_addon_statsyst,rand.Uniform(10,100));
		      }
		    }
	       	  }
                firstXBinFound = true;
	        }
	      for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
	        {
		//the if(m > yMaxBin) statement is to avoid duplicating the bin numbers
		if(i >= 1 && i <= hContStatSyst[h]->GetNbinsX() && m >= 1 && m <= hContStatSyst[h]->GetNbinsY() && m > yMaxBin)
		  {
		    StatSystBins << i+nbinx_addon_statsyst << "\t" << m+nbiny_addon_statsyst <<"\n";
		    visualise_output_statsyst->SetBinContent(i+nbinx_addon_statsyst,m+nbiny_addon_statsyst,rand.Uniform(10,100));
		  }
	        }
              yMaxBin = yZeroBin+bandWidth;	            
	    }//end of if((hContStatSyst[h]->GetBinContent(i, j) > 0 && hContStatSyst[h]->GetBinContent(i, j-1) < 0).....
          }//end of for(int j=1; j<=hContStatSyst[h]->GetNbinsY(); j++)
        }//end of for(int i=1; i<=hContStatSyst[h]->GetNbinsX(); i++)
  
      StatSystBins << "\n";

      }//end of for(int h=0; h<4; h++)
    }//end of else to if(write_single_list)
  }//end of if(write_contours_with_systematics)
  
  TCanvas * c = new TCanvas();
  if(write_contours_without_systematics) {
    visualise_output_statonly->SetStats(0);
    visualise_output_statonly->Draw("colz");
    for(int i = 0; i < 4; i++) {
      if((i==0 && !write_68cl) || (i==1 && !write_90cl) || (i==2 && !write_95cl) ||(i==3 && !write_99cl))
	continue;
      hContStatOnly[i]->SetContour(1);
      hContStatOnly[i]->SetContourLevel(0,0);
      hContStatOnly[i]->Draw("cont3same");
    }
    c->SaveAs("visual_statonly.pdf");
  }
  if(write_contours_with_systematics) {
    visualise_output_statsyst->SetStats(0);
    visualise_output_statsyst->Draw("colz");
    for(int i = 0; i < 4; i++) {
      if((i==0 && !write_68cl) || (i==1 && !write_90cl) || (i==2 && !write_95cl) ||(i==3 && !write_99cl))
	continue;
      hContStatSyst[i]->SetContour(1);
      hContStatSyst[i]->SetContourLevel(0,0);
      hContStatSyst[i]->Draw("cont3same");
    }
    c->SaveAs("visual_statsyst.pdf");
  }
}
