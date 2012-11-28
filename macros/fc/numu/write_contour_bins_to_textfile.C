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
  const char * hstname_cl99_statsyst = "hCL99StatSyst"
)
{
  cout << "Writing bin numbers in a band " << bandWidth
       << ((bandWidth == 1) ? " bin" : " bins")
       << " wide around confidence contours from " << infile << endl;
 
  // Read input contours
  TFile *contourfile = TFile::Open(infile, "READ");
  contourfile->cd();

  TH2D *hContStatOnly[4];
  hContStatOnly[0] = (TH2D*) contourfile->Get(hstname_cl68_statonly);
  hContStatOnly[1] = (TH2D*) contourfile->Get(hstname_cl90_statonly);
  hContStatOnly[2] = (TH2D*) contourfile->Get(hstname_cl95_statonly);
  hContStatOnly[3] = (TH2D*) contourfile->Get(hstname_cl99_statonly);

  TH2D *hContStatSyst[4];
  hContStatSyst[0] = (TH2D*) contourfile->Get(hstname_cl68_statsyst);
  hContStatSyst[1] = (TH2D*) contourfile->Get(hstname_cl90_statsyst);
  hContStatSyst[2] = (TH2D*) contourfile->Get(hstname_cl95_statsyst);
  hContStatSyst[3] = (TH2D*) contourfile->Get(hstname_cl99_statsyst);

  TString outfilename;

  bool firstXBinFound;
  int yZeroBin;
  int yMaxBin;

  bool equalBinsStatOnly = true;
  for(int i=0; i<4; i++)
    {
    for(int j=0; j<4; j++)
      {
      if(hContStatOnly[i]->GetNbinsX() != hContStatOnly[j]->GetNbinsX() || hContStatOnly[i]->GetNbinsY() != hContStatOnly[j]->GetNbinsY())
	equalBinsStatOnly = false;
      }
    }  

  bool equalBinsStatSyst = true;
  for(int i=0; i<4; i++)
    {
    for(int j=0; j<4; j++)
      {
      if(hContStatSyst[i]->GetNbinsX() != hContStatSyst[j]->GetNbinsX() || hContStatSyst[i]->GetNbinsY() != hContStatSyst[j]->GetNbinsY())
	equalBinsStatSyst = false;
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
		      if(k >=1 && k <= hContStatOnly[h]->GetNbinsX() && l >=1 && l <= hContStatOnly[h]->GetNbinsY())
		        StatOnlyBins << k << "\t" << l << "\n";
                      }
                    }
                  firstXBinFound = true;
                  }
                for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
                  {
                  //the if(m > yMaxBin) condition is to avoid duplicating the bin numbers                                                                                                   
		  if(i >= 1 && i <= hContStatOnly[h]->GetNbinsX() && m >= 1 && m <= hContStatOnly[h]->GetNbinsY() && m > yMaxBin)
                    StatOnlyBins << i << "\t" << m <<"\n";
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
		      if(k >=1 && k <= hContStatOnly[h]->GetNbinsX() && l >=1 && l <= hContStatOnly[h]->GetNbinsY())
		        StatOnlyBins << k << "\t" << l << "\n";
		      }
		    }
		  firstXBinFound = true;
		  }
		for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
		  {
		  //the if(m > yMaxBin) condition is to avoid duplicating the bin numbers 
		  if(i >= 1 && i <= hContStatOnly[h]->GetNbinsX() && m >= 1 && m <= hContStatOnly[h]->GetNbinsY() && m > yMaxBin)
		    StatOnlyBins << i << "\t" << m <<"\n";
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
		      StatSystBins << k << "\t" << l << "\n";
		    }
		  }
		  firstXBinFound = true;
		}
                for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
		  {
		  //the if(m > yMaxBin) condition is to avoid duplicating the bin numbers
		  if(i >= 1 && i <= hContStatSyst[h]->GetNbinsX() && m >= 1 && m <= hContStatSyst[h]->GetNbinsY() && m > yMaxBin)
		    StatSystBins << i << "\t" << m <<"\n";
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
		      StatSystBins << k << "\t" << l << "\n";  
		    }
	       	  }
                firstXBinFound = true;
	        }
	      for(int m=yZeroBin-bandWidth; m<=yZeroBin+bandWidth; m++)
	        {
		//the if(m > yMaxBin) statement is to avoid duplicating the bin numbers
		if(i >= 1 && i <= hContStatSyst[h]->GetNbinsX() && m >= 1 && m <= hContStatSyst[h]->GetNbinsY() && m > yMaxBin)
                  StatSystBins << i << "\t" << m <<"\n";
	        }
              yMaxBin = yZeroBin+bandWidth;	            
	    }//end of if((hContStatSyst[h]->GetBinContent(i, j) > 0 && hContStatSyst[h]->GetBinContent(i, j-1) < 0).....
          }//end of for(int j=1; j<=hContStatSyst[h]->GetNbinsY(); j++)
        }//end of for(int i=1; i<=hContStatSyst[h]->GetNbinsX(); i++)
  
      StatSystBins << "\n";

      }//end of for(int h=0; h<4; h++)
    }//end of else to if(write_single_list)
  }//end of if(write_contours_with_systematics)
  
}
