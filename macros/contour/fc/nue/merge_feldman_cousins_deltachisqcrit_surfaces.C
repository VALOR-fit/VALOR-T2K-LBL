//
// --------------------------------------------------------------------------------------------------------
// Merge the output of make_feldman_cousins_deltachisqcrit_surface.C
// in case the task was parallelized to many jobs submitted to a batch farm for faster results.
// Each such job would have computed \Delta\chi^{2} for part of the oscillation parameter space.
//
// Inputs
//  - nruns      : number of runs submitted by
//  - inpname    : part of the input file name (<inpname><job_number>.root)
//  - outputfile : name of output ROOT file with merged \Delta\chi^{2} surfaces
//              
// C.Andreopoulos
// --------------------------------------------------------------------------------------------------------
//

void merge_feldman_cousins_deltachisqcrit_surfaces(
  const int nruns = ,  
  const char * inpname = ,
  const char * outputfile = )
{
  const int nConfidenceLevel = 4;
  int confidenceLevel [nConfidenceLevel] = {68, 90, 95, 99 };

  //
  // Define grid
  //
    
  using namespace t2kana;//just for kpi
      
  // deltacp, mixing13 range and step
  double mixing13_min   = 0.0;
  double mixing13_max   = 0.4;
  double mixing13_step  = 0.004;
  double deltacp_min    = -1 * conv::kpi;
  double deltacp_max    = conv::kpi;
  double deltacp_step   = 2*conv::kpi/90.;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing13_nbins  = 1 + (int) ((mixing13_max       - mixing13_min      ) / mixing13_step);
  int deltacp_nbins   = 1 + (int) ((deltacp_max        - deltacp_min       ) / deltacp_step    );
  
  cout << "mixing13_nbins = " << mixing13_nbins << ", deltacp_nbins = " << deltacp_nbins<<endl;
  
  TH2D * hDeltaChi2Surface[nConfidenceLevel];
  for(int icl = 0; icl < nConfidenceLevel; icl++) { 
    const char * histname = Form("hDeltaChi2Surface%02d", confidenceLevel[icl]);
    hDeltaChi2Surface[icl] =
       new TH2D(histname, "",
          mixing13_nbins,
          mixing13_min - 0.5*mixing13_step,
          mixing13_max + 0.5*mixing13_step,
          deltacp_nbins,
	  deltacp_min - 0.5*deltacp_step,
          deltacp_max + 0.5*deltacp_step);
  }


  //
  // Load inputs
  //
  TFile * fin[nruns];
  TH2D * hDeltaChi2SurfacePartial[nruns][nConfidenceLevel];
  for(int irun = 0; irun < nruns; irun++) {
    fin[irun] = new TFile( Form("%s%d.root",inpname,irun), "read");
    cout << "Getting histograms from file: " <<  Form("%s%d.root",inpname,irun) << endl;
    fin[irun]->cd();
    hDeltaChi2SurfacePartial[irun][0] = (TH2D*) fin[irun]->Get("hDeltaChi2Surface68");
    hDeltaChi2SurfacePartial[irun][1] = (TH2D*) fin[irun]->Get("hDeltaChi2Surface90");
    hDeltaChi2SurfacePartial[irun][2] = (TH2D*) fin[irun]->Get("hDeltaChi2Surface95");
    hDeltaChi2SurfacePartial[irun][3] = (TH2D*) fin[irun]->Get("hDeltaChi2Surface99");
  }

  cout << "Combining histograms" << endl;

  //
  // Add histograms (each covers part of the mixing13, deltacp space)
  //
  for(int imixing13_bin = 1; imixing13_bin <= mixing13_nbins; imixing13_bin++) {
    for(int ideltacp_bin = 1; ideltacp_bin <= deltacp_nbins; ideltacp_bin++) {
      for(int icl = 0; icl < nConfidenceLevel; icl++) {
          double dchisq = 0;
          for(int irun = 0; irun < nruns; irun++) {
             dchisq += hDeltaChi2SurfacePartial[irun][icl]->GetBinContent(imixing13_bin,ideltacp_bin);
          }//irun
          hDeltaChi2Surface[icl]->SetBinContent(imixing13_bin,ideltacp_bin,dchisq);
      }//icl
    }//ideltacp_bin
  }//imixing13_bin

  //
  // Save combined histograms
  //

  cout << "Saving histograms to " << outputfile << endl;

  TFile * fout = new TFile(outputfile, "RECREATE");
  for(int icl=0; icl<nConfidenceLevel; icl++) {
     hDeltaChi2Surface[icl]->Write();
  }
  fout->Close();

}
