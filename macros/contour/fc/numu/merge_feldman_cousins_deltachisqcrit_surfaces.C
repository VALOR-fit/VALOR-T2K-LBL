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
  const int nruns, 
  const char * inpname, 
  const char * outputfile)
{
  const int nConfidenceLevel = 4;
  int confidenceLevel [nConfidenceLevel] = {68, 90, 95, 99 };

  //
  // Define grid
  //
          
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.5;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
    
  TH2D * hDeltaChi2Surface[nConfidenceLevel];
  for(int icl = 0; icl < nConfidenceLevel; icl++) { 
    const char * histname = Form("hDeltaChi2Surface%02d", confidenceLevel[icl]);
    hDeltaChi2Surface[icl] =
       new TH2D(histname, "",
          mixing23_nbins,
          mixing23_min - 0.5*mixing23_step,
          mixing23_max + 0.5*mixing23_step,
          dm23sq_nbins,
          dm23sq_min - 0.5*dm23sq_step,
          dm23sq_max + 0.5*dm23sq_step);
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
  // Add histograms (each covers part of the sinsq2theta, dmsq space)
  //
  for(int imixing23_bin = 1; imixing23_bin <= mixing23_nbins; imixing23_bin++) {
    for(int idm23sq_bin = 1; idm23sq_bin <= dm23sq_nbins; idm23sq_bin++) {
      for(int icl = 0; icl < nConfidenceLevel; icl++) {
          double dchisq = 0;
          for(int irun = 0; irun < nruns; irun++) {
             dchisq += hDeltaChi2SurfacePartial[irun][icl]->GetBinContent(imixing23_bin,idm23sq_bin);
          }//irun
          hDeltaChi2Surface[icl]->SetBinContent(imixing23_bin,idm23sq_bin,dchisq);
      }//icl
    }//idm23sq_bin
  }//imixing23_bin

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
