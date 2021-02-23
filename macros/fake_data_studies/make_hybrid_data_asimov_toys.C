
//___________________________________________________________________________
void make_hybrid_data_asimov_toys(
 const char * config_file      = "analysis/t2k_analysis_configuration_oa2020_run1-10_v1_postfit.xml",
 const char * tag              = "run1-10_hybrid_4s_data_1s_AsimovBFNO",
 const bool use_simple_binning = false
 )
{
  using namespace valor;

  CompoundFitElements* elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file));
  FitModelABC* model = elements->FitModel();
  FitDataABC * data  = elements->FitData ();
  CompoundFitModel* cm = dynamic_cast<CompoundFitModel*>(model);
  CompoundFitData * cd = dynamic_cast<CompoundFitData *>(data );

  const char* opt_binning = "";
  if(use_simple_binning) { opt_binning = "+SB"; }

  // The data is stored as a tree (to retain event-by-event info), the location of which is specified in the
  // analysis config file, and is read at runtime. This just checks that it actually exists.
  if(!cd->DataExists(opt_binning)) {
    std::cerr << "ERROR: no data specified in the analysis configuration file" << std::endl;
    util::Die(-1);
  }
  std::cout << "Pure data PDF integrals:\n" << cd->GetPDFIntegralAsString(opt_binning) << std::endl;

  // Generate the Asimov BF NO toy according to hardcoded osc param values
  FitParamList * fp = util::GetFitParameterList("3f_osc");
  if (fp->Exists(conv::kPhypMixing12)) { fp->Update(conv::kPhypMixing12, 0.307   ); }
  if (fp->Exists(conv::kPhypMixing13)) { fp->Update(conv::kPhypMixing13, 0.021825); }
  if (fp->Exists(conv::kPhypMixing23)) { fp->Update(conv::kPhypMixing23, 0.555   ); }
  if (fp->Exists(conv::kPhypDm12sq  )) { fp->Update(conv::kPhypDm12sq  , 7.53E-5 ); }
  if (fp->Exists(conv::kPhypDm23sq  )) { fp->Update(conv::kPhypDm23sq  , 2.488E-3); }
  if (fp->Exists(conv::kPhypCPV     )) { fp->Update(conv::kPhypCPV     , -2.0106 ); }
  cm->Generate(*fp, opt_binning);
  std::cout << "Pure Asimov PDF integrals:\n" << cm->GetPDFIntegralAsString(opt_binning) << std::endl;

  // Write to VALORArrays
  TObjArray * array_data = new TObjArray();
  array_data->SetOwner(true);
  const char* write_opt_data = "data_vf";
  if(use_simple_binning) { write_opt_data = "data_vf+SB"; }
  cd->Write(array_data, "data", write_opt_data);

  TObjArray * array_asimov = new TObjArray();
  array_asimov->SetOwner(true);
  const char* write_opt_asimov = "pdft_vf";
  if(use_simple_binning) { write_opt_asimov = "pdft_vf+SB"; }
  cm->Write(array_asimov, "asimov", write_opt_asimov);

  // Create hybrids made of N-1 data samples and 1 Asimov sample
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel* this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));

    TString name_asimov_ob = TString::Format("asimov_pdft_ob_vf_%s", this_model->Name().c_str());
    TString name_asimov_sb = TString::Format("asimov_pdft_sb_vf_%s", this_model->Name().c_str());
    VALORArray* pdf_asimov_ob = (VALORArray*) array_asimov->FindObject(name_asimov_ob.Data());
    VALORArray* pdf_asimov_sb = nullptr;
    if(use_simple_binning) { pdf_asimov_sb = (VALORArray*) array_asimov->FindObject(name_asimov_sb.Data()); }
    if(!pdf_asimov_ob) {
      std::cerr << "ERROR: Couldn't retrieve the Asimov array [" << name_asimov_ob.Data() << "]" << std::endl;
      util::Die(-1);
    }
    if(use_simple_binning && !pdf_asimov_sb) {
      std::cerr << "ERROR: Couldn't retrieve the Asimov array [" << name_asimov_sb.Data() << "]" << std::endl;
      util::Die(-1);
    }

    for(unsigned int idata = 0; idata < cd->GetNSubData(); idata++) {
      GenericSingleSampleFitData* this_data = dynamic_cast<GenericSingleSampleFitData*>(cd->GetSubData(idata));

      TString name_data_ob = TString::Format("data_data_ob_vf_%s",this_data->Name().c_str());
      TString name_data_sb = TString::Format("data_data_sb_vf_%s",this_data->Name().c_str());
      VALORArray* pdf_data_ob = (VALORArray*) array_data->FindObject(name_data_ob.Data());
      VALORArray* pdf_data_sb = nullptr;
      if(use_simple_binning) { pdf_data_sb = (VALORArray*) array_data->FindObject(name_data_sb.Data()); }
      if(!pdf_data_ob) {
        std::cerr << "ERROR: Couldn't retrieve the data array [" << name_data_ob.Data() << "]" << std::endl;
        util::Die(-1);
      }
      if(use_simple_binning && !pdf_data_sb) {
        std::cerr << "ERROR: Couldn't retrieve the data array [" << name_data_sb.Data() << "]" << std::endl;
        util::Die(-1);
      }

      VALORArray* pdf_overwrite_ob = pdf_data_ob;
      VALORArray* pdf_overwrite_sb = pdf_data_sb;
      if(idata == imodel) {
        pdf_overwrite_ob = pdf_asimov_ob;
        pdf_overwrite_sb = pdf_asimov_sb;
      }

      this_data->OverwriteData(pdf_overwrite_ob, pdf_overwrite_sb);
    }

    // Save the hybrid toy to a binary file
    std::string filename = (std::string)tag + "_" + this_model->Name() + ".bin";
    std::cout << "Made hybrid fake data toy with " << cm->GetNSubModels()-1
              << " data samples and 1 Asimov sample (" << this_model->Name()
              << "). \n Saving to file " << filename
              << "\n Hybrid toy PDF integrals:\n"
              << cd->GetPDFIntegralAsString(opt_binning) << std::endl;
    util::ToyStream* toy_stream = cd->InitToyStream(filename, "write", write_opt_data);
    cd->Write(toy_stream);
    delete toy_stream;
  }

  std::cout << "Done!" << std::endl;
}
    
