
void convert_nd_dst_tree_to_varray_format(
 const char * config_file = "dune_ndfit.xml",
 const char * output_filename = "valor_templates.root",
 const char * tag = "nd"
 )
{
  using namespace std;
  using namespace valor;
  
  TString config = TString::Format("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), config_file);
  
  // Get ND model & configure
  FitModelABC * model = new SDModel;
  model->Read(config.Data(),"forcedst");
  
  // Get array with all VArrays in
  TObjArray output_arr;
  model->Write(&output_arr,tag,"mct_vf+SB");
  
  // Open an output file
  
  TFile outfile(output_filename,"RECREATE");
  if (!outfile.IsOpen()) {
    cerr<<"Failed to open outfile file "<<output_filename<<endl;
    exit(-1);
  }
  
  cout<<"Writing arrays to "<<output_filename<<endl;
  
  // Loop over the objects, select the nominal templates only
  for (int i = 0 ; i < output_arr.GetEntries(); i++) {
    TObject * obj = output_arr.At(i);
    TString name = obj->GetName();
    if (name.Contains("mctnom_ob_vf") || name.Contains("mctnom_sb_vf")) {
      VALORArray * va = dynamic_cast<VALORArray*>(obj);
      if (!va) {
        std::cerr<<"Warning, not able to convert "<<name<<" to VALORArray"<<endl;
        continue;
      }
      outfile.cd();
      va->Write(name);
      cout<<"Written: "<<name<<endl;
    }
  }
  
  outfile.Close();
  return;
}
    
