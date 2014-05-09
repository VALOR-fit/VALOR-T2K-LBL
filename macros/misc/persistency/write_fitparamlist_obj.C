void write_fitparamlist_obj(
   const char * filename = "./fp.root", const char * objname = "fit_param_list")
{
  using namespace valor;

  // get a param list
  FitParamList * fparm = util::NamedParameterList("nuosc3f_numudisapp_fix4+nuisance");

  // tweak a couple of params
  fparm->Update(conv::kPhypDm23sq,   2.345E-3);
  fparm->Update(conv::kPhypMixing23, 0.98765);

  // save
  TFile f(filename,"recreate");
  fparm->Write(objname);
  f.Close();
}
