void write_physparamlist_obj(
   const char * filename = "./physlist.root", const char * objname = "phys_param_list")
{
  using namespace t2kana;

  // create a physics parameter list
  PhysParamList pl;

  // set a few params
  pl.SetValue(conv::kPhypDm23sq,                2.345E-3);
  pl.SetValue(conv::kPhypMixing23,              0.987654);
  pl.SetValue(conv::kNuisSuperkEff1RmuCCQErc0u, 1.010101);
  pl.SetValue(conv::kNuisSuperkEff1RmuCCQErc1u, 1.234567);
  pl.SetValue(conv::kNuisSuperkEff1RmuCCQErc2u, 0.987654);
  pl.SetValue(conv::kNuisBANFF00,               1.020202);
  pl.SetValue(conv::kNuisBANFF01,               1.030303);
  pl.SetValue(conv::kNuisBANFF02,               1.040404);
  pl.SetValue(conv::kNuisBANFF03,               0.999999);

  // save
  TFile f(filename,"recreate");
  pl.Write(objname);
  f.Close();
}
