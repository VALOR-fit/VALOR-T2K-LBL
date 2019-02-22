double get_mass(double h, double r)
{
  double v = TMath::Pi() * r * r * h; //m^3
  const double density = 997; //kg/m^3
  double mass = v * density / 1E6;

  double fid_v = TMath::Pi() * (r-2) * (r-2) * (h-4);
  double fid_mass = fid_v * density / 1E6;
  double surface_area = 2 * TMath::Pi() * r * h + 2 * TMath::Pi() * r * r;

  cout << "\tDiameter " << r*2 << " m Height " << h << " m" << endl;
  cout << "\tTotal mass " << mass << " kton" << endl;
  cout << "\t2 m cut fiducial mass " << fid_mass << " kton" << endl;
  cout << "\tSurface area " << surface_area << " m^2" << endl;
  return mass;
}

void detector_mass()
{
  double d, h;

  cout << "SK total (total should be 50 kton)" << endl;
  get_mass(42,39/2);
  cout << "SK ID (total should be 32 kton, FV should be 22.5 kton)" << endl;
  get_mass(36.2,33.8/2);
  cout << "HK total (total should be 258 kton)" << endl;
  get_mass(60,74/2);
  cout << "HK ID (FV should be 187 kton)" << endl;
  get_mass(54.8,70.8/2);
}
