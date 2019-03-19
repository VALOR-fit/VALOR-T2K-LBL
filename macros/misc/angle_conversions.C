#include <iostream>
#include <string>
#include "TMath.h"

double sinsqtheta_to_deg(double x)
{
  x = TMath::Sqrt(x);
  x = TMath::ASin(x);
  //cout << "sinsqtheta_to_deg\t" << x << " rad\t" << TMath::RadToDeg() * x << " deg" << endl;
  return TMath::RadToDeg() * x;
}

double sinsqtheta_to_rad(double x)
{
  x = TMath::Sqrt(x);
  x = TMath::ASin(x);
  //cout << "sinsqtheta_to_rad\t" << x << " rad\t" << TMath::RadToDeg() * x << " deg" << endl;
  return x;
}

double sinsq2theta_to_deg(double x)
{
  x = TMath::Sqrt(x);
  x = TMath::ASin(x);
  x /= 2;
  //cout << "sinsq2theta_to_deg\t" << x << " rad\t" << TMath::RadToDeg() * x << " deg" << endl;
  return TMath::RadToDeg() * x;
}

double sinsq2theta_to_rad(double x)
{
  x = TMath::Sqrt(x);
  x = TMath::ASin(x);
  x /= 2;
  //cout << "sinsq2theta_to_rad\t" << x << " rad\t" << TMath::RadToDeg() * x << " deg" << endl;
  return x;
}

double deg_to_sinsqtheta(double x)
{
  x *= TMath::DegToRad();
  x = TMath::Sin(x);
  x = TMath::Power(x, 2);
  //cout << "deg_to_sinsqtheta\t" << x << endl;
}

double rad_to_sinsqtheta(double x)
{
  x = TMath::Sin(x);
  x = TMath::Power(x, 2);
  //cout << "rad_to_sinsqtheta\t" << x << endl;
}

double deg_to_sinsq2theta(double x)
{
  x *= TMath::DegToRad();
  x *= 2;
  x = TMath::Sin(x);
  x = TMath::Power(x, 2);
  //cout << "deg_to_sinsq2theta\t" << x << endl;
}

double rad_to_sinsq2theta(double x)
{
  x *= 2;
  x = TMath::Sin(x);
  x = TMath::Power(x, 2);
  //cout << "rad_to_sinsq2theta\t" << x << endl;
}

void angle_conversions(double x = 0.1, int mode = 0)
{
  const int nmodes = 4;
  assert(mode >= 0 && mode < nmodes);
  string mode2str[nmodes] = {"sinsq2theta", "sinsqtheta", "degrees", "radians"};
  cout << "Input is " << x << " " << mode2str[mode] << endl;
  double sinsqtheta, sinsq2theta, deg, rad;
  switch(mode) {
  case 0:
    sinsq2theta = x;
    deg = sinsq2theta_to_deg(sinsq2theta);
    rad = sinsq2theta_to_rad(sinsq2theta);
    sinsqtheta = rad_to_sinsqtheta(rad);
    break;
  case 1:
    sinsqtheta = x;
    deg = sinsqtheta_to_deg(sinsqtheta);
    rad = sinsqtheta_to_rad(sinsqtheta);
    sinsq2theta = rad_to_sinsq2theta(rad);
    break;
  case 2:
    deg = x;
    rad = TMath::DegToRad() * deg;
    sinsqtheta  = deg_to_sinsqtheta(deg);
    sinsq2theta = deg_to_sinsq2theta(deg);
    break;
  case 3:
    rad = x;
    deg = TMath::RadToDeg() * rad;
    sinsqtheta  = rad_to_sinsqtheta(rad);
    sinsq2theta = rad_to_sinsq2theta(rad);
    break;
  }

  cout << "sinsq2theta    " << sinsq2theta << endl
       << "sinsqtheta     " << sinsqtheta << endl
       << "degrees        " << deg << endl
       << "radians        " << rad << endl;
}
