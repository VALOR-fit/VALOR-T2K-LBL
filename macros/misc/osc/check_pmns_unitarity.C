
void check_pmns_unitarity
(
 double sinsq2theta12 = 0.38,
 double sinsq2theta13 = 0.96,
 double sinsq2theta23 = 0.13,
 double deltaCP = 1.6
)
{
  enum {
    re,
    im
  };

  //calculate sin(mixingAngle) from sin^2(2*mixingAngle) 
  double s12 = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta12)));
  double s13 = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta13)));
  double s23 = sqrt(0.5 * (1 - sqrt(1 - sinsq2theta23)));
  //calculate cos(mixingAngle) from sin(mixingAngle)     
  double c12 = sqrt(1 - s12 * s12);
  double c13 = sqrt(1 - s13 * s13);
  double c23 = sqrt(1 - s23 * s23);

  double U[3][3][2];
  double Udagger[3][3][2];
  double UtimesUdagger[3][3][2];

  //U is the PMNS matrix
  //first row has indices [0][0], [0][1], [0][2]
  //first column has indices [0][0], [1][0], [2][0]
  U[0][0][re] = c12*c13;
  U[0][0][im] = 0.0;
  U[0][1][re] = c13*s12;
  U[0][1][im] = 0.0;
  U[0][2][re] = s13*cos(deltaCP);
  U[0][2][im] = -s13*sin(deltaCP);
  U[1][0][re] = -c23*s12 - s13*s23*c12*cos(deltaCP);
  U[1][0][im] = -s13*s23*c12*sin(deltaCP);
  U[1][1][re] = c12*c23 - s12*s13*s23*cos(deltaCP);
  U[1][1][im] = -s12*s13*s23*sin(deltaCP);
  U[1][2][re] = c13*s23;
  U[1][2][im] = 0.0;
  U[2][0][re] = s12*s23 - s13*c12*c23*cos(deltaCP);
  U[2][0][im] = -s13*c12*c23*sin(deltaCP);
  U[2][1][re] = -s23*c12 - s12*c23*s13*cos(deltaCP);
  U[2][1][im] = -s12*c23*s13*sin(deltaCP);
  U[2][2][re] = c13*c23;
  U[2][2][im] = 0.0;
  
  //Udagger is the transpose of the complex conjugate of the PMNS matrix
  for(int i=0; i<3; i++)
    {
    for(int j=0; j<3; j++)
      {
      Udagger[i][j][re] = U[j][i][re];
      Udagger[i][j][im] = -U[j][i][im];
      }
    }

  //check unitarity of PMNS matrix
  for(int i=0; i<3; i++)
    {
    for(int j=0; j<3; j++)
      {
      UtimesUdagger[i][j][re] = U[i][0][re]*Udagger[0][j][re] - U[i][0][im]*Udagger[0][j][im]
                              + U[i][1][re]*Udagger[1][j][re] - U[i][1][im]*Udagger[1][j][im]
	                      + U[i][2][re]*Udagger[2][j][re] - U[i][2][im]*Udagger[2][j][im];
      UtimesUdagger[i][j][im] = U[i][0][re]*Udagger[0][j][im] + U[i][0][im]*Udagger[0][j][re]
                              + U[i][1][re]*Udagger[1][j][im] + U[i][1][im]*Udagger[1][j][re]
	                      + U[i][2][re]*Udagger[2][j][im] + U[i][2][im]*Udagger[2][j][re];
      }
    }
 
  cout<<endl;
  cout<<"U is the PMNS matrix"<<endl;
  cout<<endl;
  cout<<"These are the real parts of U times Udagger"<<endl;
  cout<<endl;

  for(int i=0; i<3; i++)
    {
    for(int j=0; j<3; j++)
      {
	cout<<UtimesUdagger[i][j][re]<<"          ";
      }
    cout<<endl;
    }

  cout<<endl;

  cout<<"These are the imaginary parts of U times Udagger"<<endl;
  cout<<endl;

  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  cout<<UtimesUdagger[i][j][im]<<"          ";
	}
      cout<<endl;
    }
  
cout<<endl;

}//end of check_pmns_unitarity()
