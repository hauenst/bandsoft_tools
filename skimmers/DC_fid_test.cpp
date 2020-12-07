//This is a simple code to test the DC_Fid
#include "DC_fiducial.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"

#include "reader.h"
#include "bank.h"
#include "clas12fiducial.h"

#include "RCDB/Connection.h"

#include "constants.h"


bool inbending = true;
bool outbending = false; 


using namespace std;
int main()  {

  int sec =0;
  int lay =0;

  double x =0;  //in cm
  double y =0;  //in cm

  bool bending = false;

  cout << "Pick a point to test, x, y, sec, lay, bending "<< endl;

  cin >> x >> y >> sec >> lay >> bending;

  DCFiducial dc_acc;

  bool DC_fid = dc_acc.DC_e_fid(x, y, sec, lay, bending);


  cout << "Check: "<< DC_fid << endl;


}
