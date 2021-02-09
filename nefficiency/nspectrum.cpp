#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TChain.h"
#include "TClonesArray.h"

#include "reader.h"
#include "bank.h"
#include "clas12fiducial.h"
#include "e_pid.h"
#include "constants.h"
#include "bandhit.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;

bool pointsToBand(double theta,double phi,double z_m);

int main(int argc, char** argv) {

  if( argc < 3 ){
    cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [inputFiles]... \n\n";
    return -1;
  }
  
  //Define Variables
  bool goodneutron;
  int Runno, pMult, scinHits, nleadindex;
  double Ebeam, gated_charge, livetime, starttime, current;
  int pIndex[100];
  double pChi2pid[100], p_vtx[100], p_vty[100], p_vtz[100], p_p[100], theta_p[100], phi_p[100],pBeta[100], pPid[100], pCharge[100], hit_pindex[100], hit_detid[100];
  clashit * eHit = new clashit;
  TClonesArray* nHits = new TClonesArray("bandhit");
  char temp[100];

  TChain inChain("skim");
 //Creat input tree
  for(int k = 2; k < argc; k++){
    inChain.Add(argv[k]);
  }
  //TFile * inFile = new TFile(argv[2]);
  //TTree * inTree = (TTree*)inFile->Get("skim");
  // 	Event branches:
  inChain.SetBranchAddress("Runno"              ,&Runno         );
  inChain.SetBranchAddress("Ebeam"              ,&Ebeam         );
  inChain.SetBranchAddress("gated_charge"       ,&gated_charge  );
  inChain.SetBranchAddress("livetime"           ,&livetime      );
  inChain.SetBranchAddress("starttime"          ,&starttime     );
  inChain.SetBranchAddress("current"            ,&current       );
  //	Electron branches:
  inChain.SetBranchAddress("eHit"		,&eHit	        );
  //    Charged particle branches:
  inChain.SetBranchAddress("pMult"   		,&pMult  	);
  inChain.SetBranchAddress("pIndex",pIndex);
  inChain.SetBranchAddress("pPid"    		,pPid   	);
  inChain.SetBranchAddress("pCharge" 		,pCharge	);
  inChain.SetBranchAddress("pChi2pid"		,pChi2pid	);
  inChain.SetBranchAddress("p_vtx"   		,p_vtx  	);
  inChain.SetBranchAddress("p_vty"   		,p_vty  	);
  inChain.SetBranchAddress("p_vtz"   		,p_vtz  	);
  inChain.SetBranchAddress("p_p"     		,p_p    	);
  inChain.SetBranchAddress("theta_p" 		,theta_p	);
  inChain.SetBranchAddress("phi_p"   		,phi_p  	);
  inChain.SetBranchAddress("pBeta"   		,pBeta  	);
  //    Scintilator branches:
  inChain.SetBranchAddress("scinHits"  		,&scinHits  	);
  inChain.SetBranchAddress("hit_pindex"		,hit_pindex  	);
  inChain.SetBranchAddress("hit_detid"		,hit_detid  	);
  //    Neutron branches:
  inChain.SetBranchAddress("nHits"		,&nHits	        );
  inChain.SetBranchAddress("goodneutron"        ,&goodneutron   );
  inChain.SetBranchAddress("nleadindex"         ,&nleadindex    );

  // Create output tree
  TFile * outFile = new TFile(argv[1],"RECREATE");

  //Creating TH2D from tree
  vector<TH1*> hist_list_1;
  TH1D * h_ToF = new TH1D("ToF","ToF;ToF;Events",100,0,60);
  hist_list_1.push_back(h_ToF);
  TH1D * h_ToFpM = new TH1D("ToFpM","ToFpM;ToFpM;Events",400,-40,120);
  hist_list_1.push_back(h_ToFpM);
  TH1D * h_Beta = new TH1D("Beta","Beta;Beta;Events",100,0.5,1);
  hist_list_1.push_back(h_Beta);
  TH1D * h_EDep = new TH1D("EDep","EDep;EDep;Events",100,0,10);
  hist_list_1.push_back(h_EDep);
  TH1D * h_Mom = new TH1D("Mom","Mom;Mom;Events",100,0,2);
  hist_list_1.push_back(h_Mom);

  vector<TH2*> hist_list_2;
  TH2D * h_ToF_EDep = new TH2D("ToF_EDep","ToF_EDep;ToF;EDep;Events",100,0,60,100,0,10);
  hist_list_2.push_back(h_ToF_EDep);
  TH2D * h_ToFpM_EDep = new TH2D("ToFpM_EDep","ToFpM_EDep;ToFpM;EDep;Events",400,-40,120,100,0,10);
  hist_list_2.push_back(h_ToFpM_EDep);
  TH2D * h_ToFpM_Mom = new TH2D("ToFpM_Mom","ToFpM_Mom;ToFpM;Mom;Events",400,-40,120,100,0,2);
  hist_list_2.push_back(h_ToFpM_Mom);
  TH2D * h_ToFpM_i_Beta = new TH2D("ToFpM_Beta","ToFpM_Beta;ToFpM;Beta;Events",400,-40,120,100,0,1);
  hist_list_2.push_back(h_ToFpM_i_Beta);
  TH2D * h_Mom_Mom = new TH2D("Mom_Mom","Mom_Mom;Mom;Mom;Events",100,0,2,100,0,2);
  TH2D * h_Mom_EDep = new TH2D("Mom_EDep","Mom_EDep;Mom;EDep;Events",100,0,2,100,0,10);
  hist_list_2.push_back(h_Mom_EDep);
  TH2D * h_X_Y = new TH2D("X_Y","X_Y;X;Y;Events",100,-150,150,100,-100,100);
  hist_list_2.push_back(h_X_Y);

  int fin = inChain.GetEntries();
  for(int i = 0; i < fin; i++){

    //Display completed  
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }
    
    eHit->Clear();
    nHits->Clear();
    inChain.GetEntry(i);

    if(!goodneutron){ continue; }
    bandhit * this_nHit = (bandhit*)nHits->At(nleadindex);
    double ToF = this_nHit->getTof();
    double DL = this_nHit->getDL().Mag();
    double ToM = ToF * 100 / DL;
    double Beta = this_nHit->getBeta();
    double EDep = this_nHit->getEdep() / 1100;
    TVector3 Mom = this_nHit->getMomentumN();
    double Gamma = 1/sqrt(1-(Beta*Beta));
    double calcMom = Beta * Gamma * mN;
    double X = this_nHit->getX();
    double Y = this_nHit->getY();


    if(ToM< 6){continue;}
    if(ToM>15){continue;}
    h_ToF_EDep->	Fill(ToF	,EDep	);
    h_ToFpM_EDep->	Fill(ToM	,EDep	);
    h_ToFpM_Mom->	Fill(ToM	,Mom.Mag()	);
    h_ToFpM_i_Beta->	Fill(ToM	,1/Beta	);
    h_Mom_Mom-> 	Fill(Mom.Mag()	,calcMom	);
    h_Mom_EDep->	Fill(Mom.Mag()	,EDep	);
    h_X_Y->		Fill(X		,Y	);

    if(EDep < 2.5){continue;}

    h_ToF->		Fill(ToF		);
    h_ToFpM->		Fill(ToM		);
    h_Beta->		Fill(Beta		);
    h_EDep->		Fill(EDep		);
    h_Mom->		Fill(Mom.Mag()		);

  }


  //inFile->Close();
  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
  }
  outFile->Close();
  cout<<"Finished making file: "<<argv[1]<<"\n";
  return 0;
}

bool pointsToBand(double theta,double phi,double z_m){
	//double z = z_m*100; // from m to cm
	double z = z_m;

	// Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
	double thickness  = 7.2;                                // thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

	// Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = (zUpst + 5*thickness) - z_m;

	double rho   = zDown/cos(theta);
	double xDown = rho*sin(theta)*cos(phi);
	double yDown = rho*sin(theta)*sin(phi);

	double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	// Sector boundaries
	double topSec1  = globalY + 13*thickness;
	double topSec2  = globalY + 10*thickness;
	double topSec34 = globalY +  3*thickness;
	double topSec5  = globalY -  3*thickness;
	double downSec5 = globalY -  5*thickness;

	if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

	if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
			(yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
			(yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
			(yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
			){
	  return 1;
	}
	return 0;
}
