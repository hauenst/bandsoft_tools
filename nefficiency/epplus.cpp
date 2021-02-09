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
#include "TFitResult.h"
#include "TFitResultPtr.h"

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

  if( argc < 4 ){
    cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [num Par] [inputFiles]... \n\n";
    return -1;
  }
  
  //Define Variables
  bool goodneutron;
  int Runno, pMult, scinHits, nleadindex, aiP, nMult, eventnumber;
  double Ebeam, gated_charge, livetime, starttime, current, p_pmiss, theta_pmiss, phi_pmiss, mmiss;
  int pIndex[100];
  double pChi2pid[100], p_vtx[100], p_vty[100], p_vtz[100], p_p[100], theta_p[100], phi_p[100],pBeta[100], pPid[100], pCharge[100], hit_pindex[100], hit_detid[100];
  clashit * eHit = new clashit;
  TClonesArray* nHits = new TClonesArray("bandhit");
  char temp[100];

  int numpar = atoi(argv[2]);

  TChain inChain("skim");
 //Creat input tree
  for(int k = 3; k < argc; k++){
    inChain.Add(argv[k]);
  }
  // 	Event branches:
  inChain.SetBranchAddress("Runno"              ,&Runno         );
  inChain.SetBranchAddress("eventnumber"        ,&eventnumber   );
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
  inChain.SetBranchAddress("nMult"		,&nMult	        );
  inChain.SetBranchAddress("goodneutron"        ,&goodneutron   );
  inChain.SetBranchAddress("nleadindex"         ,&nleadindex    );
  //   Missing momentum branches in (e,e'p)n
  inChain.SetBranchAddress("aiP"		,&aiP	        );
  inChain.SetBranchAddress("p_pmiss"		,&p_pmiss	);
  inChain.SetBranchAddress("theta_pmiss"       	,&theta_pmiss	);
  inChain.SetBranchAddress("phi_pmiss"		,&phi_pmiss     );
  inChain.SetBranchAddress("mmiss"		,&mmiss	        );

  // Create output tree
  TFile * outFile = new TFile(argv[1],"RECREATE");

  //Creating TH2D from tree
  vector<TH1*> hist_list_1;
  TH1D * h_vtz_11 = new TH1D("vtz_11","vtz_11;vtz;Events",100,-10,10);
  hist_list_1.push_back(h_vtz_11);
  TH1D * h_vtz_0 = new TH1D("vtz_0","vtz_0;vtz;Events",100,-10,10);
  hist_list_1.push_back(h_vtz_0);
  TH1D * h_vtz_E = new TH1D("vtz_E","vtz_E;vtz;Events",100,-10,10);
  hist_list_1.push_back(h_vtz_E);
  TH1D * h_vtz_diff_11 = new TH1D("vtz_diff_11","vtz_diff_11;vtz_diff;Events",100,-10,10);
  hist_list_1.push_back(h_vtz_diff_11);
  TH1D * h_vtz_diff_0 = new TH1D("vtz_diff_0","vtz_diff_0;vtz_diff;Events",100,-10,10);
  hist_list_1.push_back(h_vtz_diff_0);
  TH1D * h_vtz_diff_E = new TH1D("vtz_diff_E","vtz_diff_E;vtz_diff;Events",100,-10,10);
  hist_list_1.push_back(h_vtz_diff_E);

  TH1D * h_theta_e_11 = new TH1D("theta_e_11","theta_e_11;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_e_11);
  TH1D * h_theta_e_0 = new TH1D("theta_e_0","theta_e_0;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_e_0);
  TH1D * h_theta_e_E = new TH1D("theta_e_E","theta_e_E;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_e_E);

  TH1D * h_theta_p_11 = new TH1D("theta_p_11","theta_p_11;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_p_11);
  TH1D * h_theta_p_0 = new TH1D("theta_p_0","theta_p_0;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_p_0);
  TH1D * h_theta_p_E = new TH1D("theta_p_E","theta_p_E;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_p_E);

  TH1D * h_mom_11 = new TH1D("mom_11","mom_11;mom;Events",100,0,5);
  hist_list_1.push_back(h_mom_11);
  TH1D * h_mom_0 = new TH1D("mom_0","mom_0;mom;Events",100,0,5);
  hist_list_1.push_back(h_mom_0);
  TH1D * h_mom_E = new TH1D("mom_E","mom_E;mom;Events",100,0,5);
  hist_list_1.push_back(h_mom_E);
  TH1D * h_theta_11 = new TH1D("theta_11","theta_11;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_11);
  TH1D * h_theta_0 = new TH1D("theta_0","theta_0;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_0);
  TH1D * h_theta_E = new TH1D("theta_E","theta_E;theta;Events",180,0,180);
  hist_list_1.push_back(h_theta_E);
  TH1D * h_phi_11 = new TH1D("phi_11","phi_11;phi;Events",360,-180,180);
  hist_list_1.push_back(h_phi_11);
  TH1D * h_phi_0 = new TH1D("phi_0","phi_0;phi;Events",360,-180,180);
  hist_list_1.push_back(h_phi_0);
  TH1D * h_phi_E = new TH1D("phi_E","phi_E;phi;Events",360,-180,180);
  hist_list_1.push_back(h_phi_E);
  TH1D * h_mom_rho_11 = new TH1D("mom_rho_11","mom_rho_11;mom_rho;Events",100,0,5);
  hist_list_1.push_back(h_mom_rho_11);
  TH1D * h_mom_rho_0 = new TH1D("mom_rho_0","mom_rho_0;mom_rho;Events",100,0,5);
  hist_list_1.push_back(h_mom_rho_0);
  TH1D * h_mom_rho_E = new TH1D("mom_rho_E","mom_rho_E;mom_rho;Events",100,0,5);
  hist_list_1.push_back(h_mom_rho_E);

  TH1D * h_1 = new TH1D("h_1","h1;1;Events",1,0,1);
  hist_list_1.push_back(h_1);
  TH1D * h_pid = new TH1D("h_pid","h_pid;pid;Events",6000,-3000,3000);
  hist_list_1.push_back(h_pid);
  TH1D * h_mult = new TH1D("h_mult","h_mult;mult;Events",10,0,10);
  hist_list_1.push_back(h_mult);

  vector<TH2*> hist_list_2;

  int fin = inChain.GetEntries();
  for(int i = 0; i < fin; i++){
    eHit->Clear();
    nHits->Clear();
    inChain.GetEntry(i);

    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }    
    if(mmiss < 0.85){ continue; }
    if(mmiss > 1.05){ continue; }
    //cout<<"here\n";

    h_mult->Fill(pMult);
    if(pMult!=numpar){continue;}
    h_1->Fill(1);

    TVector3 vbeam(0,0,Ebeam);
    TVector3 ve;
    ve.SetMagThetaPhi(eHit->getMomentum(),eHit->getTheta(),eHit->getPhi());
    TVector3 vp;
    vp.SetMagThetaPhi(p_p[aiP],theta_p[aiP],phi_p[aiP]);
    TVector3 vmiss;
    vmiss.SetMagThetaPhi(p_pmiss,theta_pmiss,phi_pmiss);

    for(int j = 0; j < pMult; j++){
      if(j != aiP){
	
	h_pid->Fill(pPid[j]);
	//cout<<j<<" "<<pPid[j]<<"  "<<pBeta[j]<<"\n";
	TVector3 vX;
	vX.SetMagThetaPhi(p_p[j],theta_p[j],phi_p[j]);
	double theta_e_X = ve.Angle(vX) * 180 / M_PI;
	double theta_p_X = vp.Angle(vX) * 180 / M_PI;
	double z = eHit->getVtz()-p_vtz[j];
	if(pPid[j]==11){
	  h_vtz_diff_11->Fill(z);
	  h_vtz_11->Fill(p_vtz[j]);
	  h_theta_e_11->Fill(theta_e_X);
	  h_theta_p_11->Fill(theta_p_X);
	  h_mom_11->Fill(p_p[j]);
	  h_theta_11->Fill(vX.Theta() * 180 /M_PI);
	  h_phi_11->Fill(vX.Phi() * 180 /M_PI);
	  h_mom_rho_11->Fill(vX.Pt());
	}
	else if(pPid[j]==0){
	  h_vtz_diff_0->Fill(z);
	  h_vtz_0->Fill(p_vtz[j]);
	  h_theta_e_0->Fill(theta_e_X);
	  h_theta_p_0->Fill(theta_p_X);
	  h_mom_0->Fill(p_p[j]);
	  h_theta_0->Fill(vX.Theta() * 180 /M_PI);
	  h_phi_0->Fill(vX.Phi() * 180 /M_PI);
	  h_mom_rho_0->Fill(vX.Pt());
	}
	else{
	  h_vtz_diff_E->Fill(z);
	  h_vtz_E->Fill(p_vtz[j]);
	  h_theta_e_E->Fill(theta_e_X);
	  h_theta_p_E->Fill(theta_p_X);
	  h_mom_E->Fill(p_p[j]);
	  h_theta_E->Fill(vX.Theta() * 180 /M_PI);
	  h_phi_E->Fill(vX.Phi() * 180 /M_PI);
	  h_mom_rho_E->Fill(vX.Pt());
	}
      }
    }

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
