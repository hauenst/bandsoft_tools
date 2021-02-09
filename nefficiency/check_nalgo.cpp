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
    cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [out text file] [inputFiles]... \n\n";
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

  std::ofstream textfile(argv[2]);
  
  TChain inChain("skim");
 //Creat input tree
  for(int k = 3; k < argc; k++){
    inChain.Add(argv[k]);
  }
  //TFile * inFile = new TFile(argv[2]);
  //TTree * inTree = (TTree*)inFile->Get("skim");
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
  TH1D * h_mmiss_band = new TH1D("mmiss_band","Missing Mass;mmiss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_band);
  TH1D * h_mmiss_neut = new TH1D("mmiss_neut","Missing Mass;mmiss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_neut);
  TH1D * h_theta_mn = new TH1D("theta_mn","Theta Between miss and neutron;theta_mn;Events",90,0,45);
  hist_list_1.push_back(h_theta_mn);
  TH1I * h_nMult = new TH1I("h_nMult","nMult;nMult;Events",5,0,5);
  hist_list_1.push_back(h_nMult);
  TH1D * h_nres = new TH1D("h_nres","nres;nres;Events",20,-2,2);
  hist_list_1.push_back(h_nres);

  vector<TH2*> hist_list_2;
  TH2D * h_mom_m_mom_n = new TH2D("mom_m_mom_n","Missing Momentum vs Neutron Momentum;mom_m;mom_n;Events",40,0,0.8,40,0,0.8);
  hist_list_2.push_back(h_mom_m_mom_n);
  TH2D * h_theta_m_theta_n = new TH2D("theta_m_theta_n","Missing Theta vs Neutron Theta;theta_m;theta_n;Events",40,150,180,40,150,180);
  hist_list_2.push_back(h_theta_m_theta_n);
  TH2D * h_phi_m_phi_n = new TH2D("phi_m_phi_n","Missing Phi vs Neutron Phi;phi_m;phi_n;Events",40,-180,180,40,-180,180);
  hist_list_2.push_back(h_phi_m_phi_n);
  TH2D * h_x_m_x_n = new TH2D("x_m_x_n","Missing X vs Neutron X;x_m;x_n;Events",200,-200,200,200,-200,200);
  hist_list_2.push_back(h_x_m_x_n);
  TH2D * h_y_m_y_n = new TH2D("y_m_y_n","Missing Y vs Neutron Y;y_m;y_n;Events",200,-50,150,200,-50,150);
  hist_list_2.push_back(h_y_m_y_n);
  TH2D * h_ToM_Edep = new TH2D("ToM_Edep","ToM vs Edep;ToM;Edep;Events",200,-50,100,200,0,10);
  hist_list_2.push_back(h_ToM_Edep);
  TH2D * h_tpm_mm = new TH2D("tpm_mm","Theta Pmiss vs Mmiss;Theta;MMiss;Events",180,0,180,200,0,5);
  hist_list_2.push_back(h_tpm_mm);


  TH2D * h_theta_mn_x_m = new TH2D("theta_mn_x_m",";theta_mn;x_m;Events",90,0,45,200,-150,150);
  hist_list_2.push_back(h_theta_mn_x_m);
  TH2D * h_theta_mn_y_m = new TH2D("theta_mn_y_m",";theta_mn;y_m;Events",90,0,45,200,-50,150);
  hist_list_2.push_back(h_theta_mn_y_m);
  TH2D * h_theta_mn_mom_m = new TH2D("theta_mn_mom_m",";theta_mn;mom_m;Events",90,0,45,40,0,0.8);
  hist_list_2.push_back(h_theta_mn_mom_m);
  TH2D * h_theta_mn_theta_m = new TH2D("theta_mn_theta_m",";theta_mn;theta_m;Events",90,0,45,40,150,180);
  hist_list_2.push_back(h_theta_mn_theta_m);
  TH2D * h_theta_mn_phi_m = new TH2D("theta_mn_phi_m",";theta_mn;phi_m;Events",90,0,45,40,-180,180);
  hist_list_2.push_back(h_theta_mn_phi_m);

  TH2D * h_theta_mn_x_n = new TH2D("theta_mn_x_n",";theta_mn;x_n;Events",90,0,45,200,-150,150);
  hist_list_2.push_back(h_theta_mn_x_n);
  TH2D * h_theta_mn_y_n = new TH2D("theta_mn_y_n",";theta_mn;y_n;Events",90,0,45,200,-50,150);
  hist_list_2.push_back(h_theta_mn_y_n);
  TH2D * h_theta_mn_mom_n = new TH2D("theta_mn_mom_n",";theta_mn;mom_n;Events",90,0,45,40,0,0.8);
  hist_list_2.push_back(h_theta_mn_mom_n);
  TH2D * h_theta_mn_theta_n = new TH2D("theta_mn_theta_n",";theta_mn;theta_n;Events",90,0,45,40,150,180);
  hist_list_2.push_back(h_theta_mn_theta_n);
  TH2D * h_theta_mn_phi_n = new TH2D("theta_mn_phi_n",";theta_mn;phi_n;Events",90,0,45,40,-180,180);
  hist_list_2.push_back(h_theta_mn_phi_n);

  TH2D * h_theta_mn_EDep = new TH2D("theta_mn_EDep",";theta_mn;EDep;Events",90,0,45,40,0,10);
  hist_list_2.push_back(h_theta_mn_EDep);
  TH2D * h_theta_mn_ToM = new TH2D("theta_mn_ToM",";theta_mn;ToM;Events",90,0,45,40,0,30);
  hist_list_2.push_back(h_theta_mn_ToM);
  TH2D * h_theta_mn_ToMadc = new TH2D("theta_mn_ToMadc",";theta_mn;ToMadc;Events",90,0,45,40,0,30);
  hist_list_2.push_back(h_theta_mn_ToMadc);
  int fin = inChain.GetEntries();
  for(int i = 0; i < fin; i++){
    eHit->Clear();
    nHits->Clear();
    inChain.GetEntry(i);
    if(pMult!=1){continue;}
    if(mmiss<0.85){ continue; }
    if(mmiss>1.05){ continue; }

    textfile << Runno <<"\t"<< eventnumber<<"\n";

    
    //Make Calculation for pMiss and mMiss
    TVector3 vbeam(0,0,Ebeam);
    TVector3 ve;
    ve.SetMagThetaPhi(eHit->getMomentum(),eHit->getTheta(),eHit->getPhi());
    TVector3 vp;
    vp.SetMagThetaPhi(p_p[aiP],theta_p[aiP],phi_p[aiP]);
    TVector3 vmiss;
    vmiss.SetMagThetaPhi(p_pmiss,theta_pmiss,phi_pmiss);

    h_mmiss_band->Fill(mmiss);
    

    ////////////////////////////////////////////////
    //(e,e'pn) Cuts
    ////////////////////////////////////////////////
    
    int nMult=nHits->GetEntries();
    h_nMult->Fill(nMult);
    if(nMult==1){
      bandhit * this_nHit = (bandhit*)nHits->At(0);
      double ToF = this_nHit->getTof();
      double ToFadc = this_nHit->getTofFadc();
      double DL = this_nHit->getDL().Mag();
      double ToM = ToF * 100 / DL;
      double ToMadc = ToFadc * 100 / DL;
      double Beta = this_nHit->getBeta();
      double EDep = this_nHit->getEdep() / 1100;
      TVector3 Mom = this_nHit->getMomentumN();
      double thetamn=vmiss.Angle(Mom)*180/M_PI;
      h_ToM_Edep->Fill(ToM,EDep);
      h_theta_mn_ToM->Fill(thetamn,ToM);
      h_theta_mn_ToMadc->Fill(thetamn,ToMadc);
      if((ToM >  6) && (ToM < 15)){
	h_theta_mn_EDep->Fill(thetamn,EDep);
      }
    }
    //if(nMult!=1){continue;}
    if(!goodneutron){ continue; }
    int k=0;
    //if(k<0){ continue; }
    bandhit * this_nHit = (bandhit*)nHits->At(nleadindex);
    double ToF = this_nHit->getTof();
    double ToFadc = this_nHit->getTofFadc();
    double DL = this_nHit->getDL().Mag();
    double ToM = ToF * 100 / DL;
    double ToMadc = ToFadc * 100 / DL;
    double Beta = this_nHit->getBeta();
    double EDep = this_nHit->getEdep() / 1100;
    TVector3 Mom = this_nHit->getMomentumN();
    //if(ToM > 15){ continue; }
    //if(thetamn>15){continue;}
    double thetamn=vmiss.Angle(Mom)*180/M_PI;
    double nX = this_nHit->getX();
    double nY = this_nHit->getY();
    double mX = this_nHit->getZ()*tan(theta_pmiss)*cos(phi_pmiss);
    double mY = this_nHit->getZ()*tan(theta_pmiss)*sin(phi_pmiss);
    int sector = this_nHit->getSector();
    int barID = this_nHit->getBarID();
    //if(sector!=3){continue;}
    //if((nY<-25) || (nY>25)){continue;}
    //if(EDep < 2){ continue; }
    //if(ToM <  6){ continue; }
    //if(ToM > 15){ continue; }
    h_theta_mn->Fill(thetamn);
    h_mmiss_neut->Fill(mmiss);
    h_mom_m_mom_n->Fill(p_pmiss,Mom.Mag());
    h_theta_m_theta_n->Fill(theta_pmiss*180/M_PI,Mom.Theta()*180/M_PI);
    h_phi_m_phi_n->Fill(phi_pmiss*180/M_PI,Mom.Phi()*180/M_PI);
    h_x_m_x_n->Fill(mX,nX);
    h_y_m_y_n->Fill(mY,nY);
    h_nres->Fill((p_pmiss-Mom.Mag())/p_pmiss);
    h_theta_mn_x_m->Fill(thetamn,mX);
    h_theta_mn_y_m->Fill(thetamn,mY);
    h_theta_mn_mom_m->Fill(thetamn,p_pmiss);
    h_theta_mn_theta_m->Fill(thetamn,theta_pmiss*180/M_PI);
    h_theta_mn_phi_m->Fill(thetamn,phi_pmiss*180/M_PI);
    h_theta_mn_x_n->Fill(thetamn,nX);
    h_theta_mn_y_n->Fill(thetamn,nY);
    h_theta_mn_mom_n->Fill(thetamn,Mom.Mag());
    h_theta_mn_theta_n->Fill(thetamn,Mom.Theta()*180/M_PI);
    h_theta_mn_phi_n->Fill(thetamn,Mom.Phi()*180/M_PI);
 

  }
  //inFile->Close();
  textfile.close();
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
