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

bool lowThetaCut(double theta, double chi2PID, double vtzDiff){
  
  if(theta > (50 * M_PI / 180)){
    return false;
  }
  if(abs(chi2PID-0.459179)>(3*1.2085)){
    return false;
  }
  if(abs(vtzDiff-0.484268)>(3*1.30286)){
    return false;
  }
  
  return true;
}

bool highThetaCut(double theta, double chi2PID, double vtzDiff){
  
  if(theta < (30 * M_PI / 180)){
    return false;
  }
  if(abs(chi2PID-0.821486)>(3*2.10464)){
    return false;
  }
  if(abs(vtzDiff+0.0163482)>(3*0.781676)){
    return false;
  }
  
  return true;
}


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
  TH1D * h_e_vtx = new TH1D("e_vtx","X Vertex of e;vtx;Events",100,-10,10);
  hist_list_1.push_back(h_e_vtx);
  TH1D * h_e_vty = new TH1D("e_vty","Y Vertex of e;vty;Events",100,-10,10);
  hist_list_1.push_back(h_e_vty);
  TH1D * h_e_vtz = new TH1D("e_vtz","Z Vertex of e;vtz;Events",100,-10,10);
  hist_list_1.push_back(h_e_vtz);
  TH1D * h_p_chi = new TH1D("p_chi","Chi2 pid of p; Chi2;Events",100,-10,10);
  hist_list_1.push_back(h_p_chi);
  TH1D * h_p_theta = new TH1D("p_theta","theta of p;theta;Events",155,0,155);
  hist_list_1.push_back(h_p_theta);
  TH1D * h_mmiss = new TH1D("mmiss","Missing Mass;mmiss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss);
  TH1D * h_mmiss_band = new TH1D("mmiss_band","Missing Mass;mmiss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_band);
  TH1D * h_pmiss = new TH1D("pmiss","Missing Momentum;pmiss;Events",100,0,1);
  hist_list_1.push_back(h_pmiss);
  TH1D * h_theta_mn = new TH1D("theta_mn","Theta Between miss and neutron;theta_mn;Events",90,0,45);
  hist_list_1.push_back(h_theta_mn);
  TH1D * h_Mom[6];
  for(int i = 0; i<6; i++){
    sprintf(temp,"Mom_sec%d",i);
    h_Mom[i] = new TH1D(temp,"p;p;Events",250,0,4.5);
    hist_list_1.push_back(h_Mom[i]);
  } 

  vector<TH2*> hist_list_2;
  TH2D * h_e_vtz_p_vtz = new TH2D("e_vtz_p_vtz","Z Vertex of e vs Z Vertex of p;e_vtx;p_vtx;Events",100,-10,10,100,-10,10);
  hist_list_2.push_back(h_e_vtz_p_vtz);
  TH2D * h_p_beta_p_mom = new TH2D("p_beta_p_mom","Beta of p vs Momentum of p;Beta;Momentum;Events",100,0,1,100,0,4.5);
  hist_list_2.push_back(h_p_beta_p_mom);
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","mmiss vs pmiss;mmiss;pmiss;Events",100,0,2,100,0,1);
  hist_list_2.push_back(h_mmiss_pmiss);
  TH2D * h_mmiss_deltaphi = new TH2D("mmiss_deltaphi","mmiss vs deltaphi;mmiss;deltaphi;Events",100,0,2,100,90,270);
  hist_list_2.push_back(h_mmiss_deltaphi);
  TH2D * h_mmiss_p_theta = new TH2D("mmiss_p_theta","mmiss vs p_theta;mmiss;p_theta;Events",100,0,2,100,0,100);
  hist_list_2.push_back(h_mmiss_p_theta);
  TH2D * h_mmiss_pmiss_theta = new TH2D("mmiss_pmiss_theta","mmiss vs p_theta;mmiss;p_theta;Events",100,0,2,180,0,180);
  hist_list_2.push_back(h_mmiss_pmiss_theta);
  TH2D * h_mmiss_deltaphi_band = new TH2D("mmiss_deltaphi_band","mmiss vs deltaphi;mmiss;deltaphi;Events",100,0,2,100,90,270);
  hist_list_2.push_back(h_mmiss_deltaphi_band);
  TH2D * h_mmiss_p_theta_band = new TH2D("mmiss_p_theta_band","mmiss vs p_theta;mmiss;p_theta;Events",100,0,2,100,0,100);
  hist_list_2.push_back(h_mmiss_p_theta_band);
  TH2D * h_mmiss_pmiss_theta_band = new TH2D("mmiss_pmiss_theta_band","mmiss vs p_theta;mmiss;p_theta;Events",100,0,2,180,0,180);
  hist_list_2.push_back(h_mmiss_pmiss_theta_band);
  TH2D * h_mom_m_mom_n = new TH2D("mom_m_mom_n","Missing Momentum vs Neutron Momentum;mom_m;mom_n;Events",40,0,0.8,40,0,0.8);
  hist_list_2.push_back(h_mom_m_mom_n);
  TH2D * h_theta_m_theta_n = new TH2D("theta_m_theta_n","Missing Theta vs Neutron Theta;theta_m;theta_n;Events",40,150,180,40,150,180);
  hist_list_2.push_back(h_theta_m_theta_n);
  TH2D * h_phi_m_phi_n = new TH2D("phi_m_phi_n","Missing Phi vs Neutron Phi;phi_m;phi_n;Events",40,-180,180,40,-180,180);
  hist_list_2.push_back(h_phi_m_phi_n);
  TH2D * h_Mom_SF[6];
  for(int i = 0; i<6; i++){
    sprintf(temp,"Mom_SF_sec%d",i);
    h_Mom_SF[i] = new TH2D(temp,"p vs SF;p;SF;Events",250,0,4.5,100,0.0,0.4);
    hist_list_2.push_back(h_Mom_SF[i]);
  } 

  int fin = inChain.GetEntries();
  for(int i = 0; i < fin; i++){
    eHit->Clear();
    nHits->Clear();
    inChain.GetEntry(i);

    //Display completed  
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }
    

    ////////////////////////////////////////////////
    //(e,e') Cuts
    ////////////////////////////////////////////////

    //Electron Cuts
    if(eHit->getV() < 14){ continue; }
    if(eHit->getW() < 14){ continue; }
    if(eHit->getEoP() < 0.18){ continue; }
    if(eHit->getEoP() > 0.28){ continue; }
    if(eHit->getMomentum() < 1){ continue; }
    if(eHit->getMomentum() > Ebeam){ continue; }
    h_Mom_SF[eHit->getSector() - 1]->Fill(eHit->getMomentum(),eHit->getEoP());

    //Vertex Cuts
    if((eHit->getVtz() < -5) || (eHit->getVtz() > -1)){ continue; }
    h_e_vtx->Fill(eHit->getVtx());
    h_e_vty->Fill(eHit->getVty());
    h_e_vtz->Fill(eHit->getVtz());

    ////////////////////////////////////////////////
    //(e,e'p) Cuts
    ////////////////////////////////////////////////
    //Find protons in particle branches
    bool isProton = false;
    int numProton = 0;
    int pindexProton;
    int aiP; //array index of the Proton
    for(int j = 0; j < pMult; j++){
      if((abs(pPid[j]-2212)<0.001) && (abs(pCharge[j]-1)<0.001)){
	aiP = j;
	pindexProton=pIndex[j];
	numProton++;
      }
    }
    if(numProton!=1){ continue; }

    //Check that the proton is a FD proton
    bool FToF = false;
    for(int k = 0; k < scinHits; k++){
      if(abs((double)pindexProton - hit_pindex[k])<0.01){
	if(abs(hit_detid[k]-12)<0.01){
	  FToF=true;
	}
      }
    }
    if(!FToF){continue;}
    if(!lowThetaCut(theta_p[aiP],pChi2pid[aiP],abs(eHit->getVtz()-p_vtz[aiP]))){continue;}


    //Make Calculation for pMiss and mMiss
    TVector3 vbeam(0,0,Ebeam);
    TVector3 ve;
    ve.SetMagThetaPhi(eHit->getMomentum(),eHit->getTheta(),eHit->getPhi());
    double Ee = ve.Mag();

    TVector3 vp;
    vp.SetMagThetaPhi(p_p[aiP],theta_p[aiP],phi_p[aiP]);
    double Ep = sqrt((mP * mP) + vp.Mag2());

    TVector3 vMiss = vbeam - ve - vp;
    double eMiss = Ebeam + mD - Ee - Ep;
    double pMiss = vMiss.Mag();
    double mMiss = sqrt((eMiss * eMiss) - vMiss.Mag2());

    h_p_chi->Fill(pChi2pid[aiP]);
    h_p_theta->Fill(theta_p[aiP] * 180 / M_PI);
    h_mmiss->Fill(mMiss);    
    h_Mom[eHit->getSector() - 1]->Fill(eHit->getMomentum());
    h_e_vtz_p_vtz->Fill(eHit->getVtz(),p_vtz[aiP]);
    h_p_beta_p_mom->Fill(pBeta[aiP],p_p[aiP]);
    h_mmiss_pmiss->Fill(mMiss,pMiss);    
    h_mmiss_deltaphi->Fill(mMiss,abs(eHit->getPhi() - phi_p[aiP]) * 180 / M_PI);
    h_mmiss_p_theta->Fill(mMiss,theta_p[aiP] * 180 / M_PI);
    h_mmiss_pmiss_theta->Fill(mMiss,vMiss.Theta() * 180 / M_PI);

    ////////////////////////////////////////////////
    //(e,e'p)n Cuts
    ////////////////////////////////////////////////
    //High theta and pMiss cut
    //double aveVtz = (eHit->getVtz()+p_vtz[aiP]) / 2;
    if(pMiss<0.25){ continue; }
    if( vMiss.Theta() < (M_PI/2)){ continue; }
    if(!pointsToBand(vMiss.Theta(),vMiss.Phi(),p_vtz[aiP])){ continue; }

    h_mmiss_pmiss_theta_band->Fill(mMiss,vMiss.Theta() * 180 / M_PI);
    h_mmiss_band->Fill(mMiss);
    h_mmiss_deltaphi_band->Fill(mMiss,abs(eHit->getPhi() - phi_p[aiP]) * 180 / M_PI);
    h_mmiss_p_theta_band->Fill(mMiss,theta_p[aiP] * 180 / M_PI);
    h_pmiss->Fill(pMiss);

    ////////////////////////////////////////////////
    //(e,e'pn) Cuts
    ////////////////////////////////////////////////
    if(!goodneutron){ continue; }
    bandhit * this_nHit = (bandhit*)nHits->At(nleadindex);
    double ToF = this_nHit->getTof();
    double DL = this_nHit->getDL().Mag();
    double ToM = ToF * 100 / DL;
    double Beta = this_nHit->getBeta();
    double EDep = this_nHit->getEdep() / 1100;
    TVector3 Mom = this_nHit->getMomentumN();
    if(EDep < 2){ continue; }
    if(ToM <  6){ continue; }
    if(ToM > 15){ continue; }
    if(mMiss < 0.85){ continue; }
    if(mMiss > 1.05){ continue; }
    //if( Mom.Theta() < (M_PI/2)){ continue; }
    //if(!pointsToBand(Mom.Theta(),Mom.Phi(),aveVtz)){ continue; }

    h_mom_m_mom_n->Fill(pMiss,Mom.Mag());
    h_theta_m_theta_n->Fill(vMiss.Theta()*180/M_PI,Mom.Theta()*180/M_PI);
    h_phi_m_phi_n->Fill(vMiss.Phi()*180/M_PI,Mom.Phi()*180/M_PI);
    h_theta_mn->Fill(vMiss.Angle(Mom));

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
