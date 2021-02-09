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
#include "genpart.h"

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

  if( argc < 4 ){
    cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputRootFile] [outputPDFFile] [inputFiles]... \n\n";
    return -1;
  }
  
  //Define Variables
  bool goodneutron;
  int Runno, pMult, scinHits, nleadindex, nMult, eventnumber, genMult;
  double Ebeam, gated_charge, livetime, starttime, current;
  int pIndex[100];
  double pChi2pid[100], p_vtx[100], p_vty[100], p_vtz[100], p_p[100], theta_p[100], phi_p[100],pBeta[100], pPid[100], pCharge[100], hit_pindex[100], hit_detid[100];
  clashit * eHit = new clashit;
  TClonesArray* nHits = new TClonesArray("bandhit");
  TClonesArray * mcParts = new TClonesArray("genpart");
  char temp[100];

  ////////////////////////////////////////////////////
  //Start with defining input TChain
  ////////////////////////////////////////////////////
  TChain inChain("skim");
  for(int k = 3; k < argc; k++){
    inChain.Add(argv[k]);
  }
  // 	Run branches:
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
  inChain.SetBranchAddress("pIndex"		,pIndex);
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
  //    Simulation branches:
  inChain.SetBranchAddress("genMult"		,&genMult       );
  inChain.SetBranchAddress("mcParts"		,&mcParts       );

  ////////////////////////////////////////////////////
  //Now define output TFile and TTree
  ////////////////////////////////////////////////////
  TFile * outFile = new TFile(argv[1],"RECREATE");
  //Creating TH from tree
  vector<TH1*> hist_list_1;
  TH1D * h_mmiss_band_rec = new TH1D("mmiss_band_rec","Missing Mass of Reconstructed Particles (e,ep);m_miss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_band_rec);
  TH1D * h_mmiss_neut_rec = new TH1D("mmiss_neut_rec","Missing Mass of Reconstructed Particles (e,ep)n;m_miss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_neut_rec);


  TH1D * h_n_g_g = new TH1D("h_n_g_g","Momentum Resolution of Neutron; Missing Momentum Generated - Neutron Momentum Generated / Momentum;Events",100,-0.5,0.5);
  hist_list_1.push_back(h_n_g_g);
  TH1D * h_n_g_r = new TH1D("h_n_g_r","Momentum Resolution of Neutron; Missing Momentum Generated - Neutron Momentum Reconstructed / Momentum;Events",100,-0.5,0.5);
  hist_list_1.push_back(h_n_g_r);
  TH1D * h_n_r_g = new TH1D("h_n_r_g","Momentum Resolution of Neutron; Missing Momentum Reconstructed - Neutron Momentum Generated / Momentum;Events",100,-0.5,0.5);
  hist_list_1.push_back(h_n_r_g);
  TH1D * h_n_r_r = new TH1D("h_n_r_r","Momentum Resolution of Neutron; Missing Momentum Reconstructed - Neutron Momentum Reconstructed / Momentum;Events",100,-0.5,0.5);
  hist_list_1.push_back(h_n_r_r);


  //Begin Loop
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
    
    //Vertex Cuts
    if((eHit->getVtz() < -5) || (eHit->getVtz() > -1)){ continue; }
    

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
    //if(pMult==1){ continue; }
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

    genpart * e_Part = (genpart*)mcParts->At(0);
    genpart * p_Part = (genpart*)mcParts->At(1);
    genpart * n_Part = (genpart*)mcParts->At(2);

    TVector3 ve_g;
    ve_g.SetMagThetaPhi(e_Part->getMomentum(),e_Part->getTheta(),e_Part->getPhi());
    double Ee_g = ve_g.Mag();

    TVector3 vp_g;
    vp_g.SetMagThetaPhi(p_Part->getMomentum(),p_Part->getTheta(),p_Part->getPhi());
    double Ep_g = sqrt((mP * mP) + vp_g.Mag2());

    TVector3 vn_g;
    vn_g.SetMagThetaPhi(n_Part->getMomentum(),n_Part->getTheta(),n_Part->getPhi());
    double En_g = sqrt((mN * mN) + vn_g.Mag2());

    double Ebeam_2 = Ee_g + Ep_g + En_g - mD;
    cout<<Ebeam_2<<"\n";
    TVector3 vbeam(0,0,Ebeam_2);
    TVector3 vMiss_g = vbeam - ve_g - vp_g;

    double p_n_g = vn_g.Mag();
    double p_pmiss_g = vMiss_g.Mag();

    
    //Make Calculation for pMiss and mMiss
    TVector3 ve_r;
    ve_r.SetMagThetaPhi(eHit->getMomentum(),eHit->getTheta(),eHit->getPhi());
    double Ee_r = ve_r.Mag();

    TVector3 vp_r;
    vp_r.SetMagThetaPhi(p_p[aiP],theta_p[aiP],phi_p[aiP]);
    double Ep_r = sqrt((mP * mP) + vp_r.Mag2());

    TVector3 vMiss_r = vbeam - ve_r - vp_r;
    double eMiss_r = Ebeam_2 + mD - Ee_r - Ep_r;
    double p_pmiss_r = vMiss_r.Mag();
    double theta_pmiss_r = vMiss_r.Theta();
    double phi_pmiss_r = vMiss_r.Phi();
    double mmiss_r = sqrt((eMiss_r * eMiss_r) - vMiss_r.Mag2());

    ////////////////////////////////////////////////
    //(e,e'p)n Cuts
    ////////////////////////////////////////////////

    //High theta and pMiss cut
    //double aveVtz = (eHit->getVtz()+p_vtz[aiP]) / 2;
    if(p_pmiss_r<0.2){ continue; }
    if( theta_pmiss_r < (M_PI/2)){ continue; }
    if(!pointsToBand(theta_pmiss_r,phi_pmiss_r,p_vtz[aiP])){ continue; }    

    h_mmiss_band_rec->Fill(mmiss_r);

    if(!goodneutron){ continue; }
    bandhit * this_nHit = (bandhit*)nHits->At(nleadindex);
    double ToF = this_nHit->getTof();
    double ToFadc = this_nHit->getTofFadc();
    double DL = this_nHit->getDL().Mag();
    double ToM = ToF * 100 / DL;
    double ToMadc = ToFadc * 100 / DL;
    double Beta = this_nHit->getBeta();
    double EDep = this_nHit->getEdep() / 1100;
    TVector3 Mom = this_nHit->getMomentumN();
    double p_n_r= Mom.Mag();
    /*
    double thetamn=vmiss.Angle(Mom);
    double nX = this_nHit->getX();
    double nY = this_nHit->getY();
    double mX = this_nHit->getZ()*tan(theta_pmiss)*cos(phi_pmiss);
    double mY = this_nHit->getZ()*tan(theta_pmiss)*sin(phi_pmiss);
    */
    int sector = this_nHit->getSector();
    int barID = this_nHit->getBarID();

    if(mmiss_r < 0.85){ continue; }
    if(mmiss_r > 1.05){ continue; }
    if(EDep < 2){ continue; }
    if(ToM <  6){ continue; }
    if(ToM > 15){ continue; }

    h_mmiss_neut_rec->Fill(mmiss_r);

    h_n_g_g->Fill((p_pmiss_g-p_n_g)/p_pmiss_g);
    h_n_g_r->Fill((p_pmiss_g-p_n_r)/p_pmiss_g);
    h_n_r_g->Fill((p_pmiss_r-p_n_g)/p_pmiss_r);
    h_n_r_r->Fill((p_pmiss_r-p_n_r)/p_pmiss_r);


  }

  //Plot on pdf
  //Open
  char fileName[100];
  sprintf(fileName,"%s[",argv[2]);
  TCanvas * myCanvas = new TCanvas("My_Canvas","My_Canvas",1600,1000);
  myCanvas->SaveAs(fileName);
  //Write
  sprintf(fileName,"%s",argv[2]);
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Draw();
    myCanvas->SaveAs(fileName);

  }
  //Close
  sprintf(fileName,"%s]",argv[2]);
  TCanvas * endCanvas = new TCanvas("End_Canvas","End_Canvas",1600,1000);
  myCanvas->SaveAs(fileName);


  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
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
