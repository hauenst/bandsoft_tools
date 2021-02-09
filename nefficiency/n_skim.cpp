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
  for(int k = 2; k < argc; k++){
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
  //inChain.SetBranchAddress("mcParts"		,&mcParts       );
  //inChain.SetBranchAddress("genMult"		,&genMult       );

  ////////////////////////////////////////////////////
  //Now define output TFile and TTree
  ////////////////////////////////////////////////////
  int aiP;
  double p_pmiss, theta_pmiss, phi_pmiss, mmiss, eBind, p_neut, theta_neut, phi_neut, theta_nm;

  TFile * outFile = new TFile(argv[1],"RECREATE");
  TTree * outTree = new TTree("skim","Reduced Tree");

  // 	Run branches:
  outTree->Branch("Runno"       	,&Runno         ,"Runno/I"       );
  outTree->Branch("eventnumber"       	,&eventnumber   ,"eventnumber/I" );
  outTree->Branch("Ebeam"       	,&Ebeam         ,"Ebeam/D"       );
  outTree->Branch("gated_charge"	,&gated_charge  ,"gated_charge/D");
  outTree->Branch("livetime"    	,&livetime      ,"livetime/D"    );
  outTree->Branch("starttime"   	,&starttime     ,"starttime/D"   );
  outTree->Branch("current"     	,&current       ,"current/D"     );
  //	Electron branches:
  outTree->Branch("eHit"		,&eHit	        );
  //    Charged particle branches:
  outTree->Branch("pMult"   		,&pMult  	,"pMult/I"   	 );
  outTree->Branch("pIndex"		,pIndex		,"pIndex[100]/I"  	 );
  outTree->Branch("pPid"    		,pPid   	,"pPid[100]/D"    	 );
  outTree->Branch("pCharge" 		,pCharge	,"pCharge[100]/D" 	 );
  outTree->Branch("pChi2pid"		,pChi2pid	,"pChi2pid[100]/D"	 );
  outTree->Branch("p_vtx"   		,p_vtx  	,"p_vtx[100]/D"   	 );
  outTree->Branch("p_vty"   		,p_vty  	,"p_vty[100]/D"   	 );
  outTree->Branch("p_vtz"   		,p_vtz  	,"p_vtz[100]/D"   	 );
  outTree->Branch("p_p"     		,p_p    	,"p_p[100]/D"     	 );
  outTree->Branch("theta_p" 		,theta_p	,"theta_p[100]/D" 	 );
  outTree->Branch("phi_p"   		,phi_p  	,"phi_p[100]/D"   	 );
  outTree->Branch("pBeta"   		,pBeta  	,"pBeta[100]/D"   	 );
  //    Scintilator branches:
  outTree->Branch("scinHits"  		,&scinHits  	,"scinHits/I"    );
  outTree->Branch("hit_pindex"		,hit_pindex  	,"hit_pindex[100]/D"  );
  outTree->Branch("hit_detid"		,hit_detid  	,"hit_detid[100]/D"   );
  //    Neutron branches:
  outTree->Branch("nHits"		,&nHits	        );
  outTree->Branch("nMult"		,&nMult	        ,"nMult/I");
  outTree->Branch("goodneutron"         ,&goodneutron   ,"goodneutron/O" );
  outTree->Branch("nleadindex"          ,&nleadindex    ,"nleadindex/I"  );
  //   Missing momentum branches in (e,e'p)n
  outTree->Branch("aiP"			,&aiP		,"aiP/I"	 );
  outTree->Branch("p_pmiss"    		,&p_pmiss      	,"p_pmiss/D"	 );
  outTree->Branch("theta_pmiss"    	,&theta_pmiss   ,"theta_pmiss/D" );
  outTree->Branch("phi_pmiss"    	,&phi_pmiss    	,"phi_pmiss/D"	 );
  outTree->Branch("mmiss" 	   	,&mmiss    	,"mmiss/D"	 );
  outTree->Branch("eBind" 	   	,&eBind    	,"eBind/D"	 );
  //outTree->Branch("p_neut"    		,&p_neut      	,"p_neut/D"	 );
  //outTree->Branch("theta_neut"    	,&theta_neut    ,"theta_neut/D"  );
  //outTree->Branch("phi_neut"    	,&phi_neut    	,"phi_neut/D"	 );
  //outTree->Branch("theta_nm"    	,&theta_nm      ,"theta_nm/D"  );


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
    aiP; //array index of the Proton
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

    TVector3 vmiss = vbeam - ve - vp;
    double emiss = Ebeam + mD - Ee - Ep;
    double e_rec = sqrt((mN*mN) + vmiss.Mag2());
    eBind = mP + mN + Ebeam - Ep - Ee - e_rec;
    p_pmiss = vMiss.Mag();
    theta_pmiss = vMiss.Theta();
    phi_pmiss = vMiss.Phi();
    mmiss = sqrt((eMiss * eMiss) - vMiss.Mag2());

    ////////////////////////////////////////////////
    //(e,e'p)n Cuts
    ////////////////////////////////////////////////
    //High theta and pMiss cut
    //double aveVtz = (eHit->getVtz()+p_vtz[aiP]) / 2;
    if(p_pmiss<0.2){ continue; }
    //if(mmiss<0.85){ continue; }
    //if(mmiss>1.05){ continue; }
    if( theta_pmiss < (M_PI/2)){ continue; }
    if(!pointsToBand(theta_pmiss,phi_pmiss,p_vtz[aiP])){ continue; }
    
    //Now fill the new tree
    outTree->Fill();
  }

  outFile->cd();
  outTree->Write();
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
