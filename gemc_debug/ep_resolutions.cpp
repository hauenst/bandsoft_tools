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
  TH1D * h_mmiss_band_gen = new TH1D("mmiss_band_gen","Missing Mass of Generated Particles;m_miss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_band_gen);
  TH1D * h_mmiss_band_rec = new TH1D("mmiss_band_rec","Missing Mass of Reconstructed Particles;m_miss;Events",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_band_rec);
  //TH1D * h_mmiss_neut = new TH1D("mmiss_neut","Missing Mass;mmiss;Events",100,0.4,1.4);
  //hist_list_1.push_back(h_mmiss_neut);
  TH1D * h_e_mom_res = new TH1D("h_e_mom_res","Momentum Resolution of Electron;Delta Mom / Mom;Events",100,-0.05,0.05);
  hist_list_1.push_back(h_e_mom_res);
  TH1D * h_e_theta_res = new TH1D("h_e_theta_res","Theta Resolution of Electron;Delta Theta;Events",100,-0.5,0.5);
  hist_list_1.push_back(h_e_theta_res);
  TH1D * h_e_phi_res = new TH1D("h_e_phi_res","Phi Resolution of Electron;Delta Phi;Events",100,-3,3);
  hist_list_1.push_back(h_e_phi_res);
  TH1D * h_p_mom_res = new TH1D("h_p_mom_res","Momentum Resolution of Proton;Delta Mom / Mom;Events",100,-0.05,0.05);
  hist_list_1.push_back(h_p_mom_res);
  TH1D * h_p_theta_res = new TH1D("h_p_theta_res","Theta Resolution of Proton;Delta Theta;Events",100,-0.5,0.5);
  hist_list_1.push_back(h_p_theta_res);
  TH1D * h_p_phi_res = new TH1D("h_p_phi_res","Phi Resolution of Proton;Delta Phi;Events",100,-3,3);
  hist_list_1.push_back(h_p_phi_res);

  /*
  TH1D * h_n_mom_res = new TH1D("h_n_mom_res","n_mom_res;res;Events",100,-0.25,0.25);
  hist_list_1.push_back(h_n_mom_res);
  TH1D * h_n_theta_res = new TH1D("h_n_theta_res","n_theta_res;res;Events",100,-0.25,0.25);
  hist_list_1.push_back(h_n_theta_res);
  TH1D * h_n_phi_res = new TH1D("h_n_phi_res","n_phi_res;res;Events",100,-0.25,0.25);
  hist_list_1.push_back(h_n_phi_res);
  */

  vector<TH2*> hist_list_2;
  TH2D * h_e_mom_diff_mom = new TH2D("e_mom_diff_mom","Electron;Mom;Delta Mom / Mom;Events",100,2,4.5,100,-0.05,0.05);
  hist_list_2.push_back(h_e_mom_diff_mom);
  TH2D * h_e_mom_diff_theta = new TH2D("e_mom_diff_theta","Electron;Mom;Delta Theta;Events",100,2,4.5,100,-0.5,0.5);
  hist_list_2.push_back(h_e_mom_diff_theta);
  TH2D * h_e_mom_diff_phi = new TH2D("e_mom_diff_phi","Electron;Mom;Delta Phi;Events",100,2,4.5,100,-3,3);
  hist_list_2.push_back(h_e_mom_diff_phi);
  TH2D * h_e_theta_diff_mom = new TH2D("e_theta_diff_mom","Electron;Theta;Delta Mom / Mom;Events",80,10,40,100,-0.05,0.05);
  hist_list_2.push_back(h_e_theta_diff_mom);
  TH2D * h_e_theta_diff_theta = new TH2D("e_theta_diff_theta","Electron;Theta;Delta Theta;Events",80,10,40,100,-0.5,0.5);
  hist_list_2.push_back(h_e_theta_diff_theta);
  TH2D * h_e_theta_diff_phi = new TH2D("e_theta_diff_phi","Electron;Theta;Delta Phi;Events",80,10,40,100,-3,3);
  hist_list_2.push_back(h_e_theta_diff_phi);
  TH2D * h_e_phi_diff_mom = new TH2D("e_phi_diff_mom","Electron;Phi;Delta Mom / Mom;Events",360,-180,180,100,-0.05,0.05);
  hist_list_2.push_back(h_e_phi_diff_mom);
  TH2D * h_e_phi_diff_theta = new TH2D("e_phi_diff_theta","Electron;Phi;Delta Theta;Events",360,-180,180,100,-0.5,0.5);
  hist_list_2.push_back(h_e_phi_diff_theta);
  TH2D * h_e_phi_diff_phi = new TH2D("e_phi_diff_phi","Electron;Phi;Delta Phi;Events",360,-180,180,100,-3,3);
  hist_list_2.push_back(h_e_phi_diff_phi);

  TH2D * h_p_mom_diff_mom = new TH2D("p_mom_diff_mom","Proton;Mom;Delta Mom / Mom;Events",100,0,3.5,100,-0.05,0.05);
  hist_list_2.push_back(h_p_mom_diff_mom);
  TH2D * h_p_mom_diff_theta = new TH2D("p_mom_diff_theta","Proton;Mom;Delta Theta;Events",100,0,3.5,100,-0.5,0.5);
  hist_list_2.push_back(h_p_mom_diff_theta);
  TH2D * h_p_mom_diff_phi = new TH2D("p_mom_diff_phi","Proton;Mom;Delta Phi;Events",100,0,3.5,100,-3,3);
  hist_list_2.push_back(h_p_mom_diff_phi);
  TH2D * h_p_theta_diff_mom = new TH2D("p_theta_diff_mom","Proton;Theta;Delta Mom / Mom;Events",80,10,45,100,-0.05,0.05);
  hist_list_2.push_back(h_p_theta_diff_mom);
  TH2D * h_p_theta_diff_theta = new TH2D("p_theta_diff_theta","Proton;Theta;Delta Theta;Events",80,10,45,100,-0.5,0.5);
  hist_list_2.push_back(h_p_theta_diff_theta);
  TH2D * h_p_theta_diff_phi = new TH2D("p_theta_diff_phi","Proton;Theta;Delta Phi;Events",80,10,45,100,-3,3);
  hist_list_2.push_back(h_p_theta_diff_phi);
  TH2D * h_p_phi_diff_mom = new TH2D("p_phi_diff_mom","Proton;Phi;Delta Mom / Mom;Events",360,-180,180,100,-0.05,0.05);
  hist_list_2.push_back(h_p_phi_diff_mom);
  TH2D * h_p_phi_diff_theta = new TH2D("p_phi_diff_theta","Proton;Phi;Delta Theta;Events",360,-180,180,100,-0.5,0.5);
  hist_list_2.push_back(h_p_phi_diff_theta);
  TH2D * h_p_phi_diff_phi = new TH2D("p_phi_diff_phi","Proton;Phi;Delta Phi;Events",360,-180,180,100,-3,3);
  hist_list_2.push_back(h_p_phi_diff_phi);

  /*
  TH2D * h_e_theta = new TH2D("e_theta","e_theta;e_theta;theta;Events",180,0,180,180,0,180);
  hist_list_2.push_back(h_e_theta);
  TH2D * h_e_phi = new TH2D("e_phi","e_phi;phi;Events",360,-180,180,360,-180,180);
  hist_list_2.push_back(h_e_phi);
  TH2D * h_p_mom = new TH2D("p_mom","p_mom;mom;Events",100,0,4,100,0,4);
  hist_list_2.push_back(h_p_mom);
  TH2D * h_p_theta = new TH2D("p_theta","p_theta;p_theta;theta;Events",180,0,180,180,0,180);
  hist_list_2.push_back(h_p_theta);
  TH2D * h_p_phi = new TH2D("p_phi","p_phi;phi;Events",360,-180,180,360,-180,180);
  hist_list_2.push_back(h_p_phi);
  */

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

    
    //Make Calculation for pMiss and mMiss
    TVector3 vbeam(0,0,Ebeam_2);
    TVector3 ve;
    ve.SetMagThetaPhi(eHit->getMomentum(),eHit->getTheta(),eHit->getPhi());
    double Ee = ve.Mag();

    TVector3 vp;
    vp.SetMagThetaPhi(p_p[aiP],theta_p[aiP],phi_p[aiP]);
    double Ep = sqrt((mP * mP) + vp.Mag2());

    TVector3 vMiss = vbeam - ve - vp;
    double eMiss = Ebeam_2 + mD - Ee - Ep;
    double p_pmiss = vMiss.Mag();
    double theta_pmiss = vMiss.Theta();
    double phi_pmiss = vMiss.Phi();
    double mmiss = sqrt((eMiss * eMiss) - vMiss.Mag2());

    ////////////////////////////////////////////////
    //(e,e'p)n Cuts
    ////////////////////////////////////////////////

    //High theta and pMiss cut
    //double aveVtz = (eHit->getVtz()+p_vtz[aiP]) / 2;
    //if(p_pmiss<0.2){ continue; }
    //if(mmiss<0.85){ continue; }
    //if(mmiss>1.05){ continue; }
    //if( theta_pmiss < (M_PI/2)){ continue; }
    //if(!pointsToBand(theta_pmiss,phi_pmiss,p_vtz[aiP])){ continue; }

    
    TVector3 vMiss_g = vbeam - ve_g - vp_g;
    double eMiss_g = Ebeam_2 + mD - Ee_g - Ep_g;
    double p_pmiss_g = vMiss_g.Mag();
    double theta_pmiss_g = vMiss_g.Theta();
    double phi_pmiss_g = vMiss_g.Phi();
    double mmiss_g = sqrt((eMiss_g * eMiss_g) - vMiss_g.Mag2());
    

    double eg_m = e_Part->getMomentum();
    double eg_t = e_Part->getTheta()*180/M_PI;
    double eg_p = e_Part->getPhi()*180/M_PI;

    double er_m = eHit->getMomentum();    
    double er_t = eHit->getTheta()*180/M_PI;
    double er_p = eHit->getPhi()*180/M_PI;

    double pg_m = p_Part->getMomentum();
    double pg_t = p_Part->getTheta()*180/M_PI;
    double pg_p = p_Part->getPhi()*180/M_PI;

    double pr_m = p_p[aiP];
    double pr_t = theta_p[aiP]*180/M_PI;
    double pr_p = phi_p[aiP]*180/M_PI;
    /*
    h_e_phi_gen->Fill(e_Part->getPhi()*180/M_PI);
    h_e_phi_rec->Fill(eHit->getPhi()*180/M_PI);
    h_p_phi_gen->Fill(p_Part->getPhi()*180/M_PI);
    h_p_phi_rec->Fill(phi_p[aiP]*180/M_PI);
x    */

    h_mmiss_band_gen->Fill(mmiss_g);
    h_mmiss_band_rec->Fill(mmiss);

    h_e_mom_res->Fill(  (eg_m-er_m)/eg_m);
    h_e_theta_res->Fill((eg_t-er_t));
    h_e_phi_res->Fill(  (eg_p-er_p));

    h_p_mom_res->Fill(  (pg_m-pr_m)/pg_m);
    h_p_theta_res->Fill((pg_t-pr_t));
    h_p_phi_res->Fill(  (pg_p-pr_p));

  h_e_mom_diff_mom->Fill(eg_m    ,(eg_m-er_m)/eg_m);
  h_e_mom_diff_theta->Fill(eg_m  ,(eg_t-er_t));
  h_e_mom_diff_phi->Fill(eg_m    ,(eg_p-er_p));
  h_e_theta_diff_mom->Fill(eg_t  ,(eg_m-er_m)/eg_m);
  h_e_theta_diff_theta->Fill(eg_t,(eg_t-er_t));
  h_e_theta_diff_phi->Fill(eg_t  ,(eg_p-er_p));
  h_e_phi_diff_mom->Fill(eg_p    ,(eg_m-er_m)/eg_m);
  h_e_phi_diff_theta->Fill(eg_p  ,(eg_t-er_t));
  h_e_phi_diff_phi->Fill(eg_p    ,(eg_p-er_p));
    
  h_p_mom_diff_mom->Fill(pg_m    ,(pg_m-pr_m)/pg_m);
  h_p_mom_diff_theta->Fill(pg_m  ,(pg_t-pr_t));
  h_p_mom_diff_phi->Fill(pg_m    ,(pg_p-pr_p));
  h_p_theta_diff_mom->Fill(pg_t  ,(pg_m-pr_m)/pg_m);
  h_p_theta_diff_theta->Fill(pg_t,(pg_t-pr_t));
  h_p_theta_diff_phi->Fill(pg_t  ,(pg_p-pr_p));
  h_p_phi_diff_mom->Fill(pg_p    ,(pg_m-pr_m)/pg_m);
  h_p_phi_diff_theta->Fill(pg_p  ,(pg_t-pr_t));
  h_p_phi_diff_phi->Fill(pg_p    ,(pg_p-pr_p));

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
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Draw("colz");
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
