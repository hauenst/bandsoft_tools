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

#include "reader.h"
#include "bank.h"
#include "clas12fiducial.h"
#include "e_pid.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;
double diag(double x, double C){
  return C - x;
}

int main(int argc, char** argv) {

  if( argc != 4 ){
    cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [inputFile] [Ebeam]\n\n";
    return -1;
  }
  

  //Define Variables
  int Runno;
  double Ebeam, gated_charge, livetime, starttime, current;
  clashit * eHit = new clashit;
 //Creat input tree
  TFile * inFile = new TFile(argv[2]);
  TTree * inTree = (TTree*)inFile->Get("electrons");
  // 	Event branches:
  inTree->SetBranchAddress("Runno"               ,&Runno                 );
  inTree->SetBranchAddress("Ebeam"               ,&Ebeam                 );
  inTree->SetBranchAddress("gated_charge"        ,&gated_charge          );
  inTree->SetBranchAddress("livetime"            ,&livetime              );
  inTree->SetBranchAddress("starttime"           ,&starttime             );
  inTree->SetBranchAddress("current"             ,&current               );
  //	Electron branches:
  inTree->SetBranchAddress("eHit"		 ,&eHit			 );

  //Make some histograms
  vector<TH2*> hist_list;
  TH2D * h2_SF_Epcal_nc[6];
  TH2D * h2_SF_Epcal_lc[6];
  TH2D * h2_SF_Epcal_fc[6];
  TH2D * h2_SF_Mom_nc[6];
  TH2D * h2_SF_Mom_lc[6];
  TH2D * h2_SF_Mom_fc[6];

  char temp[100];
  char temp2[100];

  for(int i = 0; i<6; i++){

    sprintf(temp,"SF_Epcal_nc_sec%d",i);
    sprintf(temp2,"SF_Epcal_nc_sec%d;Epcal;SF;Events",i);
    h2_SF_Epcal_nc[i] = new TH2D(temp,temp2,250,0,1.7,150,0.05,0.40);
    hist_list.push_back(h2_SF_Epcal_nc[i]);

    sprintf(temp,"SF_Epcal_lc_sec%d",i);
    sprintf(temp2,"SF_Epcal_lc_sec%d;Epcal;SF;Events",i);
    h2_SF_Epcal_lc[i] = new TH2D(temp,temp2,250,0,1.7,150,0.05,0.40);
    hist_list.push_back(h2_SF_Epcal_lc[i]);

    sprintf(temp,"SF_Epcal_fc_sec%d",i);
    sprintf(temp2,"SF_Epcal_fc_sec%d;Epcal;SF;Events",i);
    h2_SF_Epcal_fc[i] = new TH2D(temp,temp2,250,0,1.7,150,0.05,0.40);
    hist_list.push_back(h2_SF_Epcal_fc[i]);

    sprintf(temp,"SF_Mom_nc_sec%d",i);
    sprintf(temp2,"SF_Mom_nc_sec%d;Mom;SF;Events",i);
    h2_SF_Mom_nc[i] = new TH2D(temp,temp2,250,0,10,150,0.05,0.40);
    hist_list.push_back(h2_SF_Mom_nc[i]);

    sprintf(temp,"SF_Mom_lc_sec%d",i);
    sprintf(temp2,"SF_Mom_lc_sec%d;Mom;SF;Events",i);
    h2_SF_Mom_lc[i] = new TH2D(temp,temp2,250,0,10,150,0.05,0.40);
    hist_list.push_back(h2_SF_Mom_lc[i]);

    sprintf(temp,"SF_Mom_fc_sec%d",i);
    sprintf(temp2,"SF_Mom_fc_sec%d;Mom;SF;Events",i);
    h2_SF_Mom_fc[i] = new TH2D(temp,temp2,250,0,10,150,0.05,0.40);
    hist_list.push_back(h2_SF_Mom_fc[i]);

  }
 
  int fin = inTree->GetEntries();
  
  //Initiallize PID with beam energy
  inTree->GetEntry(0);
  e_pid ePID_b;
  ePID_b.setParamsRGB(atof(argv[3]));

  for(int i = 0; i < fin; i++){
    eHit->Clear();
    inTree->GetEntry(i);
    
    //cout<<"Ebeam = " << Ebeam << "\n";
    //cout<<"Runno = " << Runno << "\n";
    //cout<<"EoP = " << eHit->getEoP() << "\n";

    //Display completed  
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }

    //Do ECal Fid and PID cuts
    if(eHit->getV() < 14){ continue; }
    if(eHit->getW() < 14){ continue; }

    int sector = eHit->getSector() - 1;

    h2_SF_Epcal_nc[sector]->Fill(eHit->getEpcal(),eHit->getEoP());
    h2_SF_Mom_nc[sector]->Fill(eHit->getMomentum(),eHit->getEoP());

    if(ePID_b.isElectronLoose(eHit)){
      h2_SF_Epcal_lc[sector]->Fill(eHit->getEpcal(),eHit->getEoP());
      h2_SF_Mom_lc[sector]->Fill(eHit->getMomentum(),eHit->getEoP());
    }

    if(ePID_b.isElectron(eHit)){
      h2_SF_Epcal_fc[sector]->Fill(eHit->getEpcal(),eHit->getEoP());
      h2_SF_Mom_fc[sector]->Fill(eHit->getMomentum(),eHit->getEoP());
    }

  }
  
  char fileName[100];
  sprintf(fileName,"%s[",argv[1]);
  TCanvas * myCanvas = new TCanvas("My_Canvas","My_Canvas",1200,1000);
  myCanvas->SaveAs(fileName);

  sprintf(fileName,"%s",argv[1]);
  for(int j = 0; j < 6; j++){
    myCanvas->cd();

    h2_SF_Epcal_nc[j]->Draw("colz");
    ePID_b.drawEpcal(j+1,myCanvas);
    myCanvas->SaveAs(fileName);

    h2_SF_Epcal_lc[j]->Draw("colz");
    ePID_b.drawEpcal(j+1,myCanvas);
    myCanvas->SaveAs(fileName);

    h2_SF_Epcal_fc[j]->Draw("colz");
    ePID_b.drawEpcal(j+1,myCanvas);
    myCanvas->SaveAs(fileName);

    h2_SF_Mom_nc[j]->Draw("colz");
    ePID_b.drawMom(j+1,myCanvas);
    myCanvas->SaveAs(fileName);

    h2_SF_Mom_lc[j]->Draw("colz");
    ePID_b.drawMom(j+1,myCanvas);
    myCanvas->SaveAs(fileName);

    h2_SF_Mom_fc[j]->Draw("colz");
    ePID_b.drawMom(j+1,myCanvas);
    myCanvas->SaveAs(fileName);
  }

  sprintf(fileName,"%s]",argv[1]);
  TCanvas * endCanvas = new TCanvas("End_Canvas","End_Canvas",1200,1000);
  myCanvas->SaveAs(fileName);

  inFile->Close();

  cout<<"Finished making files in: "<<argv[1]<<"\n";
  return 0;
}
