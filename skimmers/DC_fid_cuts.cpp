// I just want to write a simple code
// access to inc_data_skimed file
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "TApplication.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TH2.h"
#include "TRint.h"
#include "TStyle.h"
#include "TTree.h"
#include "TF1.h"

#include "bandhit.h"
#include "clas12fiducial.h"
#include "clashit.h"
#include "constants.h"
#include "taghit.h"
#include "DC_fiducial.h"


using namespace std;

int main(int argc, char **argv) {

  // These lines of code to make sure we can see the plots

#ifdef WITHRINT
  // TRint *myapp = new TRint("RootSession", &argc, argv, NULL, 0);
  TRint *myapp = new TRint("RootSession", 0, 0, NULL, 0);
#else
  TApplication *myapp = new TApplication("myapp", 0, 0);
#endif

  cout << argv[0] << "  " << argv[1] << " " << argv[3] << endl;
  // Star the program

  if (argc < 5) {
    cerr
        << "Wrong number of arguments. Instead use\n"
        << "\t./code [outfile.root] [apply fiducial (1, 0)] [bending (1- inbending/0 -outbending)] [inputDatafiles]\n";
       
    return -1;
  }

  // output rootfile
  TFile *outFile = new TFile(argv[1], "RECREATE");

  cout << "test oufile: " << outFile->GetName() << endl;
  // Define histogram to look at e-fiducial cut

  TH2D *h2_DC1_xy = new TH2D("h2_DC1_xy", "DC: ele Region1", 300, -300, 300,300, -300, 300);
  TH2D *h2_DC2_xy = new TH2D("h2_DC2_xy", "DC: ele Region2", 300, -300, 300,300, -300, 300);
  TH2D *h2_DC3_xy = new TH2D("h2_DC3_xy", "DC: ele Region3", 300, -300, 300,300, -300, 300);

  TH2D *h2_DC1_xy_acc = new TH2D("h2_DC1_xy_acc", "DC: ele Region1 acc", 300, -300, 300,300, -300, 300);
  TH2D *h2_DC2_xy_acc = new TH2D("h2_DC2_xy_acc", "DC: ele Region2 acc", 300, -300, 300,300, -300, 300);
  TH2D *h2_DC3_xy_acc = new TH2D("h2_DC3_xy_acc", "DC: ele Region3 acc", 300, -300, 300,300, -300, 300);
  
  //Histogram for each each section in each layer
  int Nsec[6] = {1, 2, 3, 4, 5, 6};
  int Nreg[3] = {6, 18, 36};

  TH2D *h2_DC_xy1_S[6];
  TH2D *h2_DC_xy2_S[6];
  TH2D *h2_DC_xy3_S[6];
 
    for (int s =1; s<7; s++){
      h2_DC_xy1_S[s-1] = new TH2D(Form("h2_DC_xy1_S_%d",s), Form("Region 1, sector %d",s),150, 0, 150, 300, -150, 150);
      h2_DC_xy2_S[s-1] = new TH2D(Form("h2_DC_xy2_S_%d",s), Form("Region 2, sector %d",s),250, 0, 250, 300, -150, 150);
      h2_DC_xy3_S[s-1] = new TH2D(Form("h2_DC_xy3_S_%d",s), Form("Region 3, sector %d",s),300, 0, 300, 300, -150, 150);  
    }
    

  // Def the fiducial cut
  int doFiducial = atoi(argv[2]);
  clas12fiducial *fid = new clas12fiducial();

  DCFiducial DCfid;
  

  // Def bending setting
  int bending = atoi(argv[3]);


  // access to the data root file

  for (int i = 4; i < argc; i++) {
    TFile *inFile = new TFile(argv[i]);
    TTree *inTree;

    // I only care electron tree
    inTree = (TTree *)inFile->Get("electrons");

    // intitial the input branches
    int Runno = 0;
    double Ebeam = 0;
    double gated_charge = 0;
    double livetime = 0;
    double starttime = 0;
    double current = 0;
    clashit *eHit = new clashit;

    inTree->SetBranchAddress("Runno", &Runno);
    inTree->SetBranchAddress("Ebeam", &Ebeam);
    inTree->SetBranchAddress("gated_charge", &gated_charge);
    inTree->SetBranchAddress("livetime", &livetime);
    inTree->SetBranchAddress("starttime", &starttime);
    inTree->SetBranchAddress("current", &current);

    inTree->SetBranchAddress("eHit", &eHit);

    // checking the file we are working on
    cout << "Working on file: " << argv[i] << endl;
    cout << "Number of event in the file" << inTree->GetEntries() << endl;

    

    //Define 3vector for each event for rotation
    TVector3 DC_xyz1, DC_xyz2, DC_xyz3;
    
    //Define accepted varibles for each region: 
    bool DC1_fid = false;
    bool DC2_fid = false;
    bool DC3_fid = false;

    int count1 =0, count1_acc =0;
    int count2 =0, count2_acc =0;
    int count3 =0, count3_acc =0;
   
    // Loop over number of event
    for (int ev = 0; ev < inTree->GetEntries(); ev++) {
      if (ev % 100000 == 0)
        cout << "Test:" << ev << endl;
      // Clear all branches before getting the entry from tree
      Runno = 0;
      Ebeam = 0;
      gated_charge = 0;
      livetime = 0;
      starttime = 0;
      current = 0;
      eHit->Clear();

      inTree->GetEntry(ev);
      // Check electron information
      if (eHit->getPID() != 11)
        continue;
      if (eHit->getCharge() != -1)
        continue;
      /*
      if (eHit->getEoP() < 0.17)
        continue;
      if (eHit->getEoP() > 0.3)
        continue;
      if (eHit->getEpcal() < 0.07)
        continue;
      if (eHit->getV() < 9)
        continue;
      if (eHit->getW() < 9)
        continue;
      if (eHit->getVtz() < -8)
        continue;
      if (eHit->getVtz() > 3)
        continue;
      */
      if (eHit->getMomentum() < 2.)
        continue;

    
      if (doFiducial) {
        int eSect = fid->GetElectronAcceptance(
            eHit->getTheta() * TMath::RadToDeg(),
            eHit->getPhi() * TMath::RadToDeg(), eHit->getMomentum());
        if (eSect < 0)
          continue;
      }
 

      //Plot global xy for all 6 sectors for 3 region

      h2_DC1_xy->Fill(eHit->getDC_x1(), eHit->getDC_y1());
      h2_DC2_xy->Fill(eHit->getDC_x2(), eHit->getDC_y2());
      h2_DC3_xy->Fill(eHit->getDC_x3(), eHit->getDC_y3());
      
      //Now need to check to plots in Sector coordinator system
      //And get plots for each sector for each layer 

      int sector_ID = eHit->getDC_sector();
      count1++;
      count2++;
      count3++;

      //Doing rotation to sector coordinator system
      DC_xyz1 = DCfid.rotate(eHit->getDC_x1(), eHit->getDC_y1(), sector_ID);
      DC_xyz2 = DCfid.rotate(eHit->getDC_x2(), eHit->getDC_y2(), sector_ID);
      DC_xyz3 = DCfid.rotate(eHit->getDC_x3(), eHit->getDC_y3(), sector_ID);
      //Filling the histogram for each sector

      h2_DC_xy1_S[sector_ID-1]->Fill(DC_xyz1.X(), DC_xyz1.Y());
      h2_DC_xy2_S[sector_ID-1]->Fill(DC_xyz2.X(), DC_xyz2.Y());
      h2_DC_xy3_S[sector_ID-1]->Fill(DC_xyz3.X(), DC_xyz3.Y());

      
      //======Checking the fiducial cut for each layer =========//
    
      //checking DC fiducial for layer 1
      DC1_fid = DCfid.DC_e_fid(eHit->getDC_x1(),eHit->getDC_y1(),sector_ID, 1, bending);   

     if(DC1_fid){
          h2_DC1_xy_acc->Fill(eHit->getDC_x1(), eHit->getDC_y1());
	  count1_acc++;
	 }

       //checking DC fiducial for layer 2
      DC2_fid = DCfid.DC_e_fid(eHit->getDC_x2(),eHit->getDC_y2(),sector_ID, 2, bending);

      if(DC2_fid){
          h2_DC2_xy_acc->Fill(eHit->getDC_x2(), eHit->getDC_y2());
	  count2_acc++;
      }

       //checking DC fiducial for layer 3
      DC3_fid = DCfid.DC_e_fid(eHit->getDC_x3(),eHit->getDC_y3(),sector_ID, 3, bending);

      if(DC3_fid){
          h2_DC3_xy_acc->Fill(eHit->getDC_x3(), eHit->getDC_y3());
	  count3_acc++;
      }   

    } // end loop event

    cout << "count1, count2, count3: "<< count1<< "  "<< count2<< "  "<< count3<< endl;
    cout << "count1_acc, count2_acc, count3_acc: "<< count1_acc<< "  "<< count2_acc<< "  "<< count3_acc<< endl;
    cout << "acc: "<< (double)count1_acc/count1<< "   "<< (double)count2_acc/count2<< "  "<< (double)count3_acc/count3<< endl;

    // inFile->Close();
  } // end loop over file

  //=========Doing plots============//


  //Plot for region1 
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  h2_DC1_xy->Draw();
  h2_DC1_xy->SetMarkerColor(2);
  h2_DC1_xy->SetMarkerSize(2);

  h2_DC1_xy->GetXaxis()->SetTitle("X [cm]");
  h2_DC1_xy->GetYaxis()->SetTitle("Y [cm]");
  h2_DC1_xy_acc->Draw("colz same");
  
  TCanvas *c11 = new TCanvas("c11", "", 1200, 800);
   c11->Divide(3,2);
   for(int i =1; i<7; i++){
    c11->cd(i);
    h2_DC_xy1_S[i-1]->Draw("colz");
    h2_DC_xy1_S[i-1]->GetXaxis()->SetTitle("X [cm]");
    h2_DC_xy1_S[i-1]->GetYaxis()->SetTitle("Y [cm]");

    DCfid.get_fmin_in(i,1)->Draw("same");
    DCfid.get_fmin_in(i,1)->SetLineColor(2);

    DCfid.get_fmax_in(i,1)->Draw("same");
    DCfid.get_fmax_in(i,1)->SetLineColor(2);
  }
  
  //Plot for region2
   TCanvas *c2 = new TCanvas("c2", "", 800, 600);
   h2_DC2_xy->Draw("");
   h2_DC2_xy->SetMarkerColor(2);
   h2_DC2_xy->SetMarkerSize(2);

   h2_DC2_xy->GetXaxis()->SetTitle("X [cm]");
   h2_DC2_xy->GetYaxis()->SetTitle("Y [cm]");
   h2_DC2_xy_acc->Draw("colz same");

   TCanvas *c21 = new TCanvas("c21", "", 1200, 800);
   c21->Divide(3,2);
   for(int i =1; i<7; i++){
    c21->cd(i);
    h2_DC_xy2_S[i-1]->Draw("colz");
    h2_DC_xy2_S[i-1]->GetXaxis()->SetTitle("X [cm]");
    h2_DC_xy2_S[i-1]->GetYaxis()->SetTitle("Y [cm]");
    
    DCfid.get_fmin_in(i,2)->Draw("same");
    DCfid.get_fmin_in(i,2)->SetLineColor(2);

    DCfid.get_fmax_in(i,2)->Draw("same");
    DCfid.get_fmax_in(i,2)->SetLineColor(2);
 
   }

   
    //Plot for region 3
   TCanvas *c3 = new TCanvas("c3", "", 800, 600);
   h2_DC3_xy->Draw();
   h2_DC3_xy->SetMarkerColor(2);
   h2_DC3_xy->SetMarkerSize(2);

   h2_DC3_xy->GetXaxis()->SetTitle("X [cm]");
   h2_DC3_xy->GetYaxis()->SetTitle("Y [cm]");
   h2_DC3_xy_acc->Draw("colz same");

   TCanvas *c31 = new TCanvas("c31", "", 1200, 800);
   c31->Divide(3,2);
   for(int i =1; i<7; i++){
    c31->cd(i);
    h2_DC_xy3_S[i-1]->Draw("colz");
    
    DCfid.get_fmin_in(i,3)->Draw("same");
    DCfid.get_fmin_in(i,3)->SetLineColor(2);

    DCfid.get_fmax_in(i,3)->Draw("same");
    DCfid.get_fmax_in(i,3)->SetLineColor(2);
    
   }
  

   // outFile->cd();
   // outFile->Close();

   c1->Print("DC_fid_plot.pdf[");
   c1->Print("DC_fid_plot.pdf");
   c11->Print("DC_fid_plot.pdf");
   c2->Print("DC_fid_plot.pdf");
   c21->Print("DC_fid_plot.pdf");
   c3->Print("DC_fid_plot.pdf");
   c31->Print("DC_fid_plot.pdf");
   c31->Print("DC_fid_plot.pdf]");

  cout << "All Done: " << endl;

  myapp->Run();
  return 0;
}
