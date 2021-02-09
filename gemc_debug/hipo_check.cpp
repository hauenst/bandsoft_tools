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
#include "TClonesArray.h"

#include "reader.h"
#include "bank.h"
#include "clas12fiducial.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 3 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
		return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	vector<TH2*> hist_list_2;
	TH2D * h_ept_g_a = new TH2D("e_ept_g_a","Electron;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ept_g_a);
	TH2D * h_ept_g_0 = new TH2D("e_ept_g_0","Electron;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ept_g_0);
	TH2D * h_ept_g_1 = new TH2D("e_ept_g_1","Electron;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ept_g_1);
	TH2D * h_ept_r_a = new TH2D("e_ept_r_a","Electron;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ept_r_a);
	TH2D * h_ept_r_0 = new TH2D("e_ept_r_0","Electron;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ept_r_0);
	TH2D * h_ept_r_1 = new TH2D("e_ept_r_1","Electron;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ept_r_1);
	TH2D * h_ppt_g_a = new TH2D("e_ppt_g_a","Proton;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ppt_g_a);
	TH2D * h_ppt_g_0 = new TH2D("e_ppt_g_0","Proton;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ppt_g_0);
	TH2D * h_ppt_g_1 = new TH2D("e_ppt_g_1","Proton;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ppt_g_1);
	TH2D * h_ppt_r_a = new TH2D("e_ppt_r_a","Proton;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ppt_r_a);
	TH2D * h_ppt_r_0 = new TH2D("e_ppt_r_0","Proton;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ppt_r_0);
	TH2D * h_ppt_r_1 = new TH2D("e_ppt_r_1","Proton;Phi;Theta;Events",360,-180,180,60,10,40);
	hist_list_2.push_back(h_ppt_r_1);


	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Get CLAS12 fiducials
	clas12fiducial* fFiducial = new clas12fiducial();

	// Load input file
	for( int i = 3 ; i < argc ; i++ ){

	  int runNum = 1;
	  Runno = runNum;

	  // Setup hipo reading for this file
	  TString inputFile = argv[i];
	  hipo::reader reader;
	  reader.open(inputFile);
	  hipo::dictionary  factory;      
	  hipo::schema	  schema;
	  reader.readDictionary(factory); 
	  BEvent		event_info		(factory.getSchema("REC::Event"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::bank      DC_Track                (factory.getSchema("REC::Track"         ));
		hipo::bank      DC_Traj                 (factory.getSchema("REC::Traj"          ));
		hipo::event 	readevent;
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank	mc_event_info		(factory.getSchema("MC::Event"		));
		hipo::bank	mc_particle		(factory.getSchema("MC::Particle"	));
		
		// Loop over all events in file
		int event_counter = 0;
		double gated_charge = 0;
		double livetime	= 0;

		while(reader.next()==true){
		  gated_charge	= 0;
		  livetime	= 0;
		  
		  // Count events
			if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			readevent.getStructure(mc_event_info);
			readevent.getStructure(mc_particle);
			// electron struct
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(DC_Track);
			readevent.getStructure(DC_Traj);
	
			if( event_info.getRows() == 0 ) continue;
			
			// electron variables
			TVector3	momentum = particles.getV3P(0);
			TVector3	vertex	 = particles.getV3v(0);
			double ETot = calorimeter.getTotE(0);
			double EoP = ETot / momentum.Mag();
			double U = scintillator.getLU(0);
			double V = scintillator.getLV(0);
			double W = scintillator.getLW(0);
			// For simulated events, get the weight for the event
			getMcInfo( mc_particle , mc_event_info , mcPart , starttime, weight, Ebeam , genMult );

			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , DC_Track, DC_Traj, eHit , starttime , Runno , Ebeam );



		} // end loop over events
	}// end loop over files
	
	outFile->cd();
	for(int i=0; i<hist_list_2.size(); i++){
	  hist_list_2[i]->Write();
	}
	outFile->Close();

	return 0;
}



