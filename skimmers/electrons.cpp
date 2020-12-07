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
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [apply fiducial: 0 = n, 1 = y] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[apply fiducial = <0, 1> \n";
		cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
		return -1;
	}

	bool applyFiducial = false;
	if (atoi(argv[2]) == 1) {
		applyFiducial = true;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("electrons","CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	//	Electron info:
	clashit eHit;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Get CLAS12 fiducials
	clas12fiducial* fFiducial = new clas12fiducial();

	// Load input file
	for( int i = 3 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		int runNum = getRunNumber(argv[i]);
	        Runno = runNum;
		auto cnd = connection.GetCondition(runNum, "beam_energy");
		Ebeam = cnd->ToDouble() / 1000.; // [GeV]
		current = connection.GetCondition( runNum, "beam_current") ->ToDouble(); // [nA]

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
		
		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;

                //int count = 0;

                //int count1 =0;
                //int count_mul =0;

		while(reader.next()==true){
			// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			eHit.Clear();

			// Count events
			if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
			event_counter++;
			//if( event_counter > 100000 ) break;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			// electron struct
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(DC_Track);
			readevent.getStructure(DC_Traj);
	
			// Get integrated charge, livetime and start-time from REC::Event
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , DC_Track, DC_Traj, eHit , starttime , Runno , Ebeam );


			//I want to get the track information here to test
                        //int count_track = 0;

			//int test_Ntrack = DC_Track.getRows();
		        //for(int i =0; i <test_Ntrack; i++){
	                //    int pindex = DC_Track.getInt(1,i);
	                //    int detector = DC_Track.getInt(2,i);

	                //    if (pindex ==0 && detector ==6)
			//      {	count_track ++;}

			//}
                        //
                        //if(count_track ==1){count1++;}
                        //if(count_track >1){count_mul++;}

			//Done checking tracking information

			bool eAccept = true;
			if(applyFiducial) {
				int sect = fFiducial->GetElectronAcceptance(eHit.getTheta(), eHit.getPhi(), eHit.getMomentum());
				if(sect < 0) {
					eAccept = false;
				}	
			}

			if(eAccept) {
			    outTree->Fill();
			}

			//count ++;

			 //if (count == 1000000) break;

		} // end loop over events
	}// end loop over files
	
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}



