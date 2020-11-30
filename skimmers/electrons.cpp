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
	if( argc < 4 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[<MC,DATA> = <0, 1> \n";
		cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
		return -1;
	}

	int MC_DATA_OPT = atoi(argv[2]);

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
	double weight		= 0;
	//	MC info:
	int genMult		= 0;
	TClonesArray * mcParts = new TClonesArray("genpart");
	TClonesArray &saveMC = *mcParts;
	//	Electron info:
	clashit eHit;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("weight"	,&weight		);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);
	//	MC branches:
	outTree->Branch("genMult"	,&genMult		);
	outTree->Branch("mcParts"	,&mcParts		);
	
	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Get CLAS12 fiducials
	clas12fiducial* fFiducial = new clas12fiducial();

	// Load input file
	for( int i = 3 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		if( MC_DATA_OPT == 0){
			int runNum = 11;
			Runno = runNum;
			Ebeam = 10.2;
		}
		else if( MC_DATA_OPT == 1){
			int runNum = getRunNumber(argv[i]);
			auto cnd = connection.GetCondition(runNum, "beam_energy");
			Ebeam = cnd->ToDouble() / 1000.; // [GeV]
			current = connection.GetCondition( runNum, "beam_current") ->ToDouble(); // [nA]
		}
		else{
			exit(-1);
		}


		// Setup hipo reading for this file
		TString inputFile = argv[i];
		hipo::reader reader;
		reader.open(inputFile);
		hipo::dictionary  factory;      
		hipo::schema	  schema;
		reader.readDictionary(factory); 
		BEvent		event_info		(factory.getSchema("REC::Event"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank	mc_event_info		(factory.getSchema("MC::Event"		));
		hipo::bank	mc_particle		(factory.getSchema("MC::Particle"	));
		
		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		while(reader.next()==true){
			// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			weight		= 1;
			eHit.Clear();
			// MC
			genMult = 0;
			genpart mcPart[maxGens];
			mcParts->Clear();

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			event_counter++;
			if( event_counter > 1000000 ) break;

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
	
			// Get integrated charge, livetime and start-time from REC::Event
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );
			
			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , eHit , starttime , Runno , Ebeam );

			// For simulated events, get the weight for the event		
			if( MC_DATA_OPT == 0){
				getMcInfo( mc_particle , mc_event_info , mcPart , starttime, weight, Ebeam , genMult );
			}

			bool eAccept = true;
			//if(applyFiducial) {
			//	int sect = fFiducial->GetElectronAcceptance(eHit.getTheta(), eHit.getPhi(), eHit.getMomentum());
			//	if(sect < 0) {
			//		eAccept = false;
			//	}	
			//}

			// Store the mc particles in TClonesArray
			for( int n = 0 ; n < maxGens ; n++ ){
				new(saveMC[n]) genpart;
				saveMC[n] = &mcPart[n];
			}

			if(eAccept) {
				outTree->Fill();
			}

		} // end loop over events
	}// end loop over files
	
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}



