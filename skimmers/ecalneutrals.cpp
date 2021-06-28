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

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"
#include "e_pid.h"
#include "DC_fiducial.h"

using namespace std;

//from readhipo_helper maxParticles	= 100;
//from readhipo_helper maxEcalhits = 100;
//from readhipo_helper maxNeutrons = 200;
//TODO:  Check for memory leaks(memset??).  Update RUN::config read to BConfig read


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

	if (MC_DATA_OPT!=1) {

	 cerr << "Option for MC_DATA_OPT not implemented. Only option 1 works. Current value is " << MC_DATA_OPT << endl;
	 return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");

	//	Event info:
	int Runno = 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	int eventnumber = 0;

	//	Electron info:
	clashit eHit;

 // Information from REC::Calorimeter for PCAL/ECAL for neutral Particles
  int ecalHits = 0;
	int hit_pindex [maxEcalhits]= {0};
	int hit_pid [maxEcalhits]= {0};
	int hit_sector [maxEcalhits]= {0};
	int hit_layer [maxEcalhits]= {0};
	int hit_detid [maxEcalhits]= {0};
	double hit_energy [maxEcalhits]= {0.};
	double hit_time [maxEcalhits]= {0.};
	double hit_x [maxEcalhits]= {0.};
	double hit_y [maxEcalhits]= {0.};
	double hit_z [maxEcalhits]= {0.};
	double hit_path [maxEcalhits]= {0.};


	//Branches:
	//Event
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("eventnumber",&eventnumber);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);


 //Calorimeter hits
  outTree->Branch("ecalHits"		,&ecalHits	);
  outTree->Branch("hit_pindex"	,&hit_pindex		,"hit_pindex[ecalHits]/I"	);
  outTree->Branch("hit_pid"	,&hit_pid		,"hit_pid[ecalHits]/I"	);
  outTree->Branch("hit_detid"	,&hit_detid		,"hit_detid[ecalHits]/I"	);
  outTree->Branch("hit_sector"	,&hit_sector		,"hit_sector[ecalHits]/I"	);
  outTree->Branch("hit_layer"	,&hit_layer		,"hit_layer[ecalHits]/I"	);
  outTree->Branch("hit_energy"	,&hit_energy		,"hit_energy[ecalHits]/D"	);
	outTree->Branch("hit_time"	,&hit_time		,"hit_time[ecalHits]/D"	);
	outTree->Branch("hit_x"	,&hit_x		,"hit_x[ecalHits]/D"	);
  outTree->Branch("hit_y"	,&hit_y		,"hit_y[ecalHits]/D"	);
	outTree->Branch("hit_z"	,&hit_z		,"hit_z[ecalHits]/D"	);
	outTree->Branch("hit_path"	,&hit_path		,"hit_path[ecalHits]/D"	);


	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");


	// Load the electron PID class:
	e_pid ePID;
	// Load the DC fiducial class for electrons;
	DCFiducial DCfid_electrons;

	// Load input file
	for( int i = 3 ; i < argc ; i++ ){
		if( MC_DATA_OPT == 1){ //Data
			// Using run number of current file, grab the beam energy from RCDB
			int runNum = getRunNumber(argv[i]);
			Runno = runNum;
			auto cnd = connection.GetCondition(runNum, "beam_energy");
			Ebeam = cnd->ToDouble() / 1000.;// [GeV] -- conversion factor
			if (runNum >= 11286 && runNum < 11304)
			{
				//Ebeam *= 1.018; //fudge factor for Low energy run due to miscalibration in RCDB
				Ebeam = 4.244; //fix beam energy for low energy run to currently known number 02/08/21
			}
			current = connection.GetCondition( runNum, "beam_current") ->ToDouble(); // [nA]
		}
		else{
			cout << "Wrong option for MC_DATA. Option is " << MC_DATA_OPT  << ".Exit program" << endl;
			exit(-1);
		}

		//Set cut parameters for electron PID. This only has 10.2 and 10.6 implemented
		ePID.setParamsRGB(Ebeam);

		// Setup hipo reading for this file
		TString inputFile = argv[i];
		hipo::reader reader;
		reader.open(inputFile);
		hipo::dictionary  factory;
		reader.readDictionary(factory);
		BEvent		event_info		(factory.getSchema("REC::Event"		));
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank      DC_Track      (factory.getSchema("REC::Track"         ));
		hipo::bank      DC_Traj      (factory.getSchema("REC::Traj"          ));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::bank  run_config (factory.getSchema("RUN::config"));
		hipo::event 	readevent;

		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		int run_number_from_run_config = 0;
		double torussetting = 0;
		while(reader.next()==true){
				// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			eventnumber = 0;
			// Electron
			eHit.Clear();




			ecalHits = 0;
			memset(	hit_pindex	,0	,sizeof(hit_pindex		)	);
			memset(	hit_pid	,0	,sizeof(hit_pid		)	);
			memset(	hit_detid		,0	,sizeof(hit_detid		)	);
			memset(	hit_sector		,0	,sizeof(hit_sector		)	);
			memset(	hit_layer		,0	,sizeof(hit_layer		)	);
			memset(	hit_energy		,0	,sizeof(hit_energy		)	);
			memset( hit_time	,0	,sizeof(hit_time	)	);
			memset( hit_x	,0	,sizeof(hit_x	)	);
			memset( hit_y	,0	,sizeof(hit_y	)	);
			memset( hit_z	,0	,sizeof(hit_z	)	);
			memset( hit_path	,0	,sizeof(hit_path	)	);




			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if( event_counter > 100 ) continue;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			readevent.getStructure(run_config);
			//Electrons and Other particles
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(DC_Track);
			readevent.getStructure(DC_Traj);

			// Get integrated charge, livetime and start-time from REC::Event
			//Currently, REC::Event has uncalibrated livetime / charge, so these will have to work
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			//Get Event number and run number from RUN::config
			run_number_from_run_config = run_config.getInt( 0 , 0 );
			eventnumber = run_config.getInt( 1 , 0 );
			if (run_number_from_run_config != Runno && event_counter < 100) {
				cout << "Run number from RUN::config and file name not the same!! File name is " << Runno << " and RUN::config is " << run_number_from_run_config << endl;
			}

			//from first event get RUN::config torus Setting
		 // inbending = negative torussetting, outbending = torusseting
			torussetting = run_config.getFloat( 7 , 0 );

			//if (event_counter < 100) {
			//	cout << "event number " << eventnumber << " , runnumebr " << run_number_from_run_config << endl;
			//}

			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , DC_Track, DC_Traj, 0, eHit , starttime , Runno , Ebeam );

			//check electron PID in EC with Andrew's class
			if( !(ePID.isElectron(&eHit)) ) continue;


			//bending field of torus for DC fiducial class ( 1 = inbeding, 0 = outbending	)
			int bending;
			//picking up torussetting from RUN::config, inbending = negative torussetting, outbending = positive torusseting
			if (torussetting > 0 && torussetting <=1.0) { //outbending
				bending = 0;
			}
			else if (torussetting < 0 && torussetting >=-1.0) { //inbending
				bending = 1;
			}
			else {
				cout << "WARNING: Torus setting from RUN::config is " << torussetting << ". This is not defined for bending value for DC fiducials. Please check " << endl;
			}
			if (eHit.getDC_sector() == -999 || eHit.getDC_sector() == -1  ) {
				cout << "Skimmer Error: DC sector is  " << eHit.getDC_sector() << " . Skipping event "<< event_counter << endl;
				eHit.Print();
				continue;
			}

			//checking DC Fiducials
			//Region 1, true = pass DC Region 1
			bool DC_fid_1  = DCfid_electrons.DC_e_fid(eHit.getDC_x1(),eHit.getDC_y1(),eHit.getDC_sector(), 1, bending);
			//checking DC Fiducials
			//Region 2, true = pass DC Region 2
			bool DC_fid_2  = DCfid_electrons.DC_e_fid(eHit.getDC_x2(),eHit.getDC_y2(),eHit.getDC_sector(), 2, bending);
			//checking DC Fiducials
			//Region 3, true = pass DC Region 3
			bool DC_fid_3  = DCfid_electrons.DC_e_fid(eHit.getDC_x3(),eHit.getDC_y3(),eHit.getDC_sector(), 3, bending);

			//check if any of the fiducials is false i.e. electron does not pass all DC fiducials
			if (!DC_fid_1 || !DC_fid_2 || !DC_fid_3) continue;


			//loop over calorimeter bank and store hits if associated particle in REC::Particle bank is photon or neutron
			for( int row = 0 ; row < calorimeter.getRows() ; row++ ){ // loop over calorimeter bank

				int calo_pindex = calorimeter.getPindex(row);
				if (particles.getPid(calo_pindex) == 22 || particles.getPid(calo_pindex) == 2112 ) { //check if calo hit is correlated with photon or neutron in particle PID
					hit_pindex[ecalHits] = calo_pindex;
					hit_pid[ecalHits] = particles.getPid(calo_pindex);
					hit_sector[ecalHits] = calorimeter.getSector(row);
					hit_detid[ecalHits] = calorimeter.getDetector(row);
					hit_layer[ecalHits] = calorimeter.getLayer(row);
					hit_energy[ecalHits] = calorimeter.getEnergy(row);
					hit_time[ecalHits] = calorimeter.getTime(row);
					hit_path[ecalHits] = calorimeter.getPath(row);
					hit_x[ecalHits] = calorimeter.getX(row);
					hit_y[ecalHits] = calorimeter.getY(row);
					hit_z[ecalHits] = calorimeter.getZ(row);

					ecalHits++;

				}
			}



			//Print out of calo hits for debug
			/*for (int ch = 0; ch < ecalHits; ch++) {
				cout << "hit " << ch << " pindex " << hit_pindex[ch] << " , pid " << hit_pid[ch] << " , detid " << hit_detid[ch] << " , layer " << 	hit_layer[ch];
				cout << " , x " << hit_x[ch] << " , y " << hit_y[ch] << endl;
			}*/
			// Fill tree based on d(e,e') for data with all neutron and other particle information
				outTree->Fill();
		//	}


		} // end loop over events
	//	cout << "Total charge collected in file: " << gated_charge << " [microC]\n";
		cout << "Total number of events written to output: " << outTree->GetEntries() << "\n";
	}// end loop over files


	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}
