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

#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"
#include "e_pid.h"
#include "DC_fiducial.h"

#include "bandreco.h"


using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 5 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [Period] [load shifts] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[<MC,DATA> = <0, 1> \n";
		cerr << "\t\t[Periods are:\n";
                cerr << "\t\t\t0 == 10.6,\n\t\t\t1 == 10.2,\n\t\t\t2 ==  10.4,\n\t\t\t3 == LER,\n\t\t\t4 == RGM 2.1,\n\t\t\t5 == RGM 4.0,\n\t\t\t6 == RGM 5.98 pt 1,\n\t\t\t7 == RGM 5.98 pt 2,\n\t\t\t8 == RGD,\n\t\t\t9 == RGE,\n\t\t\t10 == RGL 2.2,\n\t\t\t11 == RGL 6.4,\n\t\t\t12 == RGL 10.6]\n";
		cerr << "\t\t[load shifts N,Y] = 0, 1 \n";
		cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
		return -1;
	}

	int MC_DATA_OPT = atoi(argv[2]);
	int PERIOD = atoi(argv[3]);
	int loadshifts_opt = atoi(argv[4]);

	// Initialize our BAND reconstruction engine:
	BANDReco * BAND = new BANDReco();

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("neutrons","BAND Neutrons and CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	int eventnumber = 0;
	bool goodneutron = false;
	int nleadindex = -1;
	double weight		= 0;
	//	MC info:
	int genMult		= 0;
	TClonesArray * mcParts = new TClonesArray("genpart");
	TClonesArray &saveMC = *mcParts;
	// 	Neutron info:
	int nMult		= 0;
	int passed		= 0;
	TClonesArray * nHits = new TClonesArray("bandhit");
	TClonesArray &saveHit = *nHits;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("eventnumber",&eventnumber);
	outTree->Branch("weight"	,&weight		);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("passed"	,&passed		);
	outTree->Branch("nHits"		,&nHits			);
	//Branches to store if good Neutron event and leadindex
	outTree->Branch("goodneutron"		,&goodneutron	);
	outTree->Branch("nleadindex"		,&nleadindex			);
	//	MC branches:
	outTree->Branch("genMult"	,&genMult		);
	outTree->Branch("mcParts"	,&mcParts		);

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

		// Load the electron PID class:
		e_pid ePID;
		// Load the DC fiducial class for electrons;
		DCFiducial DCfid_electrons;

	// Load input file
	for( int i = 5 ; i < argc ; i++ ){
		if( MC_DATA_OPT == 0){
			int runNum = 11;
			Runno = runNum;
			// Set the initial Ebeam value so that it can be used for the PID class
			if( PERIOD == 0 ) Ebeam = 10.599; // from RCDB: 10598.6
			if( PERIOD == 1 ) Ebeam = 10.200; // from RCDB: 10199.8
			if( PERIOD == 2 ) Ebeam = 10.389; // from RCDB: 10389.4
			if( PERIOD == 3 ) Ebeam = 4.247;  // current QE-MC value, RCDB value: 4171.79. Is about ~1.018 wrong due to issues with magnet settings
			if( PERIOD == 4 ) Ebeam = 2.07052;
			if( PERIOD == 5 ) Ebeam = 4.02962;
			if( PERIOD == 6 or PERIOD == 7 ) Ebeam = 5.98636;
			if( PERIOD == 8 or PERIOD == 9) Ebeam = 10.547;
			if( PERIOD == 10 ) Ebeam = 2.239;
			if( PERIOD == 11 ) Ebeam = 6.473;
			if( PERIOD == 12 ) Ebeam = 10.676;
		}
		else if( MC_DATA_OPT == 1){
			int runNum = getRunNumber(argv[i]);
			Runno = runNum;
			auto cnd = connection.GetCondition(runNum, "beam_energy");
			Ebeam = cnd->ToDouble() / 1000.; // [GeV]
			current = connection.GetCondition( runNum, "beam_current") ->ToDouble(); // [nA]
			if (runNum >= 11286 && runNum < 11304){
			// Manual change of Ebeam for LER since RCDB is wrong by ~1.018 due to magnet setting issue
			Ebeam = 4.244; //fix beam energy for low energy run to currently known number 02/08/21
					// NOTE: this does NOT match the MC beam energy by 3MeV because it doesn't matter and
					// we don't know the exact value.
		}
		}
		else{
			exit(-1);
		}
		//Set cut parameters for electron PID. This only has 10.2 and 10.6 implemented
		ePID.setParamsRGB(Ebeam);

		// Setup hipo reading for this file
		TString inputFile = argv[i];
		hipo::reader reader;
		reader.open(inputFile);
		hipo::dictionary  factory;
		hipo::schema	  schema;
		reader.readDictionary(factory);
		BEvent		event_info		(factory.getSchema("REC::Event"		));
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::bank  run_config (factory.getSchema("RUN::config"));
		hipo::event 	readevent;
		hipo::bank	band_rawhits		(factory.getSchema("BAND::rawhits"	));
		hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
		hipo::bank	mc_event_info		(factory.getSchema("MC::Event"		));
		hipo::bank	mc_particle		(factory.getSchema("MC::Particle"	));

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
			// Neutron
			nMult		= 0;
			passed = 0;
			nleadindex = -1;
			goodneutron = false;
			bandhit nHit[maxNeutrons];
			nHits->Clear();
			// MC
			genMult = 0;
			weight = 1;
			genpart mcPart[maxGens];
			mcParts->Clear();

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			readevent.getStructure(run_config);
			// band struct
			readevent.getStructure(band_hits);
			readevent.getStructure(band_rawhits);
			readevent.getStructure(band_adc);
			readevent.getStructure(band_tdc);
			// monte carlo struct
			readevent.getStructure(mc_event_info);
			readevent.getStructure(mc_particle);

			if( event_counter == 1 ){
						int period = -1;
						BAND->setRunno(Runno);
						//Load of shifts depending on run number
						if (Runno > 6100 && Runno < 6400) { //Spring 19 data - 10.6 data
						period = 0;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 6400 && Runno < 6800) { //Spring 19 data - 10.2 data
						period = 1;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno > 11320 && Runno < 11580) { //Spring 20 data - 10.4 data
						period = 2;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 11286 && Runno < 11304) { //LER runs
						period = 3;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 15540 && Runno < 15727) { //RGM 2GeV
						period = 4;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 15731 && Runno < 15784) { //RGM 4GeV
						period = 5;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 15014 && Runno < 15491) { //RGM 6GeV pt1
						period = 6;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 15787 && Runno < 15884) { //RGM 6GeV pt2
						period = 7;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 18305 && Runno < 19131) { //RGD
						period = 8;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 20000 && Runno < 20525) { //RGE
						period = 9;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if ( (Runno >= 21278 && Runno < 21435) or (Runno >= 22709 && Runno <= 22791) or (Runno >= 22981 && Runno <= 23015) ) { //RGL 2.2GeV
						period = 10;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if ( (Runno >= 22793 && Runno < 22980) or (Runno >= 23016 && Runno <= 23065) ) { //RGL 6.4GeV
						period = 11;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}
					else if (Runno >= 21437 && Runno < 22708) { //RGL 10.6GeV
						period = 12;
						if( period != PERIOD ){ cerr << "issue setting period\n...exiting\n"; exit(-1); }
						BAND->setPeriod(period);
					}

					else if( Runno == 11 ){
						// already set the beam energy for MC runs and the period is the user input period
						period = PERIOD;
						BAND->setMC();
						BAND->setPeriod(period); // what is the simulated period (used for status table)
					}
					else {
						cout << "No bar by bar offsets loaded " << endl;
						cout << "Check shift option when starting program. Exit " << endl;
						exit(-1);
					}
					if( period == -1 ){ cerr << "invalid period\n"; exit(-1); }
						BAND->readTW();			// TW calibration values for each PMT
						BAND->readLROffset();		// (L-R) offsets for each bar
						BAND->readPaddleOffset();	// bar offsets relative to bar 2X7 in each layer X
						BAND->readLayerOffset();	// layer offsets relative to layer 5
						BAND->readGeometry();		// geometry table for each bar
						BAND->readEnergyCalib();	// energy calibration for Adc->MeVee
						BAND->readStatus();		// status table for 0,1 = bad,good bar
						if( loadshifts_opt ) BAND->readGlobalOffset();	// final global alignment relative to electron trigger
			}


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

			// For simulated events, get the weight for the event
			if( MC_DATA_OPT == 0){
				getMcInfo( mc_particle , mc_event_info , mcPart , starttime, weight, Ebeam , genMult );
			}

			/// Grab the neutron information:
				// Form the PMTs and Bars for BAND:
			BAND->createPMTs( &band_adc, &band_tdc, &run_config );
			BAND->createBars();
			BAND->storeHits( nMult , nHit , starttime , BAND->getRGBVertexOffset() ); // use average z-offset of the target for pathlength-z


			// Store the neutrons in TClonesArray
			for( int n = 0 ; n < nMult ; n++ ){
				new(saveHit[n]) bandhit;
				saveHit[n] = &nHit[n];
			}
			// Store the mc particles in TClonesArray
			for( int n = 0 ; n < maxGens ; n++ ){
				new(saveMC[n]) genpart;
				saveMC[n] = &mcPart[n];
			}

			if (nMult == 1) {
				goodneutron =  true;
				nleadindex = 0;
			}
			//If nMult > 1: Take nHit and check if good event and give back leading hit index and boolean
			if (nMult > 1) {
				//pass Nhit array, multiplicity and reference to leadindex which will be modified by function

				goodneutron = goodNeutronEvent(nHit, nMult, nleadindex, MC_DATA_OPT,passed);
			}

			// Fill tree based on d(e,n)X for data
					if(  (nMult == 1 || (nMult > 1 && goodneutron) ) && MC_DATA_OPT == 1  ){
						outTree->Fill();
					}
					else if( MC_DATA_OPT == 0 ){
						outTree->Fill();
					}// else fill tree on d(e,n)

		} // end loop over events
	}// end loop over files

	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}
