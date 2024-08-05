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

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"
#include "e_pid.h"
#include "DC_fiducial.h"

#include "bandreco.h"

using namespace std;
using namespace QA;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 6 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [Peroid] [load shifts] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[MC,DATA, MC generated info for each event] = 0, 1, 2 \n";
		cerr << "\t\t[Period 10.6, 10.2, 10.4, LER, RGM 2.1] = 0,1,2,3,4\n";
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
	TTree * outTree = new TTree("tagged","BAND Neutrons and CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	int eventnumber = 0;
	int runnumber = 0;
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
	//	Electron info:
	clashit eHit;
	//	Tagged info:
	TClonesArray * tags = new TClonesArray("taghit");
	TClonesArray &saveTags = *tags;
	//Smeared info
	clashit eHit_smeared;
	//	Smeared Tagged info:
	TClonesArray * tags_smeared = new TClonesArray("taghit");
	TClonesArray &saveTags_smeared = *tags_smeared;
	// 	Charged particles
	int charged_mult        = 0;
        TClonesArray * charged_particles = new TClonesArray("particles");
        TClonesArray &saveCharged_particles = *charged_particles;
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
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);
	//	Tagged branches:
	outTree->Branch("tag"		,&tags			);
	//	MC branches:
	outTree->Branch("genMult"	,&genMult		);
	outTree->Branch("mcParts"	,&mcParts		);
	if( MC_DATA_OPT == 0 ){ // if this is a MC file, define smeared branches
		//	Smeared Electron branches:
		outTree->Branch("eHit_smeared"		,&eHit_smeared			);
		//	Smeared Tagged branches:
		outTree->Branch("tag_smeared"		,&tags_smeared			);
	}
	//	Charged particle branches:
	outTree->Branch("charged_mult"		,&charged_mult		);
	outTree->Branch("charged_particles"	,&charged_particles	);

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Load the electron PID class:
	e_pid ePID;
	// Load the DC fiducial class for electrons;
	DCFiducial DCfid_electrons;

	// Load the QADB
	QADB * qa = new QADB();

	// Load input file
	for( int i = 5 ; i < argc ; i++ ){
		if( MC_DATA_OPT == 0 || MC_DATA_OPT == 2){
			int runNum = 11;
			Runno = runNum;
			// Set the initial Ebeam value so that it can be used for the PID class
			if( PERIOD == 0 ) Ebeam = 10.599; // from RCDB: 10598.6
			if( PERIOD == 1 ) Ebeam = 10.200; // from RCDB: 10199.8
			if( PERIOD == 2 ) Ebeam = 10.389; // from RCDB: 10389.4
			if( PERIOD == 3 ) Ebeam = 4.247;  // current QE-MC value, RCDB value: 4171.79. Is about ~1.018 wrong due to issues with magnet settings
			if( PERIOD == 4 ) Ebeam = 2.07; 
		}
		else if( MC_DATA_OPT == 1){
			int runNum = getRunNumber(argv[i]);
			Runno = runNum;
			auto cnd = connection.GetCondition(runNum, "beam_energy");
			// For data, get Ebeam from RCDB:
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
		reader.readDictionary(factory);
		hipo::bank	event_info		(factory.getSchema("REC::Event"		));
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::bank  	run_config 		(factory.getSchema("RUN::config"));
		hipo::bank      DC_Track                (factory.getSchema("REC::Track"         ));
		hipo::bank      DC_Traj                 (factory.getSchema("REC::Traj"          ));
		hipo::bank	cherenkov		(factory.getSchema("REC::Cherenkov"	));
		hipo::event 	readevent;
		hipo::bank	band_rawhits		(factory.getSchema("BAND::rawhits"	));
		hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
		BParticle	clas_particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		hipo::bank	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank	mc_event_info		(factory.getSchema("MC::Event"		));
		hipo::bank	mc_particle		(factory.getSchema("MC::Particle"	));


		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		double torussetting = 0;
		while(reader.next()==true){
			BAND->Clear();
			// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			eventnumber = 0;
			runnumber = 0;
			// Neutron
			nMult		= 0;
			passed		= 0;
			nleadindex = -1;
			goodneutron = false;
			bandhit nHit[maxNeutrons];
			nHits->Clear("C");
			saveHit.Clear("C");
			// Tag
			taghit tag[maxNeutrons];
			taghit tag_smeared[maxNeutrons];
			tags->Clear("C");
			saveTags.Clear("C");
			// Electron
			eHit.Clear();
			// MC
			genMult = 0;
			weight = 1;
			genpart mcPart[maxGens];
			mcParts->Clear("C");
			saveMC.Clear("C");
			if( MC_DATA_OPT == 0 ){ // if this is a MC file, clear smeared and input branches
				//Clear output smear branches
				eHit_smeared.Clear();
				tags_smeared->Clear();
			}

			// Charged particles
			particles charged_particle[maxParticles];
			charged_mult	 = 0;
			charged_particles->Clear();

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if( event_counter > 30 ) break;
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
			// electron struct
			readevent.getStructure(clas_particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(DC_Track);
			readevent.getStructure(DC_Traj);
			readevent.getStructure(cherenkov);
			// monte carlo struct
			readevent.getStructure(mc_event_info);
			readevent.getStructure(mc_particle);

			if( event_counter == 1 ){
				//cout << Runno << "\n";
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
				else if (Runno >= 15566 && Runno < 15700) { //RGM deuterium 2.1 runs
					period = 4;
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


			//Get Event number from RUN::config
			eventnumber = run_config.getInt( 1 , 0 );
			runnumber = run_config.getInt( 0 , 0 );
			if( runnumber == 0 ) continue; // skip empty header events
			if( runnumber != Runno ){ cerr << "Run number mistmatch! Exiting\n"; exit(-1); }

			// Do QADB cut
			if(MC_DATA_OPT == 1 && (PERIOD != 3 && PERIOD != 4)) {
				if( !qa->Golden(runnumber,eventnumber) ) continue;
			}

			//from first event get RUN::config torus Setting
		 	// inbending = negative torussetting, outbending = torusseting
			torussetting = run_config.getFloat( 7 , 0 );


			// For simulated events, get the weight for the event
			if( MC_DATA_OPT == 0 || MC_DATA_OPT == 2){
				// For MC, the above Ebeam has already been set, but we will reset it based on 
				// any new information in the header file:
				getMcInfo( mc_particle , mc_event_info , mcPart , starttime, weight, Ebeam , genMult );
			}

			//for option 2 just fill the generated infos in the output file, skip remaining part of event loop
			if( MC_DATA_OPT == 2){
				// Store the mc particles in TClonesArray
				for( int n = 0 ; n < maxGens ; n++ ){
					new(saveMC[n]) genpart;
					saveMC[n] = &mcPart[n];
				}
				outTree->Fill();
				continue;
			}

				//Till the end of the for loop is only executed for MC_DATA_OPT = 0 and 1
			// Get integrated charge, livetime and start-time from REC::Event
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );




			// Grab the electron information:
			getElectronInfo( clas_particles , calorimeter , scintillator , DC_Track, DC_Traj, cherenkov, 0, eHit , starttime , Runno , Ebeam );

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

			// Grab the neutron information:
				// Form the PMTs and Bars for BAND:
			BAND->createPMTs( &band_adc, &band_tdc, &run_config );
			BAND->createBars();
			//BAND->storeHits( nMult , nHit , starttime , eHit.getVtz() ); // use event-by-event electron vertex for pathlength-z
			BAND->storeHits( nMult , nHit , starttime , BAND->getRGBVertexOffset() ); // use average z-offset of the target for pathlength-z

			// Create the tagged information if we have neutrons appropriately aligned in time:
			getTaggedInfo(	eHit	,  nHit	 ,  tag  , Ebeam , nMult );

			//MC smearing
			if( MC_DATA_OPT == 0 ){ // if this is a MC file, do smearing and add values

				// Grab the electron information for the smeared eHit Object
				getElectronInfo( clas_particles , calorimeter , scintillator , DC_Track, DC_Traj, cherenkov, 0, eHit_smeared , starttime , Runno , Ebeam );

				//read electron vector
				TVector3 reco_electron(0,0,0);
				reco_electron.SetMagThetaPhi(eHit.getMomentum(),eHit.getTheta(),eHit.getPhi());
				//Smear Reconstructed electron in Momentum, Theta and Phi
				smearRGA(reco_electron);

				//Recalculate Electron Kinematics with smeared values
				recalculate_clashit_kinematics(eHit_smeared, Ebeam, reco_electron);

				// Create the tagged smeared information from the smeared electron and neutron information:
				getTaggedInfo(	eHit_smeared	,  nHit	 ,  tag_smeared , Ebeam , nMult );

				for( int n = 0 ; n < nMult ; n++ ){
					new(saveTags_smeared[n]) taghit;
					saveTags_smeared[n] = &tag_smeared[n];
					}
			}
			
			 // Grab the information for other charged particle:
                        getParticleInfo( clas_particles, charged_particle, scintillator, charged_mult );

                        // Store the charged particles in TClonesArray
                        for( int n = 0 ; n < charged_mult ; n++ ){
                                new(saveCharged_particles[n]) particles;
                                saveCharged_particles[n] = &charged_particle[n];
                        }

			// Store the neutrons in TClonesArray
			for( int n = 0 ; n < nMult ; n++ ){
				new(saveHit[n]) bandhit;
				saveHit[n] = &nHit[n];
			}
			for( int n = 0 ; n < nMult ; n++ ){
				new(saveTags[n]) taghit;
				saveTags[n] = &tag[n];
			}
			// Store the mc particles in TClonesArray
			for( int n = 0 ; n < maxGens ; n++ ){
				new(saveMC[n]) genpart;
				saveMC[n] = &mcPart[n];
			}

			//if (nMult == 1) {
			//	cout << nMult << " " << 0 << "\n";
			//	if() // layer 6 check
			//	goodneutron =  true;
			//	passed = 1;
			//	nleadindex = 0;
			//}
			//If nMult > 1: Take nHit and check if good event and give back leading hit index and boolean
			if (nMult > 0) {
				//cout << nMult << " ";
				//pass Nhit array, multiplicity and reference to leadindex which will be modified by function
				goodneutron = goodNeutronEvent(nHit, nMult, nleadindex, MC_DATA_OPT, passed );
			}
			
			// Fill tree based on d(e,e'n)X for data
			if( ( (nMult > 0 && goodneutron) ) && MC_DATA_OPT == 1 ){
				qa->AccumulateCharge();
				outTree->Fill();
			} // else fill tree on d(e,e')nX for MC
			else if( MC_DATA_OPT == 0 ||  MC_DATA_OPT == 2 ){
				outTree->Fill();
			}



		} // end loop over events
		cout << "Charge analyzed in file: " << argv[i] << " " << qa->GetAccumulatedCharge() << " [nC]\n";
	}// end loop over files

	TVector3 charged_analyzed;
	charged_analyzed.SetXYZ( qa->GetAccumulatedCharge() , 0 , 0 );

	outFile->cd();
	outTree->Write();
	charged_analyzed.Write("charge_analyzed");
	outFile->Close();

	return 0;
}
