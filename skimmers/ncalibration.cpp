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

using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 5 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[<MC,DATA> = <0, 1> \n";
		cerr << "\t\t[<load shifts N,Y> = <0, 1> \n";
		cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
		return -1;
	}

	int MC_DATA_OPT = atoi(argv[2]);
	int loadshifts_opt = atoi(argv[3]);

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("calib","BAND Neutrons and CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	bool goodneutron 	= false;
	int nleadindex 		= -1;
	// 	Neutron info:
	int nMult		= 0;
	TClonesArray * nHits 	= new TClonesArray("bandhit");
	TClonesArray &saveHit 	= *nHits;
	//	Electron info:
	clashit eHit;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("nHits"		,&nHits			);
	//Branches to store if good Neutron event and leadindex
	outTree->Branch("goodneutron"		,&goodneutron	);
	outTree->Branch("nleadindex"		,&nleadindex			);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	shiftsReader shifts;
	double * FADC_INITBAR;
	double * FADC_INITRUN;
	if( loadshifts_opt ){
		// Load bar shifts
		shifts.LoadInitBarFadc("../include/FADC_pass1v0_initbar.txt");
		FADC_INITBAR = (double*) shifts.getInitBarFadc();
		// Load run-by-run shifts
		shifts.LoadInitRunFadc("../include/FADC_pass1v0_initrun.txt");
		FADC_INITRUN = (double*) shifts.getInitRunFadc();
	}

	// Load input file
	for( int i = 4 ; i < argc ; i++ ){
		if( MC_DATA_OPT == 1){
			int runNum = getRunNumber(argv[i]);
			Runno = runNum;
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
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::bank      DC_Track                (factory.getSchema("REC::Track"         ));
		hipo::bank      DC_Traj                 (factory.getSchema("REC::Traj"          ));
		hipo::event 	readevent;
		hipo::bank	band_rawhits		(factory.getSchema("BAND::rawhits"	));
		hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));


		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		while(reader.next()==true){
			// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			// Neutron
			nMult		= 0;
			nleadindex = -1;
			goodneutron = false;
			bandhit nHit[maxNeutrons];
			nHits->Clear();
			// Electron
			eHit.Clear();


			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			if(event_counter > 1000000 ) break;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			// band struct
			readevent.getStructure(band_hits);
			readevent.getStructure(band_rawhits);
			readevent.getStructure(band_adc);
			readevent.getStructure(band_tdc);
			// electron struct
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(DC_Track);
			readevent.getStructure(DC_Traj);

			// Get integrated charge, livetime and start-time from REC::Event
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			// Skim the event so we only have a single electron and NO other particles:
			int nElectrons = 0;
			int nOthers = 0;
			for( int part = 0 ; part < particles.getRows() ; part++ ){
				int PID 	= particles.getPid(part);
				int charge 	= particles.getCharge(part);
				if( PID == 11 && charge == -1 ) nElectrons++;
				else{ nOthers++; }
			}
			if( nElectrons != 1 ) 	continue;
			//if( nOthers != 0 )	continue;


			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , DC_Track, DC_Traj, eHit , starttime , Runno , Ebeam );
			
			// Further skim the event so that the trigger electron passes basic fiducial requirements
			if( eHit.getPID() != 11 					) continue;
			if( eHit.getCharge() != -1					) continue; 
			if( eHit.getEoP() < 0.17 || eHit.getEoP() > 0.3			) continue;
			if( eHit.getEpcal() < 0.07					) continue;
			if( eHit.getV() < 15 || eHit.getW() < 15			) continue;
			if( eHit.getVtz() < -8 || eHit.getVtz() > 3			) continue;
			if( eHit.getMomentum() < 3 || eHit.getMomentum() > 10.6		) continue;
			if( eHit.getQ2() < 2 || eHit.getQ2() > 10			) continue;
			if( eHit.getW2() < 2*2						) continue;

			// Grab the neutron information:
			getNeutronInfo( band_hits, band_rawhits, band_adc, band_tdc, nMult, nHit , starttime , Runno);
			if( loadshifts_opt ){
				for( int n = 0 ; n < nMult ; n++ ){
					nHit[n].setTofFadc(	nHit[n].getTofFadc() - FADC_INITBAR[(int)nHit[n].getBarID()] - FADC_INITRUN[Runno] );
					//nHit[n].setTof(	nHit[n].getTof() - TDC_INITBAR[(int)nHit[n].getBarID()] - TDC_INITRUN[Runno] );
				}
			}

			// Store the neutrons in TClonesArray
			for( int n = 0 ; n < nMult ; n++ ){
				new(saveHit[n]) bandhit;
				saveHit[n] = &nHit[n];
			}

			if (nMult == 1) {
				goodneutron =  true;
				nleadindex = 0;
			}
			//If nMult > 1: Take nHit and check if good event and give back leading hit index and boolean
			if (nMult > 1) {
				//pass Nhit array, multiplicity and reference to leadindex which will be modified by function
				goodneutron = goodNeutronEvent(nHit, nMult, nleadindex, MC_DATA_OPT);
			}

			// Fill tree based on d(e,e'n)X for data
			if( (nMult == 1 || (nMult > 1 && goodneutron) )&& MC_DATA_OPT == 1 ){
				outTree->Fill();
			} 



		} // end loop over events
	}// end loop over files

	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}
