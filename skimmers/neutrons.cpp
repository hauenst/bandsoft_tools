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

#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 4 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [LoadShiftsOpt] [inputFile]\n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[LoadShiftsOpt] = 0 (don't load) 1 (load from include)\n";
		cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
		return -1;
	}
	int loadshifts_opt = atoi(argv[2]);

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
	// 	Neutron info:
	int nMult		= 0;
	bandhit nHit[maxNeutrons];
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("nHit"		,nHit			);
	
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
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;
		hipo::bank	band_rawhits		(factory.getSchema("BAND::rawhits"	));
		hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
		
		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		while(reader.next()==true){
			// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			nMult		= 0;
			for( int thishit = 0; thishit < maxNeutrons ; thishit++){
				nHit[thishit].Clear();
			}

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			event_counter++;
			//if( event_counter > 100000 ) break;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			// band struct
			readevent.getStructure(band_hits);
			readevent.getStructure(band_rawhits);
			readevent.getStructure(band_adc);
			readevent.getStructure(band_tdc);
	
			// Get integrated charge, livetime and start-time from REC::Event
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );
			
			// Grab the neutron information:
			TVector3 nMomentum[maxNeutrons], nPath[maxNeutrons];
			getNeutronInfo( band_hits, band_rawhits, band_adc, band_tdc, nMult, nHit , starttime , runNum);
			if( loadshifts_opt ){
				for( int n = 0 ; n < nMult ; n++ ){
					nHit[n].setTofFadc(	nHit[n].getTofFadc() - FADC_INITBAR[(int)nHit[n].getBarID()] - FADC_INITRUN[Runno] );
					//nHit[n].setTof(	nHit[n].getTof() - TDC_INITBAR[(int)nHit[n].getBarID()] - TDC_INITRUN[Runno] );
				}
			}

			// Fill tree based on d(e,e'n)X
			if( nMult != 0 ) outTree->Fill();

		} // end loop over events
	}// end loop over files
	
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}



