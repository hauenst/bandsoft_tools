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

//from readhipo_helper maxParticles	= 100;
//from readhipo_helper maxScinHits = 100;
//from readhipo_helper maxNeutrons = 200;
//TODO:  Check for memory leaks(memset??). IMplement Fiducial cuts. Update RUN::config read to BConfig read


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
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");

	//	Event info:
	int Runno = 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	bool goodneutron = false;
	int nleadindex = -1;
	int eventnumber = 0;

	// 	Neutron info:
	int nMult		= 0;
	TClonesArray * nHits = new TClonesArray("bandhit");
	TClonesArray &saveHit = *nHits;
	//	Electron info:
	clashit eHit;




	// 	Positive Particles info:
	int pMult		= 0;
	int pIndex		[maxParticles]= {0};
	double pPid		[maxParticles]= {0.};
	double pCharge		[maxParticles]= {0.};
	double pStatus		[maxParticles]= {0.};
	double pTime		[maxParticles]= {0.};
	double pBeta		[maxParticles]= {0.};
	double pChi2pid		[maxParticles]= {0.};
	double p_vtx		[maxParticles]= {0.};
	double p_vty		[maxParticles]= {0.};
	double p_vtz		[maxParticles]= {0.};
	double p_p		[maxParticles]= {0.};
	double theta_p		[maxParticles]= {0.};
	double phi_p		[maxParticles]= {0.};


 // Information from REC::Scintillator for positive Particles
  int scinHits = 0;
	double hit_pindex [maxScinHits]= {0.}; //pPid of associated positive particle
	double hit_detid [maxScinHits]= {0.};
	double hit_energy [maxScinHits]= {0.};
	double hit_time [maxScinHits]= {0.};
	double hit_x [maxScinHits]= {0.};
	double hit_y [maxScinHits]= {0.};
	double hit_z [maxScinHits]= {0.};
	double hit_path [maxScinHits]= {0.};
	double hit_status [maxScinHits]= {0.};

	//	MC info:
	double weight		= 0;
	int genMult		= 0;
	TClonesArray * mcParts = new TClonesArray("genpart");
	TClonesArray &saveMC = *mcParts;


	//Branches:
	//Event
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("eventnumber",&eventnumber);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("nHits"		,&nHits			);
	//Branches to store if good Neutron event and leadindex
	outTree->Branch("goodneutron"		,&goodneutron	);
	outTree->Branch("nleadindex"		,&nleadindex			);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);
	//Positive Particles
	outTree->Branch("pMult"		,&pMult			);
	outTree->Branch("pIndex"		,&pIndex			,"pIndex[pMult]/I"	);
	outTree->Branch("pPid"		,&pPid			,"pPid[pMult]/D"	);
	outTree->Branch("pCharge"	,&pCharge		,"pCharge[pMult]/D"	);
	outTree->Branch("pStatus"	,&pStatus		,"pStatus[pMult]/D"	);
	outTree->Branch("pTime"		,&pTime			,"pTime[pMult]/D"	);
	outTree->Branch("pBeta"		,&pBeta			,"pBeta[pMult]/D"	);
	outTree->Branch("pChi2pid",&pChi2pid		,"pChi2pid[pMult]/D"	);
	outTree->Branch("p_vtx"		,&p_vtx			,"p_vtx[pMult]/D"	);
	outTree->Branch("p_vty"		,&p_vty			,"p_vty[pMult]/D"	);
	outTree->Branch("p_vtz"		,&p_vtz			,"p_vtz[pMult]/D"	);
	outTree->Branch("p_p"		,&p_p			,"p_p[pMult]/D"		);
	outTree->Branch("theta_p"	,&theta_p		,"theta_p[pMult]/D"	);
	outTree->Branch("phi_p"		,&phi_p			,"phi_p[pMult]/D"	);


 //Scintillators
  outTree->Branch("scinHits"		,&scinHits		);
  outTree->Branch("hit_pindex"	,&hit_pindex		,"hit_pindex[scinHits]/D"	);
  outTree->Branch("hit_detid"	,&hit_detid		,"hit_detid[scinHits]/D"	);
  outTree->Branch("hit_energy"	,&hit_energy		,"hit_energy[scinHits]/D"	);
	outTree->Branch("hit_time"	,&hit_time		,"hit_time[scinHits]/D"	);
	outTree->Branch("hit_x"	,&hit_x		,"hit_x[scinHits]/D"	);
  outTree->Branch("hit_y"	,&hit_y		,"hit_y[scinHits]/D"	);
	outTree->Branch("hit_z"	,&hit_z		,"hit_z[scinHits]/D"	);
	outTree->Branch("hit_path"	,&hit_path		,"hit_path[scinHits]/D"	);
	outTree->Branch("hit_status"	,&hit_status		,"hit_status[scinHits]/D"	);
	if( MC_DATA_OPT == 0){
	//	MC branches:
		outTree->Branch("genMult"	,&genMult		);
		outTree->Branch("mcParts"	,&mcParts		);
	}

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	//Load Bar shifts
	//TODO: Make shifts flexible to use
	shiftsReader shifts;
	double * FADC_BARSHIFTS_LER;
	double * TDC_BARSHIFTS_LER;
	double * FADC_BARSHIFTS_SPRING19;
	double * TDC_BARSHIFTS_SPRING19;
	if( loadshifts_opt ){
		shifts.LoadInitBarFadc("../include/LER_FADC_shifts.txt");
		FADC_BARSHIFTS_LER = (double*) shifts.getInitBarFadc();
		shifts.LoadInitBar("../include/LER_TDC_shifts.txt");
		TDC_BARSHIFTS_LER = (double*) shifts.getInitBar();
		shifts.LoadInitBarFadc	("../include/FADC_pass1v0_initbar.txt");
		FADC_BARSHIFTS_SPRING19 = (double*) shifts.getInitBarFadc();
		shifts.LoadInitBar	("../include/TDC_pass1v0_initbar.txt");
		TDC_BARSHIFTS_SPRING19 = (double*) shifts.getInitBar();
	}
	/*
		double * FADC_INITRUN;
		// Load run-by-run shifts
		shifts.LoadInitRunFadc("../include/FADC_pass1v0_initrun.txt");
		FADC_INITRUN = (double*) shifts.getInitRunFadc();

*/
// Effective velocity for re-doing x- calculation
		double * FADC_EFFVEL_S6200;
		double *  TDC_EFFVEL_S6200;
		double * FADC_EFFVEL_S6291;
		double *  TDC_EFFVEL_S6291;
		double *  FADC_LROFF_S6200;
		double *   TDC_LROFF_S6200;
		double *  FADC_LROFF_S6291;
		double *   TDC_LROFF_S6291;
		shifts.LoadEffVel	("../include/EffVelocities_S6200.txt",	"../include/EffVelocities_S6291.txt");
		shifts.LoadLrOff	("../include/LrOffsets_S6200.txt",	"../include/LrOffsets_S6291.txt");
		FADC_EFFVEL_S6200	= (double*) shifts.getFadcEffVel(6200);
		TDC_EFFVEL_S6200	= (double*)  shifts.getTdcEffVel(6200);
		FADC_EFFVEL_S6291	= (double*) shifts.getFadcEffVel(6291);
		TDC_EFFVEL_S6291	= (double*)  shifts.getTdcEffVel(6291);

		FADC_LROFF_S6200	= (double*) shifts.getFadcLrOff(6200);
		TDC_LROFF_S6200		= (double*)  shifts.getTdcLrOff(6200);
		FADC_LROFF_S6291	= (double*) shifts.getFadcLrOff(6291);
		TDC_LROFF_S6291		= (double*)  shifts.getTdcLrOff(6291);


	// Load input file
	for( int i = 4 ; i < argc ; i++ ){
		if( MC_DATA_OPT == 0){
			int runNum = 11;
			Runno = runNum;
			//Ebeam is calculated later via MC readin
		}
		else if( MC_DATA_OPT == 1){ //Data
			// Using run number of current file, grab the beam energy from RCDB
			int runNum = getRunNumber(argv[i]);
			Runno = runNum;
			auto cnd = connection.GetCondition(runNum, "beam_energy");
			Ebeam = cnd->ToDouble() / 1000.;// [GeV] -- conversion factor
			if (runNum >= 11286 && runNum < 11304)
			{
				Ebeam *= 1.018; //fudge factor for Low energy run due to miscalibration in RCDB
			}
			current = connection.GetCondition( runNum, "beam_current") ->ToDouble(); // [nA]
		}
		else{
			cout << "Wrong option for MC_DATA. Option is " << MC_DATA_OPT  << ".Exit program" << endl;
			exit(-1);
		}



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
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	band_rawhits		(factory.getSchema("BAND::rawhits"	));
		hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::bank  run_config (factory.getSchema("RUN::config"));
		hipo::bank	mc_event_info		(factory.getSchema("MC::Event"		));
		hipo::bank	mc_particle		(factory.getSchema("MC::Particle"	));
		hipo::event 	readevent;

		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		int run_number_from_run_config = 0;

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
			// MC
			genMult = 0;
			weight = 1;
			genpart mcPart[maxGens];
			mcParts->Clear();


			pMult		= 0;
			memset(	pIndex		,0	,sizeof(pIndex		)	);
			memset(	pPid		,0	,sizeof(pPid		)	);
			memset(	pCharge		,0	,sizeof(pCharge		)	);
			memset(	pStatus		,0	,sizeof(pStatus		)	);
			memset(	pTime		,0	,sizeof(pTime		)	);
			memset(	pBeta		,0	,sizeof(pBeta		)	);
			memset(	pChi2pid	,0	,sizeof(pChi2pid	)	);
			memset(	p_vtx		,0	,sizeof(p_vtx		)	);
			memset(	p_vty		,0	,sizeof(p_vty		)	);
			memset(	p_vtz		,0	,sizeof(p_vtz		)	);
			memset(	p_p		,0	,sizeof(p_p		)	);
			memset(	theta_p		,0	,sizeof(theta_p		)	);
			memset(	phi_p		,0	,sizeof(phi_p		)	);

			scinHits = 0;
			memset(	hit_pindex	,0	,sizeof(hit_pindex		)	);
			memset(	hit_detid		,0	,sizeof(hit_detid		)	);
			memset(	hit_energy		,0	,sizeof(hit_energy		)	);
			memset( hit_time	,0	,sizeof(hit_time	)	);
			memset( hit_x	,0	,sizeof(hit_x	)	);
			memset( hit_y	,0	,sizeof(hit_y	)	);
			memset( hit_z	,0	,sizeof(hit_z	)	);
			memset( hit_path	,0	,sizeof(hit_path	)	);
			memset( hit_status	,0	,sizeof(hit_status	)	);



			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if( event_counter > 100 ) continue;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
			readevent.getStructure(run_config);
			//Electrons and Positive particles
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			//BAND
			readevent.getStructure(band_hits);
			readevent.getStructure(band_rawhits);
			readevent.getStructure(band_adc);
			readevent.getStructure(band_tdc);
			readevent.getStructure(DC_Track);
			readevent.getStructure(DC_Traj);
			// monte carlo struct
			readevent.getStructure(mc_event_info);
			readevent.getStructure(mc_particle);

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
			//if (event_counter < 100) {
			//	cout << "event number " << eventnumber << " , runnumebr " << run_number_from_run_config << endl;
			//}

			// For simulated events, get the weight for the event, also sets the beam energy for MC. Ebeam is reference!
			//getMcInfo has to be called before getElectronInfo because of Ebeam value
			if( MC_DATA_OPT == 0){
				getMcInfo( mc_particle , mc_event_info , mcPart , starttime, weight, Ebeam , genMult );
			}


			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , DC_Track, DC_Traj, eHit , starttime , Runno , Ebeam );

			//	Do electron PID cuts
			//		only PID (11) and charge (-1) selection on first particle
			if( eHit.getPID() != 11 					) continue;
			if( eHit.getCharge() != -1					) continue;

			//check other particles for electron or negative charge
			/*for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
				if (particles.getPid(row)==11) { //check for other electron (electron exclusivity) and skip event
					ePass = false;
				}
				if (particles.getCharge(row)==-1) { //check for other negative charge particle and skip event
					ePass = false;
				}
			}
			if( !ePass ) continue;*/


			// Grab the information for other charged particle:
			TVector3 pVertex[maxParticles], pMomentum[maxParticles];
			getParticleInfo( particles, pPid, pMomentum, pVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus, pIndex, pMult);
			//Fill the information for other charged particles
			for( int p = 0 ; p < pMult ; p++ ){
				p_vtx[p]		= pVertex[p].X();
				p_vty[p]		= pVertex[p].Y();
				p_vtz[p]		= pVertex[p].Z();
				p_p[p]			= pMomentum[p].Mag();
				theta_p[p]	= pMomentum[p].Theta();
				phi_p[p]		= pMomentum[p].Phi();

			}

			TVector3 hitVector[maxScinHits];
			getScinHits( scintillator, hit_pindex, hit_detid, hit_energy, hit_time, hitVector, hit_path, hit_status, pIndex, pMult, scinHits);
			for( int hit = 0 ; hit < scinHits ; hit++ ){
				  hit_x[hit] = hitVector[hit].X();
				  hit_y[hit] = hitVector[hit].Y();
				  hit_z[hit] = hitVector[hit].Z();
			}

			if( MC_DATA_OPT == 0 ){
				getNeutronInfo( band_hits, band_rawhits, band_adc, band_tdc, nMult, nHit , starttime , Runno );
			}
			else{
				getNeutronInfo( band_hits, band_rawhits, band_adc, band_tdc, nMult, nHit , starttime , Runno,
						1, 	FADC_LROFF_S6200,	TDC_LROFF_S6200,
							FADC_LROFF_S6291,	TDC_LROFF_S6291,
							FADC_EFFVEL_S6200,	TDC_EFFVEL_S6200,
							FADC_EFFVEL_S6291,	TDC_EFFVEL_S6291	);
			}

			if( loadshifts_opt && MC_DATA_OPT !=0){
				//Load of shifts depending on run number
				if (Runno >= 11286 && Runno < 11304)	{
					//LER corrections
					for( int n = 0 ; n < nMult ; n++ ){
						nHit[n].setTofFadc(	nHit[n].getTofFadc() 	- FADC_BARSHIFTS_LER[(int)nHit[n].getBarID()] );
						nHit[n].setTof(		nHit[n].getTof() 	- TDC_BARSHIFTS_LER[(int)nHit[n].getBarID()]  );
					}
				}
				else if (Runno > 6100 && Runno < 6800) { //Spring 19 data
					for( int n = 0 ; n < nMult ; n++ ){
						nHit[n].setTofFadc(	nHit[n].getTofFadc() 	- FADC_BARSHIFTS_SPRING19[(int)nHit[n].getBarID()] );
						nHit[n].setTof(		nHit[n].getTof() 	- TDC_BARSHIFTS_SPRING19[(int)nHit[n].getBarID()]  );
					}
				}
			}


			// Store the neutrons in TClonesArray
			for( int n = 0 ; n < nMult ; n++ ){
				new(saveHit[n]) bandhit;
				saveHit[n] = &nHit[n];
			}
			if( MC_DATA_OPT == 0){
			// Store the mc particles in TClonesArray
				for( int n = 0 ; n < maxGens ; n++ ){
					new(saveMC[n]) genpart;
					saveMC[n] = &mcPart[n];
				}
			}

			if (nMult == 1) {
				goodneutron =  true;
				nleadindex = 0;
			}
			//If nMult > 1: Take nHit and check if good event and give back leading hit index and boolean
			if (nMult > 1) {
				//pass Nhit array, multiplicity and reference to leadindex which will be modified by function
				goodneutron = goodNeutronEvent(nHit, nMult, nleadindex, 1);
			}

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
