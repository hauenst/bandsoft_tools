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

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;

const int maxPositive	= 100;
const int maxScinHits = 100;
//from readhipo_helper maxNeutrons = 200;


bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getPositiveInfo( BParticle particles, double pid[maxPositive], TVector3 momentum[maxPositive], TVector3 vertex[maxPositive],
			double time[maxPositive], double charge[maxPositive], double beta[maxPositive], double chi2pid[maxPositive], double status[maxPositive] , int index[maxPositive], int& multiplicity );
void getScinHits( BScintillator scintillator, double pindex[maxScinHits], double detid[maxScinHits], double energy[maxScinHits], double time[maxScinHits],
	    TVector3 posVector[maxScinHits], double path[maxScinHits], double status[maxScinHits], int posIndex[maxPositive], int posMult, int &scinHits);


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 3 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [inputFile]\n\n";
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

	//	Electron info:
	clashit eHit;
	int ePid		= 0;
	int eCharge		= 0;
	int eStatus		= 0;
	double eTime		= 0;
	double eBeta 		= 0;
	double eChi2pid		= 0;
	double E_tot		= 0;
	double E_pcal		= 0;
	double E_ecin		= 0;
	double E_ecout	= 0;
	double t_e		= 0;
	double dL_e		= 0;
	double lU		= 0;
	double lV		= 0;
	double lW		= 0;
	double e_vtx		= 0;
	double e_vty		= 0;
	double e_vtz		= 0;
	double p_e		= 0;
	double theta_e		= 0;
	double phi_e		= 0;
	double q		= 0;
	double theta_q		= 0;
	double phi_q		= 0;
	double nu		= 0;
	double Q2		= 0;
	double xB		= 0;
	double W2		= 0;

	// 	Positive Particles info:
	int pMult		= 0;
	int pIndex		[maxPositive]= {0};
	double pPid		[maxPositive]= {0.};
	double pCharge		[maxPositive]= {0.};
	double pStatus		[maxPositive]= {0.};
	double pTime		[maxPositive]= {0.};
	double pBeta		[maxPositive]= {0.};
	double pChi2pid		[maxPositive]= {0.};
	double p_vtx		[maxPositive]= {0.};
	double p_vty		[maxPositive]= {0.};
	double p_vtz		[maxPositive]= {0.};
	double p_p		[maxPositive]= {0.};
	double theta_p		[maxPositive]= {0.};
	double phi_p		[maxPositive]= {0.};
	double theta_pq		[maxPositive]= {0.};

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



	// 	Neutron info:
	int nMult		= 0;
	int barID		[maxNeutrons]= {0};
	int nSector [maxNeutrons]= {0};
	int nLayer [maxNeutrons]= {0};
	int nComponent [maxNeutrons]= {0};
	double dL_n		[maxNeutrons]= {0.};
	double theta_n		[maxNeutrons]= {0.};
	double phi_n		[maxNeutrons]= {0.};
	double nTof		[maxNeutrons]= {0.};
	double nTofFADC		[maxNeutrons]= {0.};
	double nEdep		[maxNeutrons]= {0.};
	double n_bandx [maxNeutrons] =  {0.};
	double n_bandy [maxNeutrons] =  {0.};
	double n_bandz [maxNeutrons] =  {0.};
	double nADCleft [maxNeutrons]= {0.};
	double nADCright [maxNeutrons]= {0.};
	int nStatus [maxNeutrons]= {0};

	//Branches:
	//Event
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	//Electron
	outTree->Branch("ePid"		,&ePid			);
	outTree->Branch("eCharge"	,&eCharge		);
	outTree->Branch("eStatus"	,&eStatus		);
	outTree->Branch("eTime"		,&eTime			);
	outTree->Branch("eBeta"		,&eBeta 		);
	outTree->Branch("eChi2pid"	,&eChi2pid		);
	outTree->Branch("E_tot"		,&E_tot			);
	outTree->Branch("E_pcal"	,&E_pcal		);
	outTree->Branch("E_ecin"	,&E_ecin		);
	outTree->Branch("E_ecout"	,&E_ecout		);
	outTree->Branch("t_e"		,&t_e			);
	outTree->Branch("dL_e"		,&dL_e			);
	outTree->Branch("lU"		,&lU			);
	outTree->Branch("lV"		,&lV			);
	outTree->Branch("lW"		,&lW			);
	outTree->Branch("e_vtx"		,&e_vtx			);
	outTree->Branch("e_vty"		,&e_vty			);
	outTree->Branch("e_vtz"		,&e_vtz			);
	outTree->Branch("p_e"		,&p_e			);
	outTree->Branch("theta_e"	,&theta_e		);
	outTree->Branch("phi_e"		,&phi_e			);
	outTree->Branch("q"		,&q			);
	outTree->Branch("theta_q"	,&theta_q		);
	outTree->Branch("phi_q"		,&phi_q			);
	outTree->Branch("nu"		,&nu			);
	outTree->Branch("Q2"		,&Q2			);
	outTree->Branch("xB"		,&xB			);
	outTree->Branch("W2"		,&W2			);
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
	outTree->Branch("theta_pq"	,&theta_pq		,"theta_pq[pMult]/D"	);

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

//BAND Neutrons
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("barID"		,&barID			,"barID[nMult]/I"	);
	outTree->Branch("nSector" ,&nSector   ,"nSector[nMult]/I");
	outTree->Branch("nLayer"  ,&nLayer    ,"nLayer[nMult]/I");
	outTree->Branch("nComponent" ,&nComponent  ,"nComponent[nMult]/I");
	outTree->Branch("dL_n"		,&dL_n			,"dL_n[nMult]/D"	);
	outTree->Branch("theta_n"	,&theta_n		,"theta_n[nMult]/D"	);
	outTree->Branch("phi_n"		,&phi_n			,"phi_n[nMult]/D"	);
	outTree->Branch("nTof"		,&nTof		,"nTof[nMult]/D"	);
	outTree->Branch("nTofFADC"		,&nTofFADC			,"nTofFADC[nMult]/D"	);
	outTree->Branch("nEdep"		,&nEdep			,"nEdep[nMult]/D"	);
	outTree->Branch("n_bandx"		,&n_bandx			,"n_bandx[nMult]/D"	);
	outTree->Branch("n_bandy"		,&n_bandy			,"n_bandy[nMult]/D"	);
	outTree->Branch("n_bandz"		,&n_bandz			,"n_bandz[nMult]/D"	);
	outTree->Branch("nADCleft"	,&nADCleft		,"nADCleft[nMult]/D"	);
	outTree->Branch("nADCright"	,&nADCright		,"nADCright[nMult]/D"	);
	outTree->Branch("nStatus"		,&nStatus			,"nStatus[nMult]/I"	);


	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	//Load Bar shifts
	//TODO: Update name of files F.H. 25/09/2020
	shiftsReader shifts;
	double * FADC_BARSHIFTS;
	double * TDC_BARSHIFTS;
	shifts.LoadInitBarFadc("../include/LER_FADC_shifts.txt");
	FADC_BARSHIFTS = (double*) shifts.getInitBarFadc();
	shifts.LoadInitBar("../include/LER_TDC_shifts.txt");
	TDC_BARSHIFTS = (double*) shifts.getInitBar();

	/*
		double * FADC_INITRUN;
		// Load run-by-run shifts
		shifts.LoadInitRunFadc("../include/FADC_pass1v0_initrun.txt");
		FADC_INITRUN = (double*) shifts.getInitRunFadc();

*/



	// Load input file
	for( int i = 2 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		int runNum = getRunNumber(argv[i]);
		Runno = runNum;
		auto cnd = connection.GetCondition(runNum, "beam_energy");
		Ebeam = cnd->ToDouble() / 1000. * 1.018; // [GeV] -- conversion factor due to miscalibration in RCDB

		// Setup hipo reading for this file
		TString inputFile = argv[i];
		hipo::reader reader;
		reader.open(inputFile);
		hipo::dictionary  factory;
		hipo::schema	  schema;
		reader.readDictionary(factory);
		BEvent		event_info		(factory.getSchema("REC::Event"		));
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank      DC_Track                (factory.getSchema("REC::Track"         ));
		hipo::bank      DC_Traj                 (factory.getSchema("REC::Traj"          ));
		//new BAND banks in file with new cook F.H. 28/09/2020
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	band_rawhits		(factory.getSchema("BAND::rawhits"	));
		hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;

		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		while(reader.next()==true){
			// Clear all branches
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			eHit.Clear();
			ePid		= 0;
			eCharge		= 0;
			eStatus		= 0;
			eTime		= 0;
			eBeta 		= 0;
			eChi2pid	= 0;
			E_tot		= 0;
			E_pcal		= 0;
			t_e		= 0;
			dL_e		= 0;
			lU		= 0;
			lV		= 0;
			lW		= 0;
			e_vtx		= 0;
			e_vty		= 0;
			e_vtz		= 0;
			p_e		= 0;
			theta_e		= 0;
			phi_e		= 0;
			q		= 0;
			theta_q		= 0;
			phi_q		= 0;
			nu		= 0;
			Q2		= 0;
			xB		= 0;
			W2		= 0;
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
			memset( theta_pq	,0	,sizeof(theta_pq	)	);

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

			bandhit nHit[maxNeutrons];
			nMult		= 0;
			memset( barID		,0	,sizeof(barID		)	);
			memset( dL_n		,0	,sizeof(dL_n		)	);
			memset( theta_n		,0	,sizeof(theta_n		)	);
			memset( phi_n		,0	,sizeof(phi_n		)	);
			memset( nTof		,0	,sizeof(nTof		)	);
			memset( nTofFADC		,0	,sizeof(nTofFADC		)	);
			memset( nEdep		,0	,sizeof(nEdep		)	);
			memset( n_bandx		,0	,sizeof(n_bandx		)	);
			memset( n_bandy		,0	,sizeof(n_bandy		)	);
			memset( n_bandz		,0	,sizeof(n_bandz		)	);

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if( event_counter > 100000 ) break;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(scaler);
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

			// Get integrated charge, livetime and start-time from REC::Event
			//Currently, REC::Event has uncalibrated livetime / charge, so these will have to work
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			// Get electron from particle bank REC::Particle
			TVector3 eVertex, eMomentum;

			// Grab the electron information:
			getElectronInfo( particles , calorimeter , scintillator , DC_Track, DC_Traj, eHit , starttime , Runno , Ebeam );

			ePid 			= eHit.getPID();
			eTime 		= eHit.getTime();
			eCharge 	= eHit.getCharge();
			eBeta			= eHit.getBeta();
			eChi2pid	= eHit.getChi2();
			eStatus		= eHit.getStatus();
			e_vtx		  = eHit.getVtx();
			e_vty 	  = eHit.getVty();
			e_vtz 		= eHit.getVtz();
			eVertex.SetXYZ(e_vtx,e_vty,e_vtz);
			p_e				= eHit.getMomentum();
			theta_e		= eHit.getTheta();
			phi_e			= eHit.getPhi();
			eMomentum.SetMagThetaPhi( p_e, theta_e, phi_e);

			//	get electron information from scint and calo banks:
			t_e 	= eHit.getTimeScint();
			dL_e	= eHit.getPathScint();
			E_tot 	= eHit.getEtot();
			E_pcal 	= eHit.getEpcal();
			E_ecin  = calorimeter.getECinE(0);
			E_ecout = calorimeter.getECoutE(0);
			lU	= eHit.getU();
			lV	= eHit.getV();
			lW	= eHit.getW();
			//	Do electron PID cuts
			//		only PID (11) and charge (-1) selection on first particle
			bool ePass = checkElectron( ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus , lV , lW , E_tot );
			//check other particles for electron or negative charge
			for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
				if (particles.getPid(row)==11) { //check for other electron (electron exclusivity) and skip event
					ePass = false;
				}
				if (particles.getCharge(row)==-1) { //check for other negative charge particle and skip event
					ePass = false;
				}
			}
			if( !ePass ) continue;

			//read remaining kinematic variables
			q				= eHit.getQ(); //q-vector magnitude
			theta_q	= eHit.getThetaQ();
			phi_q		= eHit.getPhiQ();
			nu			= eHit.getOmega();
			Q2			= eHit.getQ2();
			xB			= eHit.getXb();
			W2			= eHit.getW2();

			//Initalize q-vector from eHit
			TVector3 qMomentum;
			qMomentum.SetMagThetaPhi(q	, theta_q, phi_q);

			// Grab the information for a positive particle:
			TVector3 pVertex[maxPositive], pMomentum[maxPositive];
			getPositiveInfo( particles, pPid, pMomentum, pVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus, pIndex, pMult);
			//Fill the information for all positive particles
			for( int p = 0 ; p < pMult ; p++ ){
				p_vtx[p]		= pVertex[p].X();
				p_vty[p]		= pVertex[p].Y();
				p_vtz[p]		= pVertex[p].Z();
				p_p[p]			= pMomentum[p].Mag();
				theta_p[p]	= pMomentum[p].Theta();
				phi_p[p]		= pMomentum[p].Phi();
				theta_pq[p]	= qMomentum.Angle(pMomentum[p]);

	//			double E_p = sqrt( p_p[p]*p_p[p] + mP*mP );
			//	m_miss[p] 	= sqrt( pow( nu + mD - E_p , 2 ) - ( q*q + p_p[p]*p_p[p] - 2*q*p_p[p]*cos(theta_pq[p]) ) );
		//		if( m_miss[p] != m_miss[p] ) m_miss[p] = 0.;

	//			double W_primeSq = mD*mD - Q2 + mP*mP + 2.*mD*(nu-E_p) - 2.*nu*E_p + 2.*q*p_p[p]*cos(theta_pq[p]);
	//			Wp[p] = sqrt(W_primeSq);
	//			if( Wp[p] != Wp[p] ) Wp[p] = 0.;
			}

			TVector3 hitVector[maxScinHits];
			getScinHits( scintillator, hit_pindex, hit_detid, hit_energy, hit_time, hitVector, hit_path, hit_status, pIndex, pMult, scinHits);
			for( int hit = 0 ; hit < scinHits ; hit++ ){
				  hit_x[hit] = hitVector[hit].X();
				  hit_y[hit] = hitVector[hit].Y();
				  hit_z[hit] = hitVector[hit].Z();
			}

			// Grab the neutron information:
			// Grab the neutron information:
			getNeutronInfo( band_hits, band_rawhits, band_adc, band_tdc, nMult, nHit , starttime , runNum);
			TVector3 nPath[maxNeutrons];
			for( int n = 0 ; n < nMult ; n++ ){

				barID[n] 			= nHit[n].getBarID();
				nSector[n] 		= nHit[n].getSector();
				nLayer[n] 		= nHit[n].getLayer();
				nComponent[n] = nHit[n].getComponent();
				nEdep[n] 			= nHit[n].getEdep();
				nTof[n]				= nHit[n].getTof() -  TDC_BARSHIFTS[barID[n]];
				nTofFADC[n]		= nHit[n].getTofFadc() -  FADC_BARSHIFTS[barID[n]];
				n_bandx[n]		= nHit[n].getX();
				n_bandy[n] 		= nHit[n].getY();
				n_bandz[n] 		= nHit[n].getZ();
	      nPath[n].SetXYZ(n_bandx[n],n_bandy[n],n_bandz[n]);
				dL_n[n]				= nPath[n].Mag();
				theta_n[n]		= nPath[n].Theta();
				phi_n[n]			= nPath[n].Phi();
				nADCleft[n] 	= nHit[n].getPmtLadc();
				nADCright[n] 	= nHit[n].getPmtRadc();
				nStatus[n] 		= nHit[n].getStatus();
			}

			// Fill tree to do any more plots on the fly
			outTree->Fill();

		} // end loop over events
	//	cout << "Total charge collected in file: " << gated_charge << " [microC]\n";
		cout << "Total number of events written to output: " << outTree->GetEntries() << "\n";
	}// end loop over files


	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}


void getPositiveInfo( BParticle particles, double pid[maxPositive], TVector3 momentum[maxPositive], TVector3 vertex[maxPositive],	double time[maxPositive],
	double charge[maxPositive], double beta[maxPositive], double chi2pid[maxPositive], double status[maxPositive] , int index[maxPositive], int& multiplicity ){
	// Takes all positive particles in REC::Particles
	multiplicity = 0;
	for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
			if( particles.getCharge(row) == 1 ){
			pid[multiplicity] 		= particles.getPid(row);
			charge[multiplicity]		= particles.getCharge(row);
			momentum 	[multiplicity]	= particles.getV3P(row);
			vertex		[multiplicity]	= particles.getV3v(row);
			time		[multiplicity]	= particles.getVt(row);
			beta		[multiplicity]	= particles.getBeta(row);
			chi2pid		[multiplicity]	= particles.getChi2pid(row);
			status		[multiplicity]	= particles.getStatus(row);
			index [multiplicity] = row;
			multiplicity ++;
		}
	}
}

void getScinHits( BScintillator scintillator, double pindex[maxScinHits], double detid[maxScinHits], double energy[maxScinHits], double time[maxScinHits],
	 TVector3 posVector[maxScinHits], double path[maxScinHits], double status[maxScinHits], int posIndex[maxPositive], int posMult, int &scinHits) {
	scinHits = 0;
	for( int row = 0 ; row < scintillator.getRows() ; row++ ){
		for ( int i = 0 ; i < posMult ; i++) { //loop over all positive particles
      if (scintillator.getPindex(row) == posIndex[i]) { //check if hit Pindex corresponds to positive particle index = row from REC::Particles bank
				pindex[scinHits ] 		= scintillator.getPindex(row);
			  detid[scinHits ] 		= scintillator.getDetector(row);
				energy[scinHits ] 		= scintillator.getEnergy(row);
				time[scinHits ] 		= scintillator.getTime(row);
				path[scinHits ] 		= scintillator.getPath(row);
				posVector[scinHits].SetXYZ(	scintillator.getX(row), scintillator.getY(row), scintillator.getZ(row) 	);
				status		[scinHits]	= scintillator.getStatus(row);
				scinHits ++;
			}
		}
	}
}



bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot){
			// Anymore fiducial cuts that are needed for the electron:
			// 	E/p cut
			// 	vertex cut
			// 	minimum momentum  cut / maximum  momentum  cut
			// 	electron ToF cut
			// 	minimum W cut
			//	chi2 cut
			//	pcal energy cut?
			//

	if( pid != 11 || charge != -1 ) return false;
	//if( lV < 2 || lW < 2 ) return false;
	//if( momentum.Mag() < 1 || momentum.Mag() > 4.2 ) return false;
	//if( chi2pid == 0 || chi2pid > 1 ) return false;
	//if( E_tot/momentum.Mag() > 0.4 || E_tot/momentum.Mag() < 0.1 ) return false;
	//if( vertex.X() < -2 || vertex.X() > 2) return false;
	//if( vertex.Y() < -2 || vertex.Y() > 2) return false;
	//if( vertex.Z() < -7 || vertex.Z() > 2) return false;
	//if( time < 15 ) return false;
	//if( E_tot / momentum.Mag() < 0.15 || E_tot / momentum.Mag() > 0.3 ) return false;

	return true;
}
