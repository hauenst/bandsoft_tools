#include "readhipo_helper.h"


void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun, std::map<int,double> &bar_x, std::map<int,double> &bar_y, std::map<int,double> &bar_z,
			std::map<int,double> &bar_edep, int hotfix, double* s6200_fadc_lroffset , double* s6200_tdc_lroffset,
			double* s6291_fadc_lroffset,	double* s6291_tdc_lroffset,
			double* s6200_fadc_effvel,	double* s6200_tdc_effvel,
	       		double* s6291_fadc_effvel,	double* s6291_tdc_effvel 	){

	if( band_hits.getRows() > maxNeutrons ) return; // not interested in events with more than max # BAND hits for now
	for( int hit = 0 ; hit < band_hits.getRows() ; hit++ ){
		//if( band_hits.getStatus(hit) != 0 ) continue;	// not interested in an event that has a veto hit

		// Set the hits information
		hits[hit].setSector		(band_hits.getSector		(hit)			);
		hits[hit].setLayer		(band_hits.getLayer		(hit)			);
		hits[hit].setComponent		(band_hits.getComponent		(hit)			);
		hits[hit].setBarID		(band_hits.getBarKey		(hit)			);
		hits[hit].setEdep		(band_hits.getEnergy		(hit)			);
		hits[hit].setTof		(band_hits.getTime		(hit) - starttime	);
		hits[hit].setTofFadc		(band_hits.getTimeFadc		(hit) - starttime	);
		hits[hit].setTdiff		(band_hits.getDifftimeTdc	(hit)			);
		hits[hit].setTdiffFadc		(band_hits.getDifftimeFadc	(hit)			);
		hits[hit].setX			(band_hits.getX			(hit)			);

		//Standard calculation of y and z if maps for bars are empty
		hits[hit].setY			(band_hits.getY			(hit)			);
		//correct for offset of BAND compared to ideal position
		hits[hit].setZ			(band_hits.getZ			(hit) + BAND_OFFSET	);
		//correct for target offset
		hits[hit].setZ			(band_hits.getZ			(hit) - VERTEX_OFFSET	);
		// Fix for the Y position for layer 5 in case not fixed by bar geometry map below:
		if( band_hits.getLayer(hit) == 5 && (band_hits.getSector(hit) == 3 || band_hits.getSector(hit) == 4 ) ){
			hits[hit].setY		(band_hits.getY			(hit) + 7.2		);
		}


		//If maps for bar geometry are not empty, use it for y and z
		if ( !bar_y.empty() ) {
				hits[hit].setY			( bar_y[band_hits.getBarKey		(hit)	]		);
		}
		if ( !bar_z.empty() ) {
		  	hits[hit].setZ			( bar_z[band_hits.getBarKey		(hit)	]	 + BAND_OFFSET	- VERTEX_OFFSET	);
		}


		// If layer == 6 and sector == 4, take the idxR, otherwise take idxL for the veto bars:
		//Also modify x position for veto if x map is not empty
		if( hits[hit].getLayer() == 6 ){
			if ( !bar_x.empty() ) {
			  	hits[hit].setX			( bar_x[band_hits.getBarKey		(hit)	]	);
			}
			int rawhit_idx = -1;
			if( hits[hit].getSector() == 4 ) rawhit_idx = band_hits.getRpmtindex(hit);
			else{ rawhit_idx = band_hits.getLpmtindex(hit); }

			hits[hit].setRawLtdc		(band_rawhits.getFloat( 7 , rawhit_idx ) 		);
			hits[hit].setRawLtdccorr	(band_rawhits.getFloat( 9 , rawhit_idx ) 		);
			hits[hit].setRawLtfadc		(band_rawhits.getFloat( 8 , rawhit_idx ) 		);
			hits[hit].setRawLamp		(band_rawhits.getFloat( 6 , rawhit_idx )		);
			hits[hit].setRawLadc		(band_rawhits.getFloat( 5 , rawhit_idx )		);

			hits[hit].setRawRtdc		(band_rawhits.getFloat( 7 , rawhit_idx ) 		);
			hits[hit].setRawRtdccorr	(band_rawhits.getFloat( 9 , rawhit_idx ) 		);
			hits[hit].setRawRtfadc		(band_rawhits.getFloat( 8 , rawhit_idx ) 		);
			hits[hit].setRawRamp		(band_rawhits.getFloat( 6 , rawhit_idx )		);
			hits[hit].setRawRadc		(band_rawhits.getFloat( 5 , rawhit_idx )		);

			//Use average Edep per channel for vetos (could be either RawLadc or RawRadc)
			hits[hit].setEdep		(hits[hit].getRawLadc()	);

			// Using the rawhit struct, get the raw PMT information to use later
			int pmtTdc	= band_rawhits.getInt( 10 , rawhit_idx );
			int pmtAdc	= band_rawhits.getInt( 11 , rawhit_idx );
			//	Get the raw pmt information corresponding to the band hit above
			hits[hit].setPmtLtdc		(band_tdc.getInt( 4 	, pmtTdc )		);
			hits[hit].setPmtRtdc		(band_tdc.getInt( 4 	, pmtTdc )		);
			hits[hit].setPmtLtfadc		(band_adc.getFloat( 6 	, pmtAdc )		);
			hits[hit].setPmtRtfadc		(band_adc.getFloat( 6 	, pmtAdc )		);
			hits[hit].setPmtLamp		(band_adc.getInt( 5 	, pmtAdc )		);
			hits[hit].setPmtRamp		(band_adc.getInt( 5 	, pmtAdc )		);
			hits[hit].setPmtLadc		(band_adc.getInt( 4 	, pmtAdc )		);
			hits[hit].setPmtRadc		(band_adc.getInt( 4 	, pmtAdc )		);
			hits[hit].setPmtLped		(band_adc.getInt( 7 	, pmtAdc )		);
			hits[hit].setPmtRped		(band_adc.getInt( 7 	, pmtAdc )		);


		}
		else{
			// Using the band hit struct, get the raw hit PMT information to use later
			int rawhit_idxL = band_hits.getLpmtindex(hit);
			int rawhit_idxR = band_hits.getRpmtindex(hit);
			// 	Get the raw hit information corresponding to the band hit above
			hits[hit].setRawLtdc		(band_rawhits.getFloat( 7 , rawhit_idxL ) 		);
			hits[hit].setRawRtdc		(band_rawhits.getFloat( 7 , rawhit_idxR ) 		);
			hits[hit].setRawLtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxL ) 		);
			hits[hit].setRawRtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxR ) 		);
			hits[hit].setRawLtfadc		(band_rawhits.getFloat( 8 , rawhit_idxL ) 		);
			hits[hit].setRawRtfadc		(band_rawhits.getFloat( 8 , rawhit_idxR ) 		);
			hits[hit].setRawLamp		(band_rawhits.getFloat( 6 , rawhit_idxL )		);
			hits[hit].setRawRamp		(band_rawhits.getFloat( 6 , rawhit_idxR )		);
			hits[hit].setRawLadc		(band_rawhits.getFloat( 5 , rawhit_idxL )		);
			hits[hit].setRawRadc		(band_rawhits.getFloat( 5 , rawhit_idxR )		);

			// Using the rawhit struct, get the raw PMT information to use later
			int pmtTdcL	= band_rawhits.getInt( 10 , rawhit_idxL );
			int pmtAdcL	= band_rawhits.getInt( 11 , rawhit_idxL );
			int pmtTdcR	= band_rawhits.getInt( 10 , rawhit_idxR );
			int pmtAdcR	= band_rawhits.getInt( 11 , rawhit_idxR );
			//	Get the raw pmt information corresponding to the band hit above
			hits[hit].setPmtLtdc		(band_tdc.getInt( 4 , pmtTdcL )		);
			hits[hit].setPmtRtdc		(band_tdc.getInt( 4 , pmtTdcR )		);
			hits[hit].setPmtLtfadc		(band_adc.getFloat( 6 , pmtAdcL )	);
			hits[hit].setPmtRtfadc		(band_adc.getFloat( 6 , pmtAdcR )	);
			hits[hit].setPmtLamp		(band_adc.getInt( 5 , pmtAdcL )		);
			hits[hit].setPmtRamp		(band_adc.getInt( 5 , pmtAdcR )		);
			hits[hit].setPmtLadc		(band_adc.getInt( 4 , pmtAdcL )		);
			hits[hit].setPmtRadc		(band_adc.getInt( 4 , pmtAdcR )		);
			hits[hit].setPmtLped		(band_adc.getInt( 7 , pmtAdcL )		);
			hits[hit].setPmtRped		(band_adc.getInt( 7 , pmtAdcR )		);
		}

		// calculate x- to do a hot fix for position using fadc time due to TDC shift in our 10.2 runs, dont do it for vetos
		int id 			= band_hits.getBarKey(hit);
		if( hotfix == 1 && hits[hit].getLayer() != 6){
			double old_tdiff_fadc 	= band_hits.getDifftimeFadc(hit);
			double old_tdiff_tdc	= band_hits.getDifftimeTdc(hit);


			// Shift the TDIFF_TDC by the residual correction for 10.2 GeV data
			// and recalculate x-position to get the global offset corrections.
			// Then we can use NEW effective velocity tables to calculate a correct x-position,
			// only for data and only for 10.2 GeV data
			double old_x_fadc 		= (-1./2) * old_tdiff_fadc * s6200_fadc_effvel[id];
			double old_x_tdc		= (-1./2) * old_tdiff_tdc  * s6200_tdc_effvel[id];

			// recoX = x_tdc + globX
			// reco_offset = recoX - x_tdc = globX
			double fadc_reco_offset	= band_hits.getX(hit) - old_x_fadc;
			double tdc_reco_offset	= band_hits.getX(hit) - old_x_tdc;

			// Now calculate the real x-position with better FADC/TDC effective velocity
			// newX = x_tdc + reco_offset
			double new_tdiff_fadc 	= band_hits.getDifftimeFadc(hit) + s6200_fadc_lroffset[id] - s6291_fadc_lroffset[id];
			double new_tdiff_tdc	= band_hits.getDifftimeTdc(hit)	 + s6200_tdc_lroffset[id]  - s6291_tdc_lroffset[id];
			double new_x_fadc 	= (-1./2) * new_tdiff_fadc * s6291_fadc_effvel[id];
			double new_x_tdc	= (-1./2) * new_tdiff_tdc  * s6291_tdc_effvel[id];

			new_x_fadc 	= new_x_fadc 	+ fadc_reco_offset;
			new_x_tdc 	= new_x_tdc 	+ tdc_reco_offset;
			hits[hit].setXFadc	( new_x_fadc );
			hits[hit].setX		( new_x_tdc );

			//cout << "old tdc diff:\t" << old_tdiff_tdc << "\n";
			//cout << "old tdc xpos:\t" << old_x_tdc << "\n";
			//cout << "tdc reco off:\t" << tdc_reco_offset << "\n";
			//cout << "new tdc diff:\t" << new_tdiff_tdc << "\n";
			//cout << "new tdc xpos:\t" << new_x_tdc << "\n";
		}

		//If Edep maps are not empty, use them to calibrate Edep in MeV for each bar
		//For MC Edep other file is used with constants
		if ( !bar_edep.empty() ) {
			//Calibration values are in channel/MeV
			//cout << "applying energycorrection before value: " << hits[hit].getEdep() << endl;
				hits[hit].setEdep		(band_hits.getEnergy	(hit)	/ bar_edep[band_hits.getBarKey	(hit)	]		);
			//		cout << "applying energycorrection after value: " << hits[hit].getEdep() << endl;
		}

		// After all corrections to positions, set the pathlength and status:
		hits[hit].setStatus		(band_hits.getStatus		(hit)			);
		hits[hit].setDL			( TVector3(hits[hit].getX(),hits[hit].getY(),hits[hit].getZ()) );
		// Save how many neutron hits we have
		mult++;
	}

}

bool goodNeutronEvent(bandhit hits[maxNeutrons], int nMult, int& leadindex, int mcdataselect){

	double time_thres = time_thresBANDhit; //300

	bool vetoHit = false;
	bool accept = false;

	// Loop over all hits in BAND and get the earliest hit in time
	// but filter out hits that have low Edep
	double fastestTime = 1E5;
	int fastestTimeIdx = -1;
	for( int thishit = 0; thishit < nMult ; thishit++){
		bandhit neutron = hits[thishit];
		//if the Edep in this PMT is less than 2MeVee, do not count it
		//Edep is assumed to be calibrated to MeV for each bar
		//MC
		if (mcdataselect == 0 && neutron.getEdep() < 2) 	continue;
		//Data
		if( mcdataselect != 0 && neutron.getEdep() < 2 ) continue;

		if( neutron.getTof() < fastestTime ){
			fastestTime = neutron.getTof();
			fastestTimeIdx = thishit;
		}
	}

	// Once we have the fastest hit, ask how many other hits there
	// are in the same event and flag how far in time they are,
	// but keep this earliest hit.
	// Accept the event based on if there is another hit > 300ps
	if( fastestTimeIdx != -1 ){
		accept = true;
		for( int thishit = 0; thishit < nMult ; thishit++){
			bandhit neutron = hits[thishit];
			// if we have a veto bar hit and if the Edep is > 0.25MeVee,  need to use getPmtLADC because no Edep
			//Edep is assumed to be calibrated to MeV for each bar
			if( neutron.getLayer() == 6 && neutron.getEdep() > 0.25 ) vetoHit = true;

			//if this bar does not have > 2MeVee, do not count the hit
			//MC
			if (mcdataselect == 0 && neutron.getEdep() < 2 ) continue;
			//Data
			if( mcdataselect != 0 && neutron.getEdep() < 2 ) continue;

			double tdiff = neutron.getTof() - fastestTime;
			if( tdiff == 0 ){
				continue;
			}
			bool adjLayer = false;
			bool adjComp = false;
			bool adjSpace = false;
			bandhit lead = hits[thishit];

			int layerDiff = lead.getLayer() - neutron.getLayer();
			if( fabs(layerDiff) <= 1 ) adjLayer = true;
			if( lead.getSector() == neutron.getSector() ){
			// if hits are in the same sector, then just ask for comp diff
				int compDiff = lead.getComponent() - neutron.getComponent();
				if( fabs(compDiff) <= 1 ) adjComp = true;
			}
			else{
			// have to implement some fancy component check due to sec comp differences
				int yDiff = lead.getY() - neutron.getY();
				if( fabs(yDiff) <= 10 ) adjComp = true;
			}
			if( adjComp && adjLayer ) adjSpace = true;

			if( tdiff < time_thres && tdiff != 0 ) accept = false;
			if( tdiff < time_thres && tdiff != 0 && adjSpace ) accept = true;
		}
	}

	leadindex = fastestTimeIdx ;

	return accept && !vetoHit;

}


int getRunNumber( string filename ){
        //string parsed = filename.substr( filename.find("inc_") );
        string parsed;
        string moreparse;
        if ( filename.find("inc_") <= filename.length() )
        {
                cout << "Parsed file and found position for string inc_ at " << filename.find("inc_") << endl;
                parsed = filename.substr( filename.find("inc_") );
                moreparse = parsed.substr(4,6);
        }
        else if (filename.find("band_") <= filename.length() )
        {
                cout << "Parsed file and found position for string band_ at " << filename.find("band_") << endl;
                parsed = filename.substr( filename.find("band_") );
                moreparse = parsed.substr(5,6);
        }
        else if (filename.find("rec_clas_") <= filename.length() )
        {
                cout << "Parsed file and found position for string rec_clas_ at " << filename.find("rec_clas_") << endl;
                parsed = filename.substr( filename.find("rec_clas_") );
                moreparse = parsed.substr(9,6);
        }
				else if (filename.find("gmn_") <= filename.length() )
        {
                cout << "Parsed file and found position for string gmn_ at " << filename.find("gmn_") << endl;
                parsed = filename.substr( filename.find("gmn_") );
                moreparse = parsed.substr(4,6);
        }
        else {
                cout << "Could not parse runnumber from inputfile. Return 0 runnumber " << endl;
                return 0;
        }
        cout << "\t*Intepreted run number from file name: " << stoi(moreparse) << "\n";
        return stoi(moreparse);
}


void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime ){
	if( eventInfo.getRows() != 1 ){
		cerr << "getEventInfo::NotImplementedFunction\n";
		exit(-1);
	}
	integrated_charge       = (double)eventInfo.getBCG(0); 	// not calibrated currently
	livetime 		= (double)eventInfo.getLT(0);		// not calibrated currently
	starttime		= (double)eventInfo.getSTTime(0);
	return;
}

void getMcInfo( hipo::bank gen_particles , hipo::bank gen_info , genpart mcParts[maxGens] ,
		double &starttime, double &weight, double &Ebeam , int &genMult ){
	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec;
	bool setElectron = false;

	// Grab the weight for the event:
	weight 	= gen_info.getFloat(9,0);
	// Grab the beam energy for this generated file:
	Ebeam 	= gen_info.getFloat(6,0);
	//ugly temporary fix for the MC beam energy for the QE simulations
	if (Ebeam > 4 && Ebeam < 4.5) {
		Ebeam = 4.247;
	}

	// Loop over all generated particles and find the electron to put that one first
	for( int hit = 0 ; hit < gen_particles.getRows() ; hit++ ){
		int pid	= gen_particles.getInt( 0 , hit );
		if( pid == 11 && setElectron == false ){
			double px = gen_particles.getFloat( 1 , hit );
			double py = gen_particles.getFloat( 2 , hit );
			double pz = gen_particles.getFloat( 3 , hit );
			double vx = gen_particles.getFloat( 4 , hit );
			double vy = gen_particles.getFloat( 5 , hit );
			double vz = gen_particles.getFloat( 6 , hit );
			double vt = gen_particles.getFloat( 7 , hit );
			eVec.SetXYZ( px, py, pz );
			TVector3	vertex; vertex.SetXYZ( vx, vy, vz );
			TVector3	qVec; qVec = beamVec - eVec;
			starttime = vt;

			mcParts[genMult].setPID		( pid 							);
			mcParts[genMult].setMomentum	( eVec.Mag()					);
			mcParts[genMult].setTheta	( eVec.Theta()					);
			mcParts[genMult].setPhi		( eVec.Phi()					);

			mcParts[genMult].setQ		( qVec.Mag()						);
			mcParts[genMult].setThetaQ	( qVec.Theta()						);
			mcParts[genMult].setPhiQ	( qVec.Phi()						);

			mcParts[genMult].setOmega	( Ebeam - sqrt( pow(eVec.Mag(),2) + mE*mE )		);
			mcParts[genMult].setQ2		( qVec.Mag()*qVec.Mag() - pow(mcParts[genMult].getOmega(),2)	);
			mcParts[genMult].setXb		( mcParts[genMult].getQ2()/(2.*mP*mcParts[genMult].getOmega())		);
			mcParts[genMult].setW2		( mP*mP - mcParts[genMult].getQ2() + 2.*mcParts[genMult].getOmega()*mP	);

			mcParts[genMult].setVtx		(vertex.X()						);
			mcParts[genMult].setVty		(vertex.Y()						);
			mcParts[genMult].setVtz		(vertex.Z()						);

			setElectron = true;
			genMult++;
		}
		else if( pid != 11 && setElectron == true ){
			double px = gen_particles.getFloat( 1 , hit );
			double py = gen_particles.getFloat( 2 , hit );
			double pz = gen_particles.getFloat( 3 , hit );
			double vx = gen_particles.getFloat( 4 , hit );
			double vy = gen_particles.getFloat( 5 , hit );
			double vz = gen_particles.getFloat( 6 , hit );
			double vt = gen_particles.getFloat( 7 , hit );

			TVector3	momentum; 	momentum.SetXYZ( px, py, pz );
			TVector3	vertex; 	vertex.SetXYZ( vx, vy, vz );

			mcParts[genMult].setPID		( pid 							);
			mcParts[genMult].setMomentum	( momentum.Mag()					);
			mcParts[genMult].setTheta	( momentum.Theta()					);
			mcParts[genMult].setPhi		( momentum.Phi()					);

			mcParts[genMult].setVtx		(vertex.X()						);
			mcParts[genMult].setVty		(vertex.Y()						);
			mcParts[genMult].setVtz		(vertex.Z()						);

			genMult++;
		}
	}



	return;
}


void getElectronInfo( BParticle particles, BCalorimeter calorimeter, BScintillator scintillator, hipo::bank DC_Track, hipo::bank DC_Traj,
			int pbankIndex,
			clashit &electron,
			double starttime , int thisRun , double Ebeam ){

	TVector3	momentum = particles.getV3P(pbankIndex);
	TVector3	vertex	 = particles.getV3v(pbankIndex);

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	qVec; qVec = beamVec - momentum;

//for Calorimeter information it is the Particle bank index for the first particle (usually electron)
//for Particle bank it is bank index 0
	electron.setSector		(	calorimeter.getElectronSector(pbankIndex)				);
	electron.setPID			  (	particles.getPid(pbankIndex)					);
	electron.setCharge		(	particles.getCharge(pbankIndex)					);
	electron.setStatus		(	particles.getStatus(pbankIndex)					);

	electron.setTime		(	particles.getVt(pbankIndex)					);
	electron.setBeta		(	particles.getBeta(pbankIndex)					);
	electron.setChi2		(	particles.getChi2pid(pbankIndex)					);
	electron.setEtot		(	calorimeter.getTotE(pbankIndex)					);
	electron.setEpcal		(	calorimeter.getPcalE(pbankIndex)					);
	electron.setEecin		(	calorimeter.getECinE(pbankIndex)					);
	electron.setEecout		(	calorimeter.getECoutE(pbankIndex)   		);
	electron.setEoP			(	electron.getEtot() / momentum.Mag()			);
	if (calorimeter.getPcalRow(pbankIndex) > -1) {
		electron.setU			(	calorimeter.getLU(calorimeter.getPcalRow(pbankIndex))					);
		electron.setV			(	calorimeter.getLV(calorimeter.getPcalRow(pbankIndex))					);
		electron.setW			(	calorimeter.getLW(calorimeter.getPcalRow(pbankIndex))					);
	}
	electron.setVtx			(	vertex.X()						);
	electron.setVty			(	vertex.Y()						);
	electron.setVtz			(	vertex.Z()						);

	electron.setMomentum		(	momentum.Mag()						);
	electron.setTheta		(	momentum.Theta()					);
	electron.setPhi			(	momentum.Phi()						);

	electron.setQ			(	qVec.Mag()						);
	electron.setThetaQ		(	qVec.Theta()						);
	electron.setPhiQ		(	qVec.Phi()						);

	electron.setOmega		(	Ebeam - sqrt( pow(momentum.Mag(),2) + mE*mE )		);
	electron.setQ2			(	qVec.Mag()*qVec.Mag() - pow(electron.getOmega(),2)	);
	electron.setXb			(	electron.getQ2()/(2.*mP*electron.getOmega())		);
	electron.setW2			(	mP*mP - electron.getQ2() + 2.*electron.getOmega()*mP	);


  //Adding Tracking bank

	int Ntrack = DC_Track.getRows();
	int index = -1;
	int count =0;
	for(int i =0; i < Ntrack; i++){
	  int pindex   = DC_Track.getInt(1,i);
	  int detector = DC_Track.getInt(2,i);

	  if (pindex == pbankIndex && detector ==6 )
	    {index = i;
	     count ++;
	   }
	}

	//At this point I only use 1 track event
	//The check show around 4% mul_track events

	if (index != -1 && count ==1) {
	  electron.setDC_chi2             (       DC_Track.getFloat(6, index)                                 );
	  electron.setDC_NDF              (       DC_Track.getInt(7, index)                                   );
	  electron.setDC_sector           (       DC_Track.getInt(3, index)                                   );

	}
	else {
		electron.setDC_sector           (       electron.getSector()                                );
	}
	if (count>1) {
		cout << "getElectronInfo Warning: DC count is larger than 1. Count =  " << count << endl;
	}
	//just give out warning for electron in first particle bank row
	if (index == -1 && pbankIndex == 0) {
		//cout << "getElectronInfo Warning: DC index is -1" << endl;
	}



	//Adding Traj bank, Structure is more complicated
	//Also only consider the traj of 1-track event

	int NTraj = DC_Traj.getRows();
	int index_traj = -1;

	for(int i =0; i <NTraj; i++){
    int pindex   = DC_Traj.getInt(0, i);
	  int detector = DC_Traj.getInt(2,i);
	  int layer    = DC_Traj.getInt(3,i);

	  if(pindex == pbankIndex && detector ==6 && layer ==6 && count ==1)

	     {
	       electron.setDC_x1        (  DC_Traj.getFloat(4, i)                                           );
	       electron.setDC_y1        (  DC_Traj.getFloat(5, i)                                           );
	       electron.setDC_z1        (  DC_Traj.getFloat(6, i)                                           );

	     }

	  if(pindex == pbankIndex && detector == 6 && layer ==18 && count ==1)

	    {
	       electron.setDC_x2        (  DC_Traj.getFloat(4, i)                                           );
	       electron.setDC_y2        (  DC_Traj.getFloat(5, i)                                           );
	       electron.setDC_z2        (  DC_Traj.getFloat(6, i)                                           );


	    }

	  if(pindex == pbankIndex && detector == 6 && layer ==36 && count ==1)

	    {
	       electron.setDC_x3        (  DC_Traj.getFloat(4, i)                                           );
	       electron.setDC_y3        (  DC_Traj.getFloat(5, i)                                           );
	       electron.setDC_z3        (  DC_Traj.getFloat(6, i)                                           );

	    }


	}


	//if( pid != 11 || charge != -1 ) return false;
	//if( lV < 2 || lW < 2 ) return false;
	//if( momentum.Mag() < 1 || momentum.Mag() > 4.2 ) return false;
	//if( chi2pid == 0 || chi2pid > 1 ) return false;
	//if( E_tot/momentum.Mag() > 0.4 || E_tot/momentum.Mag() < 0.1 ) return false;
	//if( vertex.X() < -2 || vertex.X() > 2) return false;
	//if( vertex.Y() < -2 || vertex.Y() > 2) return false;
	//if( vertex.Z() < -7 || vertex.Z() > 2) return false;
	//if( time < 15 ) return false;
	//if( momentum.Mag() < 2 || momentum.Mag() > 10.6 ) return false;
	//if( E_tot / momentum.Mag() < 0.15 || E_tot / momentum.Mag() > 0.3 ) return false;

	return;
}

void getTaggedInfo( clashit eHit, bandhit nHit[maxNeutrons], taghit tag[maxNeutrons] ,
		double Ebeam , int nMult ){

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; eVec.SetMagThetaPhi( eHit.getMomentum(), eHit.getTheta(), eHit.getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;

	// Loop over all neutrons to combine with the electron
	for( int hit = 0 ; hit < nMult ; hit++ ){

		TVector3	nVec;
		nVec = (nHit[hit].getDL()).Unit();


		TVector3 norm_scatter = qVec.Cross( beamVec );
		norm_scatter 	= norm_scatter.Unit();

		TVector3 norm_reaction = qVec.Cross( nVec );
		norm_reaction 	= norm_reaction.Unit();

		double phi_nq 		= norm_scatter.Angle( norm_reaction );
		double theta_nq 	= nVec.Angle( qVec );
		double CosTheta_nq 	= cos(theta_nq);

		TVector3 direction = norm_scatter.Cross(norm_reaction);
		if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
		}
		else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
			phi_nq *= (-1);
		}

		double beta = nHit[hit].getDL().Mag() / (nHit[hit].getTofFadc()*cAir);
		double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
		nVec.SetMagThetaPhi(p_n,nHit[hit].getDL().Theta(),nHit[hit].getDL().Phi());

		double E_n 	= sqrt( mN*mN + p_n*p_n );
		double W_primeSq = mD*mD - eHit.getQ2() + mN*mN + 2.*mD*(eHit.getOmega()-E_n) - 2.*eHit.getOmega()*E_n + 2.*eHit.getQ()*p_n*cos(theta_nq);
		double Wp = sqrt(W_primeSq);
		double Xp = eHit.getQ2()/(2.*( eHit.getOmega()*(mD-E_n) + p_n*eHit.getQ()*CosTheta_nq));
		double As = (E_n - p_n*CosTheta_nq)/mN;
		double Xp2 = eHit.getQ2()/(W_primeSq - mN*mN + eHit.getQ2());

		TVector3 Pt;
		TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;
		Pt = nVec - pN_par_q;

		tag[hit].setMomentumE	(eVec 		);
		tag[hit].setMomentumN	(nVec		);
		tag[hit].setMomentumQ	(qVec		);
		tag[hit].setMomentumB	(beamVec	);

		tag[hit].setPhiNQ	(phi_nq		);
		tag[hit].setThetaNQ	(theta_nq	);
		tag[hit].setWp		(Wp		);
		tag[hit].setXp		(Xp		);
		tag[hit].setAs		(As		);
		tag[hit].setPt		(Pt		);
		tag[hit].setXp2		(Xp2		);
	}

	//beta = dL / (ToF*cAir);
	//p_n = mN / sqrt( 1./pow(beta,2) - 1. );
	//nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
	//E_n 	= sqrt( mN*mN + p_n*p_n );
	//double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
	//Wp = sqrt(W_primeSq);
	//Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
	//As = (E_n - p_n*CosTheta_nq)/mN;


	return;
}

void getScinHits( BScintillator scintillator, int pindex[maxScinHits], int detid[maxScinHits], double energy[maxScinHits], double time[maxScinHits],
	 TVector3 posVector[maxScinHits], double path[maxScinHits], int status[maxScinHits], int posIndex[maxParticles], int posMult, int &scinHits) {
	scinHits = 0;
	for( int row = 0 ; row < scintillator.getRows() ; row++ ){
		for ( int i = 0 ; i < posMult ; i++) { //loop over all other particles
      if (scintillator.getPindex(row) == posIndex[i]) { //check if hit Pindex corresponds to particle index = row from REC::Particles bank
				pindex[scinHits ] 		= scintillator.getPindex(row); //int
			  detid[scinHits ] 		= scintillator.getDetector(row); //int
				energy[scinHits ] 		= scintillator.getEnergy(row); //float
				time[scinHits ] 		= scintillator.getTime(row);//float
				path[scinHits ] 		= scintillator.getPath(row);//float
				posVector[scinHits].SetXYZ(	scintillator.getX(row), scintillator.getY(row), scintillator.getZ(row) 	);
				status		[scinHits]	= scintillator.getStatus(row); //int
				scinHits ++;
			}
		}
	}
}

void getParticleInfo( BParticle particles, int pid[maxParticles], TVector3 momentum[maxParticles], TVector3 vertex[maxParticles],	double time[maxParticles],
	int charge[maxParticles], double beta[maxParticles], double chi2pid[maxParticles], int status[maxParticles] , int index[maxParticles], int& multiplicity ){
	// Takes all particles in REC::Particles other than leading particle
	multiplicity = 0;
	for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
			if( particles.getCharge(row) == 1 || particles.getCharge(row) == -1 ){ //checks for positive and negative charge
			pid[multiplicity] 		= particles.getPid(row); //int
			charge[multiplicity]		= particles.getCharge(row); //int
			momentum 	[multiplicity]	= particles.getV3P(row); //TVector3
			vertex		[multiplicity]	= particles.getV3v(row);//TVector3
			time		[multiplicity]	= particles.getVt(row);//float
			beta		[multiplicity]	= particles.getBeta(row);//float
			chi2pid		[multiplicity]	= particles.getChi2pid(row);//float
			status		[multiplicity]	= particles.getStatus(row); //int
			index [multiplicity] = row;
			multiplicity ++;
		}
	}
}

void getBANDBarGeometry(string filename, std::map<int,double> &bar_x, std::map<int,double> &bar_y, std::map<int,double> &bar_z) {
	ifstream f;
	int sector, layer, component, barId;
	double x, y, z;
	f.open(filename);
	if (!f.is_open()) {
		cout << "ERROR BAND Bar geometry file could not be read " << endl;
		exit(-1);
	}
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> x; //units are in cm
		f >> y; //units are in cm
		f >> z; //units are in cm
		bar_x.insert(std::pair<int,double>(barId,x));
		bar_y.insert(std::pair<int,double>(barId,y));
		bar_z.insert(std::pair<int,double>(barId,z));
	}
	//140 is number of bars. Just check y (x and z should be similar)
	if (bar_y.size() != 140) {
		cout << "ERROR BAND Geometry y map length is not 140. Length is  " << bar_y.size() << endl;
		exit(-1);
	}
	f.close();
}

void getBANDEdepCalibration(string filename, std::map<int,double> &bar_edep) {
	ifstream f;
	int sector, layer, component, barId;
	double edep;
	f.open(filename);
	if (!f.is_open()) {
		cout << "ERROR BAND Edep file could not be read " << endl;
		exit(-1);
	}
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> edep; //units are in cm
		bar_edep.insert(std::pair<int,double>(barId,edep));
	//	cout << "BAND Edep " << sector << " " << layer << " " << component << " " << edep << endl;
	}
	f.close();
	//140 is number of bars
	if (bar_edep.size() != 140) {
		cout << "ERROR BAND Edep map length is not 140. Length is  " << bar_edep.size() << endl;
		exit(-1);
	}
}

void shiftsReader::LoadInitBar( string filename ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitBar[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitBarFadc( string filename ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitBarFadc[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitRun( string filename ){
	ifstream f;
	int runno;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> runno;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitRun[runno] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitRunFadc( string filename ){
	ifstream f;
	int runno;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> runno;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitRunFadc[runno] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadEffVel( string filename_S6200 , string filename_S6291 ){
	ifstream f;
	int sector, layer, component, barId;
	double tdc_ev, fadc_ev, temp;

	// Load Spring 2019 6200-6290 constants
	f.open(filename_S6200);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		f >> tdc_ev;
		f >> fadc_ev;
		barId = sector*100 + layer*10 + component;
		S6200FadcEffVel[barId] 	= fadc_ev;
		S6200TdcEffVel[barId] 	= tdc_ev;
		f >> temp;
		f >> temp;
	}
	f.close();

	// Load Spring 2019 6291-INF constants
	f.open(filename_S6291);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		f >> tdc_ev;
		f >> fadc_ev;
		barId = sector*100 + layer*10 + component;
		S6291FadcEffVel[barId] 	= fadc_ev;
		S6291TdcEffVel[barId] 	= tdc_ev;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadLrOff( string filename_S6200 , string filename_S6291 ){
	ifstream f;
	int sector, layer, component, barId;
	double tdc_lr, fadc_lr, temp;

	// Load Spring 2019 6200-6290 constants
	f.open(filename_S6200);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		f >> tdc_lr;
		f >> fadc_lr;
		barId = sector*100 + layer*10 + component;
		S6200FadcLrOffsets[barId] 	= fadc_lr;
		S6200TdcLrOffsets[barId] 	= tdc_lr;
		f >> temp;
		f >> temp;
	}
	f.close();

	// Load Spring 2019 6291-INF constants
	f.open(filename_S6291);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		f >> tdc_lr;
		f >> fadc_lr;
		barId = sector*100 + layer*10 + component;
		S6291FadcLrOffsets[barId] 	= fadc_lr;
		S6291TdcLrOffsets[barId] 	= tdc_lr;
		f >> temp;
		f >> temp;
	}
	f.close();
}

double * shiftsReader::getInitBar(void){
	return InitBar;
}
double * shiftsReader::getInitBarFadc(void){
	return InitBarFadc;
}
double * shiftsReader::getInitRun(void){
	return InitRun;
}
double * shiftsReader::getInitRunFadc(void){
	return InitRunFadc;
}
double * shiftsReader::getFadcEffVel(int Runno){
	if( Runno < 6291 ){
		return S6200FadcEffVel;
	}
	else{
		return S6291FadcEffVel;
	}
}
double * shiftsReader::getTdcEffVel(int Runno){
	if( Runno < 6291 ){
		return S6200TdcEffVel;
	}
	else{
		return S6291TdcEffVel;
	}
}
double * shiftsReader::getFadcLrOff(int Runno){
	if( Runno < 6291 ){
		return S6200FadcLrOffsets;
	}
	else{
		return S6291FadcLrOffsets;
	}
}
double * shiftsReader::getTdcLrOff(int Runno){
	if( Runno < 6291 ){
		return S6200TdcLrOffsets;
	}
	else{
		return S6291TdcLrOffsets;
	}
}
