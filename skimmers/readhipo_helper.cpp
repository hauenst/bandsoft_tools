#include "readhipo_helper.h"


void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun){
	
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
		hits[hit].setY			(band_hits.getY			(hit)			);
		hits[hit].setZ			(band_hits.getZ			(hit)			);
		hits[hit].setStatus		(band_hits.getStatus		(hit)			);


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

		// Save how many neutron hits we have
		mult++;
	}
	
}

int getRunNumber( string filename ){
	string parsed = filename.substr( filename.find("inc_") );
	string moreparse = parsed.substr(4,6);
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
void getElectronInfo( BParticle particles, BCalorimeter calorimeter, BScintillator scintillator,
			clashit &electron,
			double starttime , int thisRun , double Ebeam ){

	TVector3	momentum = particles.getV3P(0);
	TVector3	vertex	 = particles.getV3v(0);

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	qVec; qVec = beamVec - momentum;


	electron.setSector		(	calorimeter.getSector(0)				);
	electron.setPID			(	particles.getPid(0)					);
	electron.setCharge		(	particles.getCharge(0)					);
	electron.setStatus		(	particles.getStatus(0)					);

	electron.setTime		(	particles.getVt(0)					);
	electron.setBeta		(	particles.getBeta(0)					);
	electron.setChi2		(	particles.getChi2pid(0)					);
	electron.setEtot		(	calorimeter.getTotE(0)					);
	electron.setEpcal		(	calorimeter.getPcalE(0)					);
	electron.setEoP			(	electron.getEtot() / momentum.Mag()			);
	electron.setTimeScint		(	scintillator.getTime(0)-starttime			);
	electron.setPathScint		(	scintillator.getPath(0)					);
	electron.setU			(	calorimeter.getLU(0)					);
	electron.setV			(	calorimeter.getLV(0)					);
	electron.setW			(	calorimeter.getLW(0)					);
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

		tag[hit].setMomentumE	(eVec 		);
		tag[hit].setMomentumN	(nVec		);
		tag[hit].setMomentumQ	(qVec		);
		tag[hit].setMomentumB	(beamVec	);

		tag[hit].setPhiNQ	(phi_nq		);
		tag[hit].setThetaNQ	(theta_nq	);
		tag[hit].setWp		(Wp		);
		tag[hit].setXp		(Xp		);
		tag[hit].setAs		(As		);
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
