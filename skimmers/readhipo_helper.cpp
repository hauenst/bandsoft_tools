#include "readhipo_helper.h"


bool goodNeutronEvent(bandhit hits[maxNeutrons], int nMult, int& leadindex, int mcdataselect, int& nPass ){

	if( nMult > maxNeutrons ){
		cerr << "readhipo_helper::Multiplicity larger than storage container.\n\t...exiting...\n";
		exit(-1);
	}

	//cout << "ENTERING GOOD NEUTRON\n";	
	int unblocked_hits = 0;
	bool passed[maxNeutrons];
	// First let's skim the hit list to those hits that are NOT blocked:
	for( int i1 = 0; i1 < nMult; ++i1){
		bandhit checkMe = hits[i1];

		if( checkMe.getLayer() == 6 ){
			passed[i1]=false;
			continue; 
		}// this hit cannot be blocked by anyone since it's in layer 6 but we don't want to count it as a good hit

		passed[i1] = true;
		//cout << "i1 information: \n";
		//checkMe.Print();
		for( int i2 = 0; i2 < nMult; ++i2 ){
			if( i2 == i1 ) continue;
			bandhit checkAgainst = hits[i2];
			//cout << "i2 information: \n";
			//checkAgainst.Print();
			//cout << "\n";

			if( checkAgainst.getLayer() - checkMe.getLayer() != 1 ) continue; 
			// this means i2 is not DIRECTLY in front of i1 

			// if checkAgainst is a VETO, then we have some special rules:
			if( checkAgainst.getLayer() == 6 ){
				if( fabs( checkMe.getDL().Y() - checkAgainst.getDL().Y() ) > 8 ) continue; 
				// this means i2 is DIRECTLY in front of i1 -- BUT it is NOT +/- 1 bar of i1 (i.e. i2 is not blocking i1)

				if( fabs( checkAgainst.getTof() - checkMe.getTof() ) > 15 ) continue; 	
				// this means i2 does not have a similar time as i1, so it does not block i1.
			}

			else{
				if( fabs( checkMe.getDL().Y() - checkAgainst.getDL().Y() ) > 8 ) continue; 	
				// this means i2 is DIRECTLY in front of i1 -- BUT it is NOT +/- 1 bar of i1 (i.e. i2 is not blocking i1)

				if( fabs( checkAgainst.getDL().X() - checkMe.getDL().X() ) > 15 ) continue;  
				// this means i2 is DIRECTLY in front of i1 AND it is +/- 1 bar of i2, but it is NOT a similar x (i.e. i2 is not blocking i1)

				if( fabs( checkAgainst.getTof() - checkMe.getTof() ) > 3 ) continue; 		
				// this means i2 does not have a similar time as i1, so it does not block i1.
			}

			// otherwise, i2 is blocking i1
			passed[i1] = false;
		}
		//cout << "is i1 blocked?: " << passed[i1] << "\n";
		if( passed[i1] == true ) ++unblocked_hits;
	}

	nPass = unblocked_hits;
	leadindex = -1;
	//cout << unblocked_hits << " ";
	// Now we have a list of hits that have PASSED and are NOT blocked.
	if( unblocked_hits == 0 || unblocked_hits > 2 ){
		// these events are ambiguous so let's not use them
		//cout << "don't take\n";
		return false;
	}
	else if( unblocked_hits == 2 ){
		for( int i1 = 0; i1 < nMult; ++i1){
			if( passed[i1] != 1 ) continue;
			bandhit this_hit = hits[i1];

			for( int i2 = 0; i2 < nMult; ++i2 ){
				if( i2 == i1 ) continue;
				if( passed[i2] != 1 ) continue;
				bandhit other_hit = hits[i2];

				// if checkAgainst is in the SAME layer as checkMe, then let's see if all the blocking conditions are
				// true, and if they are, let's CLUSTER them and take the one that is earliest in time
				if( this_hit.getLayer() - other_hit.getLayer() == 0 ){

					if( fabs( this_hit.getDL().Y()  - other_hit.getDL().Y() ) > 8 	) 	continue; 	
					if( fabs( this_hit.getDL().X()  - other_hit.getDL().X() ) > 15 ) 	continue;  
					if( fabs( this_hit.getTof() 	- other_hit.getTof() ) > 3 	) 	continue;

					// These two hits can be clustered together so let's just take
					// one of them:

					if( this_hit.getTof() < other_hit.getTof() ){
						leadindex = i1;
						//cout << "clustered two\n";
						return true;
					}
					else{
						leadindex = i2;
						//cout << "clustered two\n";
						return true;
					}

					continue;

				} // end if on layer
			} // end loop over i2
		} // end loop over i1

		// Could not cluster so just skip ambiguous event
		//cout << "don't take\n";
		return false;
	}
	else if( unblocked_hits == 1 ){
		for( int i1 = 0; i1 < nMult; ++i1){
			bandhit this_hit = hits[i1];
			if( passed[i1] != 1 ) continue;
			//cout << "just one\n";
			leadindex = i1;
			return true;
		}
	}
	else{ std::cerr << "unexpected nalgorithm issue\n\t...exiting...\n"; exit(-1); }

	//cout <<"wtf\n";
	return false;


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


void getEventInfo( hipo::bank eventInfo, double &integrated_charge, double &livetime, double &starttime ){
	if( eventInfo.getRows() != 1 ){
		cerr << "getEventInfo::NotImplementedFunction\n";
		exit(-1);
	}
	livetime 		= eventInfo.getDouble(3,0);
	integrated_charge 	= eventInfo.getFloat(2,0);
	starttime 		= eventInfo.getFloat(4,0);

	return;
}

void getMcInfo( hipo::bank gen_particles , hipo::bank gen_info , genpart mcParts[maxGens] ,
		double &starttime, double &weight, double &Ebeam , int &genMult ){

	// Grab the weight for the event:
	weight 	= gen_info.getFloat(9,0);
	// Grab the beam energy for this generated file:
	double file_Ebeam = gen_info.getFloat(6,0);
	if( file_Ebeam == 0 ) return; // don't do anything if there is no event information
	if( fabs(Ebeam - file_Ebeam)>0.001 && file_Ebeam != -99 ) 
		std::cerr << "---WARNING-- getMcInfo shows different beam energy than the expected period energy!\n";
	if( file_Ebeam > 0 ){
		Ebeam = file_Ebeam;
	}

	// using the header Ebeam, create the beamvector:
	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec;
	bool setElectron = false;

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


void getElectronInfo( BParticle particles, BCalorimeter calorimeter, hipo::bank scintillator, hipo::bank DC_Track, hipo::bank DC_Traj, hipo::bank cherenkov,
		int pbankIndex,
		clashit &electron,
		double starttime , int thisRun , double Ebeam ){

	TVector3	momentum = particles.getV3P(pbankIndex);
	TVector3	vertex	 = particles.getV3v(pbankIndex);
	if( particles.getPid(pbankIndex) != 11 || particles.getCharge(pbankIndex) != -1 ) 	return;

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

	// Cherenkov banks
	int cher_ectr = 0;
	for( int cher_rows = 0 ; cher_rows < cherenkov.getRows() ; ++cher_rows ){
		int cher_index 	= cherenkov.getInt(0,cher_rows);
		int cher_pindex	= cherenkov.getInt(1,cher_rows);
		int cher_det		= cherenkov.getInt(2,cher_rows);
		if( cher_pindex == pbankIndex && cher_det == 15){ // our electron information with the HTCC
			++cher_ectr;
			int cher_sec		= cherenkov.getInt(3,cher_rows);
			double cher_nphe	= cherenkov.getFloat(4,cher_rows);
			double cher_time	= cherenkov.getFloat(5,cher_rows);
			double cher_path	= cherenkov.getFloat(6,cher_rows);
			double cher_chi2	= cherenkov.getFloat(7,cher_rows);
			double cher_x	= cherenkov.getFloat(8,cher_rows);
			double cher_y	= cherenkov.getFloat(9,cher_rows);
			double cher_z	= cherenkov.getFloat(10,cher_rows);
			int cher_status	= cherenkov.getInt(13,cher_rows);

			electron.setNphe	( cher_nphe 	);
			electron.setKov_x	( cher_x	);
			electron.setKov_y	( cher_y	);
			electron.setKov_z	( cher_z	);
			electron.setKov_chi2	( cher_chi2	);
			electron.setKov_time	( cher_time	);
			electron.setKov_path	( cher_path	);
			electron.setKov_det	( cher_det	);
			electron.setKov_sector	( cher_sec	);
			electron.setKov_status	( cher_status	);
		}
	}
	if( cher_ectr > 1 ){ cout << "ERROR in readhipo_helper::getElectronInfo with cherenkov having more than one entry\n"; exit(-1); }

	// Scntillator banks for electron
	for( int sci_rows = 0 ; sci_rows < scintillator.getRows() ; ++sci_rows ){
		int sci_index		= scintillator.getInt(0,sci_rows);
		int sci_pindex		= scintillator.getInt(1,sci_rows);
		int sci_det		= scintillator.getInt(2,sci_rows);
		if( sci_pindex == pbankIndex && sci_det == 12 ){ // our electron information in FTOF
			int sci_sec	= scintillator.getInt(3,sci_rows);
			int sci_lay	= scintillator.getInt(4,sci_rows);
			int sci_comp	= scintillator.getInt(5,sci_rows);

			double sci_edep = scintillator.getFloat(6,sci_rows);
			double sci_time = scintillator.getFloat(7,sci_rows);
			double sci_path = scintillator.getFloat(8,sci_rows);
			double sci_chi2 = scintillator.getFloat(9,sci_rows);

			double sci_x	= scintillator.getFloat(10,sci_rows);
			double sci_y	= scintillator.getFloat(11,sci_rows);
			double sci_z	= scintillator.getFloat(12,sci_rows);

			int sci_status	= scintillator.getFloat(16,sci_rows);

			electron.setScint_status	( sci_status	);
			electron.setScint_sector	( sci_sec	);
			electron.setScint_layer		( sci_lay	);
			electron.setScint_component	( sci_comp	);
			electron.setScint_Edep		( sci_edep	);
			electron.setScint_time		( sci_time	);
			electron.setScint_path		( sci_path	);
			electron.setScint_chi2		( sci_chi2	);
			electron.setScint_x		( sci_x		);
			electron.setScint_y		( sci_y		);
			electron.setScint_z		( sci_z		);
		}
	}

	if( electron.getScint_sector().size() == 0 ){ // if we don't have any FTOF info, skip the event
		electron.Clear();
		return;
	}

	if( electron.getScint_sector()[0] == -999 || electron.getSector() == -999 || electron.getSector() == -1 ){ // if scint sector != pcal sector
		electron.Clear();
		return;
	}

	// If the sectors do not match, do not store:
	if( 	electron.getScint_sector()[0] 	!= electron.getSector() 	||	// scint sector != pcal sector
			electron.getScint_sector()[0] 	!= electron.getDC_sector() 	||	// scint sector != DC sector
			electron.getScint_sector()[0] 	!= electron.getKov_sector() 	||	// scint sector != htcc sector

			electron.getDC_sector()		!= electron.getSector()		||	// DC sector != pcal sector
			electron.getDC_sector()		!= electron.getKov_sector()	||	// DC sector != htcc sector

			electron.getSector()		!= electron.getKov_sector()	){	// pcal sector != htcc sector


		//particles.show();
		//calorimeter.show();
		//scintillator.show();
		//cherenkov.show();

		electron.Clear(); // do not store an electron
		return;
	}


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

		double beta = nHit[hit].getDL().Mag() / (nHit[hit].getTof()*cAir);
		double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
		nVec.SetMagThetaPhi(p_n,nHit[hit].getDL().Theta(),nHit[hit].getDL().Phi());

		double E_n 	= sqrt( mN*mN + p_n*p_n );
		double W_primeSq = mD*mD - eHit.getQ2() + mN*mN + 2.*mD*(eHit.getOmega()-E_n) - 2.*eHit.getOmega()*E_n + 2.*eHit.getQ()*p_n*cos(theta_nq);
		double Wp = sqrt(W_primeSq);
		double As = (E_n - p_n*CosTheta_nq)/mN;

		// Different definitions of x'
		// Bonus definition is the default
		double Xp = eHit.getQ2()/(2.*( eHit.getOmega()*(mD-E_n) + p_n*eHit.getQ()*CosTheta_nq));
		// W' definition
		double Xp_WP  = eHit.getQ2()/(W_primeSq - mN*mN + eHit.getQ2());
		// Bjorken definition
		double Xp_Bj  = eHit.getXb()/(2. - As);
		// PRC definition
		double Ei = mD - E_n;
		double ps_plus = mD/2. * As;
		double virt = (Ei*Ei - p_n*p_n - mN*mN)/(mN*mN);
		double p_plus = mD - ps_plus;
		double q_plus = eHit.getOmega() - eHit.getQ();
		double tP = virt * mN * mN;
		double Xp_PRC = (eHit.getQ2() - (q_plus/p_plus)*tP)/(W_primeSq - mN*mN + eHit.getQ2() - (q_plus/p_plus)*tP);


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
		tag[hit].setAs		(As		);
		tag[hit].setPt		(Pt		);
		tag[hit].setXp		(Xp		);
		tag[hit].setXp_WP	(Xp_WP		);
		tag[hit].setXp_Bj	(Xp_Bj		);
		tag[hit].setXp_PRC	(Xp_PRC		);
	}

	return;
}


void getParticleInfo( hipo::bank claspart, particles part[maxParticles], hipo::bank scintillator ,int& multiplicity ){
	multiplicity = 0;
	// Loop over particle bank and store info in particle class:
	for( int row = 1; row < claspart.getRows() ; ++row ){
		int charge = claspart.getInt(8,row);
		// Only look at +/- particles (no neutrals):
		if( charge == 1 || charge == -1 ){
			int pid = claspart.getInt(0,row);
			double px = claspart.getFloat(1,row);
			double py = claspart.getFloat(2,row);
			double pz = claspart.getFloat(3,row);
			double vx = claspart.getFloat(4,row);
			double vy = claspart.getFloat(5,row);
			double vz = claspart.getFloat(6,row);
			double vt = claspart.getFloat(7,row);
			double beta = claspart.getFloat(9,row);
			double chi2 = claspart.getFloat(10,row);
			double status = claspart.getInt(11,row);
			TVector3	momentum(px,py,pz);
			TVector3	vertex(vx,vy,vz);

			part[multiplicity].setPID		(	pid 	);
			part[multiplicity].setCharge		(	charge	);
			part[multiplicity].setStatus		(	status	);
			part[multiplicity].setTime		(	vt 	);
			part[multiplicity].setBeta		( 	beta 	);
			part[multiplicity].setChi2		(	chi2	);
			part[multiplicity].setPindex		(	row				); // mapping pindex for other DST banks
			part[multiplicity].setVtx		(	vertex.X()			);
			part[multiplicity].setVty		(	vertex.Y()			);
			part[multiplicity].setVtz		(	vertex.Z()			);
			part[multiplicity].setMomentum		(	momentum.Mag()			);
			part[multiplicity].setTheta		(	momentum.Theta()		);
			part[multiplicity].setPhi		(	momentum.Phi()			);
			++multiplicity;
		}// end if for +/-
	}// end loop over particle bank

	// Loop over scintillator bank to get mapped info
	for( int partidx = 0 ; partidx < multiplicity ; ++partidx ){
		std::vector<int>    indexes	;
		std::vector<int>    dets	;
		std::vector<int>    secs	;
		std::vector<int>    lays	;
		std::vector<int>    comps	;
		std::vector<int>    statuses	;
		std::vector<double> edeps	;
		std::vector<double> times	;
		std::vector<double> paths	;
		std::vector<double> chi2s	;
		std::vector<double> xs		;
		std::vector<double> ys		;
		std::vector<double> zs		;
		int MATCH_PINDEX = part[partidx].getPindex();
		for( int sci_rows = 0 ; sci_rows < scintillator.getRows() ; ++sci_rows ){
			int sci_pindex		= scintillator.getInt(1,sci_rows);
			if( sci_pindex == MATCH_PINDEX ){ 
				int sci_index		= scintillator.getInt(0,sci_rows);
				int sci_det		= scintillator.getInt(2,sci_rows);
				int sci_sec		= scintillator.getInt(3,sci_rows);
				int sci_lay		= scintillator.getInt(4,sci_rows);
				int sci_comp		= scintillator.getInt(5,sci_rows);

				double sci_edep = scintillator.getFloat(6,sci_rows);
				double sci_time = scintillator.getFloat(7,sci_rows);
				double sci_path = scintillator.getFloat(8,sci_rows);
				double sci_chi2 = scintillator.getFloat(9,sci_rows);

				double sci_x	= scintillator.getFloat(10,sci_rows);
				double sci_y	= scintillator.getFloat(11,sci_rows);
				double sci_z	= scintillator.getFloat(12,sci_rows);

				int sci_status	= scintillator.getFloat(16,sci_rows);

				indexes		.push_back(	sci_index	);
				dets		.push_back(	sci_det		);
				secs		.push_back(	sci_sec		);
				lays		.push_back(	sci_lay		);
				comps		.push_back(	sci_comp	);
				statuses	.push_back(	sci_status	);
				edeps		.push_back(	sci_edep	);
				times		.push_back(	sci_time	);
				paths		.push_back(	sci_path	);
				chi2s		.push_back(	sci_chi2	);
				xs		.push_back(	sci_x		);
				ys		.push_back(	sci_y		);
				zs		.push_back(	sci_z		);

			}
		}// end loop over sci rows for 	this particle index
		part[partidx].setScint_detector		( dets		);
		part[partidx].setScint_status  	 	( statuses	);
		part[partidx].setScint_sector  	 	( secs		);
		part[partidx].setScint_layer   	 	( lays		);
		part[partidx].setScint_component	( comps		);
		part[partidx].setScint_Edep		( edeps		);
		part[partidx].setScint_time		( times		);
		part[partidx].setScint_path		( paths		);
		part[partidx].setScint_chi2		( chi2s		);
		part[partidx].setScint_x		( xs		);
		part[partidx].setScint_y		( ys		);
		part[partidx].setScint_z		( zs		);

	}// end loop over particles
}





//Parametrization from Giovanni (GWU) and FX based on double pion analysis in RGA

void smearRGA(TVector3 &vpar){

	TRandom3 myRand(0);

	double mom = vpar.Mag(); // GeV
	double theta = vpar.Theta(); //radians
	double theta_deg = theta * 180 / M_PI; //degree
	double phi = vpar.Phi(); //radians

	//Determine mom resolution value dependent on theta
	double mom_S1 = 0.0184291 -0.0110083*theta_deg + 0.00227667*theta_deg*theta_deg -0.000140152*theta_deg*theta_deg*theta_deg + 3.07424e-06*theta_deg*theta_deg*theta_deg*theta_deg;
	double mom_S2 = 0.02*theta_deg;
	double mom_R = 0.01 * sqrt( pow(mom_S1*mom, 2) + pow(mom_S2, 2));
	mom_R *= 2.0; // <- to match data resolution, value on % of momentum

	//Determine theta resolution dependent on mom
	double theta_S1 = 0.004*theta_deg + 0.1;
	double theta_S2 = 0;
	double theta_R = sqrt(pow(theta_S1*sqrt(mom*mom+0.13957*0.13957)/(mom*mom),2) + pow(theta_S2,2));
	theta_R *= 2.5; // <- to match data resolution. value in Degree

	//Determine phi resolution dependent on theta and momentum
	double phi_S1 = 0.85-0.015*theta_deg;
	double phi_S2 = 0.17-0.003*theta_deg;
	double phi_R = sqrt(pow(phi_S1*sqrt(mom*mom+0.13957*0.13957)/(mom*mom),2) + pow(phi_S2,2) );
	phi_R *= 3.5; // <- to match data resolution, value in degree

	//Add smearing to mom, theta and phi values based on Gaussian distributions
	mom += myRand.Gaus(0,mom_R * mom);
	theta += myRand.Gaus(0,theta_R * (M_PI / 180));
	phi += myRand.Gaus(0,phi_R * (M_PI / 180));
	vpar.SetMagThetaPhi(mom,theta,phi);
}

void recalculate_clashit_kinematics(clashit &input_ehit, double Ebeam, TVector3 &smeared_electron) {

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	qVec; qVec = beamVec - smeared_electron;

	input_ehit.setMomentum		(	smeared_electron.Mag()						);
	input_ehit.setTheta		  (	smeared_electron.Theta()					);
	input_ehit.setPhi			  (	smeared_electron.Phi()						);

	input_ehit.setQ			  (	qVec.Mag()						);
	input_ehit.setThetaQ		(	qVec.Theta()						);
	input_ehit.setPhiQ		  (	qVec.Phi()						);

	input_ehit.setOmega		(	Ebeam - sqrt( pow(smeared_electron.Mag(),2) + mE*mE )		);
	input_ehit.setQ2			  (	qVec.Mag()*qVec.Mag() - pow(input_ehit.getOmega(),2)	);
	input_ehit.setXb			  (	input_ehit.getQ2()/(2.*mP*input_ehit.getOmega())		);
	input_ehit.setW2		  	(	mP*mP - input_ehit.getQ2() + 2.*input_ehit.getOmega()*mP	);

}

bool electron_fiducials( const int period , clashit * const eHit ){
	int sector 	= eHit->getSector();
	double v	= eHit->getV();
	double w	= eHit->getW();
	//double vtz	= eHit->getVtz();
	//double pe	= eHit->getMomentum();
	//double Q2	= eHit->getQ2();
	//double W2	= eHit->getW2();
	//double EoP	= eHit->getEoP();
	//double Epcal	= eHit->getEpcal();
	std::vector<int>	scint_sec = eHit->getScint_sector();
	std::vector<int>	scint_lay = eHit->getScint_layer();
	std::vector<int>	scint_com = eHit->getScint_component();
	
	// Lose final cuts to pare down file size:
	//if( v < 10 || w < 10 ) 					return false;
	//if( vtz < -8 || vtz > 2 )				return false;
	//if( pe < 2 )						return false;
	//if( Q2 < 1.5 )						return false;
	//if( W2 < 1.6*1.6 )					return false;
	//if( EoP < 0.12 || EoP > 0.35 )				return false;
	//if( Epcal < 0.05 )					return false;


	if( period == 1 ){	// 10.2 fiducial checks

		if( sector == 2 && ( (v > 30 && v < 55) || (v>95 && v < 120) ) ) 	return false;
		if( sector == 1 && w > 70 && w < 100 ) 					return false;

		bool sector2_layer2_hit = false;
		for( int scint_hit = 0 ; scint_hit < scint_sec.size() ; ++scint_hit ){
			if( scint_sec[scint_hit] == 2 && scint_lay[scint_hit] == 2 ) sector2_layer2_hit = true;
		}
		for( int scint_hit = 0 ; scint_hit < scint_sec.size() ; ++scint_hit ){
			if( scint_sec[scint_hit] == 2 && scint_lay[scint_hit] == 1 && (scint_com[scint_hit] == 6 || scint_com[scint_hit] == 10) && !sector2_layer2_hit ){
				return false;
			}
			if( scint_sec[scint_hit] == 5 && scint_lay[scint_hit] == 2 && (scint_com[scint_hit] == 12 || scint_com[scint_hit] == 13) ){
				return false;
			}
		}
	}

	return true;
}
