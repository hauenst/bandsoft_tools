#include "bandreco.h"

void BANDReco::Clear(){
	candidate_pmts.clear();
	candidate_bars.clear();
	return;
}

void BANDReco::setPeriod( const int period ){
	(period == 0) ? SPRING2019 = true : FALL2019_WINTER2020 = true;
	
	// Set the vertex offset of the target:
	if( SPRING2019 ){
		TRGT_VERTEX_OFFSET = -3; // [cm]
	}
	else if(FALL2019_WINTER2020){
		TRGT_VERTEX_OFFSET = -3; // [cm]
	}

	return;
}

void BANDReco::readTW(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;
	// Load L PMT file:
	f.open(path+"/time_walk_amp_L.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			double f01,f02,f03,f04,f05,f06;
			double f11,f12,f13,f14,f15,f16;
			double f21,f22,f23,f24,f25,f26;
			double f31,f32,f33,f34,f35,f36;
			f01 = f02 = f03 = f04 = f05 = f06 =
			f11 = f12 = f13 = f14 = f15 = f16 =
			f21 = f22 = f23 = f24 = f25 = f26 =
			f31 = f32 = f33 = f34 = f35 = f36 = 0;
			ss >> sector >> layer >> component 
				>> f01 >> f02 >> f03 >> f04 >> f05 >> f06 
				>> f11 >> f12 >> f13 >> f14 >> f15 >> f16
				>> f21 >> f22 >> f23 >> f24 >> f25 >> f26
				>> f31 >> f32 >> f33 >> f34 >> f35 >> f36;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ f01, f02, f03, f04, f05, f06,
		       				f11, f12, f13, f14, f15, f16,
						f21, f22, f23, f24, f25, f26,
						f31, f32, f33, f34, f35, f36	};
			TWParamsAMP[PMTID] = temp;
		}
	}
	f.close();
	// Load R PMT file:
	f.open(path+"/time_walk_amp_R.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			double f01,f02,f03,f04,f05,f06;
			double f11,f12,f13,f14,f15,f16;
			double f21,f22,f23,f24,f25,f26;
			double f31,f32,f33,f34,f35,f36;
			f01 = f02 = f03 = f04 = f05 = f06 =
			f11 = f12 = f13 = f14 = f15 = f16 =
			f21 = f22 = f23 = f24 = f25 = f26 =
			f31 = f32 = f33 = f34 = f35 = f36 = 0;
			ss >> sector >> layer >> component 
				>> f01 >> f02 >> f03 >> f04 >> f05 >> f06 
				>> f11 >> f12 >> f13 >> f14 >> f15 >> f16
				>> f21 >> f22 >> f23 >> f24 >> f25 >> f26
				>> f31 >> f32 >> f33 >> f34 >> f35 >> f36;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ f01, f02, f03, f04, f05, f06,
		       				f11, f12, f13, f14, f15, f16,
						f21, f22, f23, f24, f25, f26,
						f31, f32, f33, f34, f35, f36	};
			TWParamsAMP[PMTID] = temp;
		}
	}
	f.close();
	return;
}

void BANDReco::readLROffset(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;

	// Load L-R offset file:
	f.open(path+"/lr_offsets.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			double tdc_off, tdc_veff, tdc_width;
			double ftdc_off, ftdc_veff, ftdc_width;
			tdc_off = tdc_veff = tdc_width = 
			ftdc_off = ftdc_veff = ftdc_width = 0;
			ss >> sector >> layer >> component >>
				tdc_off >> tdc_veff >> tdc_width >>
				ftdc_off >> ftdc_veff >> ftdc_width;

			int BARID = sector*100 + layer*10 + component;
			TDCOffsets[BARID] = tdc_off;
			TDCVelocity[BARID] = tdc_veff;
			FTDCOffsets[BARID] = ftdc_off;
			FTDCVelocity[BARID] = ftdc_veff;
		}
	}
	f.close();
	return;
}

void BANDReco::readPaddleOffset(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;

	// Load Paddle offset file:
	f.open(path+"/paddle_offsets.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector,layer,component;
			double tdc_amp, tdc_mean, tdc_sigma;
			double ftdc_amp, ftdc_mean, ftdc_sigma;
			sector = layer = component = 
				tdc_amp = tdc_mean = tdc_sigma =
				ftdc_amp = ftdc_mean = ftdc_sigma = 0;
			ss >> sector >> layer >> component >> 
				tdc_amp >> tdc_mean >> tdc_sigma >>
				ftdc_amp >> ftdc_mean >> ftdc_sigma;

			int BARID = sector*100 + layer*10 + component;
			TDCPaddle[BARID] = tdc_mean;
			FTDCPaddle[BARID] = ftdc_mean;
		}
	}
	f.close();
	return;
}

void BANDReco::readGeometry(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;

	// Load the Geometry positions:
	f.open(path+"/band_geometry.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer, component;
			double x,y,z;
			sector = layer = component =
				x = y = z = 0;
			ss >> sector >> layer >> component >> x >> y >> z;

			int BARID = sector*100 + layer*10 + component;
			GlobalX[BARID] = x;
			GlobalY[BARID] = y;
			GlobalZ[BARID] = z;

		}
	}
	f.close();

	return;
}

void BANDReco::readEnergyCalib(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;

	// Load the Edep calibration:
	f.open(path+"/edep_calibration.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer, component;
			double AdcToMeVee;
			sector = layer = component =
				AdcToMeVee = 0;
			ss >> sector >> layer >> component >> AdcToMeVee;

			int BARID = sector*100 + layer*10 + component;
			ADCtoMEV[BARID] = AdcToMeVee;
		}
	}
	f.close();

	return;
}

double BANDReco::getTriggerPhase( const long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}

double BANDReco::timewalk( const double *x , const double *p){
	double A = *x;
	double f0 = p[0] + p[1]/pow(A,p[2]); 				// p[3]-p[5] not used
	double f1 = p[6]/A + p[7]/(exp((A-p[8])/p[9]) + p[10] );	// p[11] not used
	double f2 = p[12] + p[13]*A;					// p[14]-p[17] not used
	double f3 = p[18]*sin((A-p[19])/p[20]) + p[21]/pow(A,p[22]);	// p[23] not used
	if( p[20] == 0 ) f3 = 0;

	if( f0 != f0 || f1 != f1 || f2 != f2 || f3 != f3 ){ cerr << "issue with timewalk\n"; exit(-1); }

	return f0+f1+f2+f3;
}

bool BANDReco::check_bar(Bar * this_bar){
	PMT left 	= this_bar->left;
	PMT right 	= this_bar->right;
	if( left.layer == 6 || right.layer == 6 ) return true;
	if( left.amp == 0 || right.amp == 0 		) return false;
	if( left.adc == 0 || right.adc == 0		) return false;
	if( left.ftdc == 0 || right.ftdc == 0		) return false;
	if( left.tdc_corr == 0 || right.tdc_corr == 0 	) return false;
	if( left.ftdc_corr == 0 || right.ftdc_corr == 0 ) return false;
	if( left.tdc == 0 || right.tdc == 0 		) return false;
	if( left.sector != right.sector			) return false;
	if( left.layer != right.layer			) return false;
	if( left.component != right.component		) return false;
	if( abs(left.PMT_ID-right.PMT_ID)>1		) return false;

	return true;
}


void BANDReco::createPMTs( const hipo::bank * band_adc , const hipo::bank * band_tdc , const hipo::bank * run_config ){

	// Get the phase correction for the TDC
	long timestamp = run_config->getLong(4 , 0 );
	double phaseCorr = getTriggerPhase(timestamp);
	
	// Loop over all the ADC hits:
	map<int,vector<PMT>> raw_pmts_adc;
	map<int,vector<PMT>> raw_pmts_tdc;
	for( int row = 0 ; row < band_adc->getRows() ; ++row ){
		int sector 	= band_adc->getInt( 0, row );
		int layer 	= band_adc->getInt( 1, row );
		int component 	= band_adc->getInt( 2, row );
		int order	= band_adc->getInt( 3, row );
		int PMT_ID	= sector*1000 + layer*100 + component*10 + order;

		int adc		= band_adc->getInt( 4, row );
		int amp		= band_adc->getInt( 5, row );
		double ftdc	= band_adc->getFloat( 6, row );
		int ped		= band_adc->getInt( 7, row );

		if( amp < 250 || amp >= 4095 ) continue;  // cut for TW 

		PMT this_pmt;
		this_pmt.PMT_ID 	= PMT_ID;
		this_pmt.sector		= sector;
		this_pmt.layer		= layer;
		this_pmt.component	= component;
		this_pmt.order		= order;
		this_pmt.adc		= adc;
		this_pmt.amp		= amp;
		this_pmt.ftdc		= ftdc;
		this_pmt.ped		= ped;
		this_pmt.idx_ftdc	= row;
		raw_pmts_adc[PMT_ID].push_back(this_pmt);
	}
	// Loop over all the TDC hits to find a match:
	for( int row = 0 ; row < band_tdc->getRows() ; ++row ){
		int sector 	= band_tdc->getInt( 0 , row );
		int layer 	= band_tdc->getInt( 1 , row );
		int component 	= band_tdc->getInt( 2 , row );
		int order	= band_tdc->getInt( 3 , row ) - 2;
		int PMT_ID	= sector*1000 + layer*100 + component*10 + order;

		double tdc	= band_tdc->getInt( 4 , row ) * 0.02345 - phaseCorr;

		// If we do not have a matching FTDC hit for the PMT, we 
		// skip this event completely
		if( raw_pmts_adc[PMT_ID].size() == 0 ) continue;

		PMT this_pmt;
		this_pmt.PMT_ID 	= PMT_ID;
		this_pmt.sector		= sector;
		this_pmt.layer		= layer;
		this_pmt.component	= component;
		this_pmt.order		= order;
		this_pmt.tdc		= tdc;
		this_pmt.trigphase	= phaseCorr;
		this_pmt.idx_tdc	= row;
		raw_pmts_tdc[PMT_ID].push_back(this_pmt);
	}	
	
	// Now let's loop over all the FTDC PMT hits to do an additional filter.
	// and then add a valid PMT to our bar list
	map<int,vector<PMT>>::iterator it;
	int multiTdc = 0;
	//map<int,PMT> candidate_pmts;
	for( it = raw_pmts_adc.begin() ; it != raw_pmts_adc.end() ; ++it){
		int PMT_ID = it->first;

		// If we don't have a TDC, we can't do anything
		if( raw_pmts_tdc[PMT_ID].size() == 0 ) continue;

		// If we have more than one ADC, let's not consider this
		if( it->second.size() > 1 ){
			cerr << "unexpected behavior with double adc hit for single pmt. exiting...\n";
			exit(-1);
		}
		
		// If we have more than one TDC, let's not deal with this for the moment
		if( raw_pmts_tdc[PMT_ID].size() > 1 ){
			++multiTdc;
			continue;
		}

		// If we only have 1 of each, so let's store this information and
		// create a BAND PMT candidate:
		if( it->second.size() == 1 && raw_pmts_tdc[PMT_ID].size() == 1){
			PMT this_pmt	= it->second[0];
			PMT candidate;
			candidate.sector 	= this_pmt.sector;
			candidate.layer		= this_pmt.layer;
			candidate.component	= this_pmt.component;
			candidate.order		= this_pmt.order;
			candidate.adc		= this_pmt.adc;
			candidate.amp		= this_pmt.amp;
			candidate.ftdc		= this_pmt.ftdc;
			candidate.ped		= this_pmt.ped;
			candidate.idx_ftdc	= this_pmt.idx_ftdc;
			candidate.tdc		= raw_pmts_tdc[PMT_ID][0].tdc;
			candidate.trigphase	= raw_pmts_tdc[PMT_ID][0].trigphase;
			candidate.idx_tdc	= raw_pmts_tdc[PMT_ID][0].idx_tdc;

			double tdc 		= candidate.tdc;
			double * tw_pars 	= &TWParamsAMP[PMT_ID][0];
			double amp		= candidate.amp;
			double tdc_ampcorr 	= tdc - timewalk(&amp,tw_pars);
			candidate.tdc_corr	= tdc_ampcorr;				// store TW corrected time
			candidate.ftdc_corr	= candidate.ftdc;			// no TW correction for FTDC

			candidate_pmts[PMT_ID] = candidate;
		}
	}
	
	
	return;
}


void BANDReco::createBars(  ){
	map<int,PMT>::const_iterator pmt_i;
	map<int,PMT>::const_iterator pmt_j;

	for( pmt_i = candidate_pmts.begin() ; pmt_i != candidate_pmts.end() ; ++pmt_i ){
		int PMT_ID1 	= pmt_i->first;
		PMT candidate1 	= pmt_i->second;

		if( candidate1.layer == 6 ){
			// set all bar quantities as just 1 sided, same for each
			Bar newbar;
			int Bar_ID = candidate1.sector*100 + candidate1.layer*10 + candidate1.component;
			newbar.left 	= candidate1;
			newbar.right 	= candidate1;
			if( check_bar(&newbar) == false ) continue;	// make sure all the values for left and right are set
			newbar.sector 		= candidate1.sector;
			newbar.layer		= candidate1.layer;
			newbar.component	= candidate1.component;
			newbar.Bar_ID		= Bar_ID;

			double tdiff		= (newbar.left.tdc_corr - newbar.right.tdc_corr) - TDCOffsets[Bar_ID];
			double tdiff_ftdc	= (newbar.left.ftdc_corr - newbar.right.ftdc_corr) - FTDCOffsets[Bar_ID];

			double meantime		= 0.5*(newbar.left.tdc_corr + newbar.right.tdc_corr - TDCOffsets[Bar_ID]	) - TDCPaddle[Bar_ID];
			double meantime_ftdc	= 0.5*(newbar.left.ftdc_corr + newbar.right.ftdc_corr - FTDCOffsets[Bar_ID]	) - FTDCPaddle[Bar_ID];

			double x = GlobalX[Bar_ID];
			double y = GlobalY[Bar_ID];
			double z = GlobalZ[Bar_ID] - TRGT_VERTEX_OFFSET;

			double sec_length = bar_lengths[newbar.sector-1];
			double mu = 1e10;
			double adcL_corr = newbar.left.adc * exp( (sec_length/2.-x)/mu );
			double adcR_corr = newbar.right.adc * exp( (sec_length/2.-x)/mu );
			double ampL_corr = newbar.left.amp * exp( (sec_length/2.-x)/mu );
			double ampR_corr = newbar.right.amp * exp( (sec_length/2.-x)/mu );

			double edep = sqrt( adcL_corr * adcR_corr );
			edep /= ADCtoMEV[Bar_ID]; // [ADC] / ([Adc/MeV]) = [MeV]

			newbar.Tof 		= meantime;
			newbar.TofFtdc 		= meantime_ftdc;
			newbar.Tdiff		= tdiff;
			newbar.TdiffFtdc	= tdiff_ftdc;
			newbar.Edep		= edep;
			newbar.AdcL		= adcL_corr;
			newbar.AdcR		= adcR_corr;
			newbar.AmpL		= ampL_corr;
			newbar.AmpR		= ampR_corr;
			newbar.X		= x;
			newbar.Y		= y;
			newbar.Z		= z;
			candidate_bars[Bar_ID] = newbar;
			continue;
		}

		int sector 	= candidate1.sector;
		int layer	= candidate1.layer;
		int component	= candidate1.component;
		int order	= candidate1.order;
		int Bar_ID	= sector*100 + layer*10 + component;

		// Loop over the remaining candidates and try to match the PMTs
		for( pmt_j = pmt_i ; pmt_j != candidate_pmts.end() ; ++pmt_j ){
			int PMT_ID2	= pmt_j->first;
			PMT candidate2	= pmt_j->second;
			if( abs(PMT_ID1 - PMT_ID2) != 1 ) continue;

			Bar newbar;
			if( candidate1.order == 0 ){ // candidate 1 is left, candidate 2 is right
				newbar.left 	= candidate1;
				newbar.right 	= candidate2;
			}
			else if( candidate2.order == 1 ){ // candidate 1 is right, candidate 2 is left
				newbar.left 	= candidate2;
				newbar.right 	= candidate1;
			}
			else{ cerr << "what\n"; continue;}
			newbar.sector 		= sector;
			newbar.layer		= layer;
			newbar.component	= component;
			newbar.Bar_ID		= Bar_ID;
			if( check_bar(&newbar) == false ) continue;	// make sure all the values for left and right are set

			// Now we can check to make sure we have a valid bar and set the quantities correctly
			double tdiff		= (newbar.left.tdc_corr - newbar.right.tdc_corr) - TDCOffsets[Bar_ID];
			double tdiff_ftdc	= (newbar.left.ftdc_corr - newbar.right.ftdc_corr) - FTDCOffsets[Bar_ID];
			double maxTdiff		= bar_lengths[sector-1] / TDCVelocity[Bar_ID];
			double maxTdiff_ftdc	= bar_lengths[sector-1] / FTDCVelocity[Bar_ID];
			if( fabs(tdiff) > maxTdiff 		) continue;
			if( fabs(tdiff_ftdc) > maxTdiff_ftdc 	) continue;
	
			// Create the ToF:
			double meantime		= 0.5*(newbar.left.tdc_corr + newbar.right.tdc_corr - TDCOffsets[Bar_ID]	) - TDCPaddle[Bar_ID];
			double meantime_ftdc	= 0.5*(newbar.left.ftdc_corr + newbar.right.ftdc_corr - FTDCOffsets[Bar_ID]	) - FTDCPaddle[Bar_ID];

			// Create the x,y,z:
			double x = GlobalX[Bar_ID] - (0.5 * tdiff * TDCVelocity[Bar_ID]);
			double y = GlobalY[Bar_ID];
			double z = GlobalZ[Bar_ID] - TRGT_VERTEX_OFFSET;

			// Create the Edep:
			double sec_length = bar_lengths[newbar.sector-1];
			double mu = 1e10;
			double adcL_corr = newbar.left.adc * exp( (sec_length/2.-x)/mu );
			double adcR_corr = newbar.right.adc * exp( (sec_length/2.-x)/mu );
			double ampL_corr = newbar.left.amp * exp( (sec_length/2.-x)/mu );
			double ampR_corr = newbar.right.amp * exp( (sec_length/2.-x)/mu );

			double edep = sqrt( adcL_corr * adcR_corr );
			edep /= ADCtoMEV[Bar_ID]; // [ADC] / ([Adc/MeV]) = [MeV]

			newbar.Tof 		= meantime;
			newbar.TofFtdc 		= meantime_ftdc;
			newbar.Tdiff		= tdiff;
			newbar.TdiffFtdc	= tdiff_ftdc;
			newbar.Edep		= edep;
			newbar.AdcL		= adcL_corr;
			newbar.AdcR		= adcR_corr;
			newbar.AmpL		= ampL_corr;
			newbar.AmpR		= ampR_corr;
			newbar.X		= x;
			newbar.Y		= y;
			newbar.Z		= z;
			candidate_bars[Bar_ID] = newbar;
		}
	}

	return;
}

void BANDReco::storeHits( int& mult , bandhit * hits , const double starttime ){

	mult = 0;
	map<int,Bar>::iterator bar_it;
	for( bar_it = candidate_bars.begin() ; bar_it != candidate_bars.end() ; ++bar_it){

		int Bar_ID 	= bar_it->first;
		Bar this_bar 	= bar_it->second;

		// If Edep < 2, let's not store the hit:
		if( this_bar.Edep < 2 && this_bar.layer != 6 ) continue;

		// Set the hits information
		hits[mult].setSector		(this_bar.sector);
		hits[mult].setLayer		(this_bar.layer);
		hits[mult].setComponent		(this_bar.component);
		hits[mult].setBarID		(Bar_ID);
		hits[mult].setEdep		(this_bar.Edep);
		hits[mult].setTof		(this_bar.Tof 		- starttime );
		hits[mult].setTofFadc		(this_bar.TofFtdc 	- starttime );
		hits[mult].setTdiff		(this_bar.Tdiff);
		hits[mult].setTdiffFadc		(this_bar.TdiffFtdc);
		hits[mult].setX			(this_bar.X);
		hits[mult].setY			(this_bar.Y);
		hits[mult].setZ			(this_bar.Z);

		// 	Get the raw hit information corresponding to the band hit above
		hits[mult].setPmtLtdc		(this_bar.left.tdc_corr);
		hits[mult].setPmtRtdc		(this_bar.right.tdc_corr);
		hits[mult].setPmtLtfadc		(this_bar.left.ftdc_corr);
		hits[mult].setPmtRtfadc		(this_bar.right.ftdc_corr);
		hits[mult].setPmtLamp		(this_bar.left.amp);
		hits[mult].setPmtRamp		(this_bar.right.amp);
		hits[mult].setPmtLadc		(this_bar.left.adc);
		hits[mult].setPmtRadc		(this_bar.right.adc);
		hits[mult].setPmtLped		(this_bar.left.ped);
		hits[mult].setPmtRped		(this_bar.right.ped);

		hits[mult].setStatus		(0.);
		hits[mult].setDL		( TVector3(hits[mult].getX(),hits[mult].getY(),hits[mult].getZ()) );

		// If layer == 6 and sector == 4, take the idxR, otherwise take idxL for the veto bars:
		//if( hits[mult].getLayer() == 6 ){
			//hits[hit].setX			( bar_x[band_hits.getBarKey		(hit)	]	);
			//int rawhit_idx = -1;
			//if( hits[hit].getSector() == 4 ) rawhit_idx = band_hits.getRpmtindex(hit);
			//else{ rawhit_idx = band_hits.getLpmtindex(hit); }

			//hits[hit].setRawLtdc		(band_rawhits.getFloat( 7 , rawhit_idx ) 		);
			//hits[hit].setRawLtdccorr	(band_rawhits.getFloat( 9 , rawhit_idx ) 		);
			//hits[hit].setRawLtfadc		(band_rawhits.getFloat( 8 , rawhit_idx ) 		);
			//hits[hit].setRawLamp		(band_rawhits.getFloat( 6 , rawhit_idx )		);
			//hits[hit].setRawLadc		(band_rawhits.getFloat( 5 , rawhit_idx )		);

			//hits[hit].setRawRtdc		(band_rawhits.getFloat( 7 , rawhit_idx ) 		);
			//hits[hit].setRawRtdccorr	(band_rawhits.getFloat( 9 , rawhit_idx ) 		);
			//hits[hit].setRawRtfadc		(band_rawhits.getFloat( 8 , rawhit_idx ) 		);
			//hits[hit].setRawRamp		(band_rawhits.getFloat( 6 , rawhit_idx )		);
			//hits[hit].setRawRadc		(band_rawhits.getFloat( 5 , rawhit_idx )		);

			////Use average Edep per channel for vetos (could be either RawLadc or RawRadc)
			//hits[hit].setEdep		(hits[hit].getRawLadc()	);

			//// Using the rawhit struct, get the raw PMT information to use later
			//int pmtTdc	= band_rawhits.getInt( 10 , rawhit_idx );
			//int pmtAdc	= band_rawhits.getInt( 11 , rawhit_idx );
			////	Get the raw pmt information corresponding to the band hit above
			//hits[hit].setPmtLtdc		(band_tdc.getInt( 4 	, pmtTdc )		);
			//hits[hit].setPmtRtdc		(band_tdc.getInt( 4 	, pmtTdc )		);
			//hits[hit].setPmtLtfadc		(band_adc.getFloat( 6 	, pmtAdc )		);
			//hits[hit].setPmtRtfadc		(band_adc.getFloat( 6 	, pmtAdc )		);
			//hits[hit].setPmtLamp		(band_adc.getInt( 5 	, pmtAdc )		);
			//hits[hit].setPmtRamp		(band_adc.getInt( 5 	, pmtAdc )		);
			//hits[hit].setPmtLadc		(band_adc.getInt( 4 	, pmtAdc )		);
			//hits[hit].setPmtRadc		(band_adc.getInt( 4 	, pmtAdc )		);
			//hits[hit].setPmtLped		(band_adc.getInt( 7 	, pmtAdc )		);
			//hits[hit].setPmtRped		(band_adc.getInt( 7 	, pmtAdc )		);
		//}
		//else{
			// Using the band hit struct, get the raw hit PMT information to use later
			//int rawhit_idxL = band_hits.getLpmtindex(hit);
			//int rawhit_idxR = band_hits.getRpmtindex(hit);
			//// 	Get the raw hit information corresponding to the band hit above
			//hits[hit].setRawLtdc		(band_rawhits.getFloat( 7 , rawhit_idxL ) 		);
			//hits[hit].setRawRtdc		(band_rawhits.getFloat( 7 , rawhit_idxR ) 		);
			//hits[hit].setRawLtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxL ) 		);
			//hits[hit].setRawRtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxR ) 		);
			//hits[hit].setRawLtfadc		(band_rawhits.getFloat( 8 , rawhit_idxL ) 		);
			//hits[hit].setRawRtfadc		(band_rawhits.getFloat( 8 , rawhit_idxR ) 		);
			//hits[hit].setRawLamp		(band_rawhits.getFloat( 6 , rawhit_idxL )		);
			//hits[hit].setRawRamp		(band_rawhits.getFloat( 6 , rawhit_idxR )		);
			//hits[hit].setRawLadc		(band_rawhits.getFloat( 5 , rawhit_idxL )		);
			//hits[hit].setRawRadc		(band_rawhits.getFloat( 5 , rawhit_idxR )		);

			//// Using the rawhit struct, get the raw PMT information to use later
			//int pmtTdcL	= band_rawhits.getInt( 10 , rawhit_idxL );
			//int pmtAdcL	= band_rawhits.getInt( 11 , rawhit_idxL );
			//int pmtTdcR	= band_rawhits.getInt( 10 , rawhit_idxR );
			//int pmtAdcR	= band_rawhits.getInt( 11 , rawhit_idxR );
			////	Get the raw pmt information corresponding to the band hit above
			//hits[hit].setPmtLtdc		(band_tdc.getInt( 4 , pmtTdcL )		);
			//hits[hit].setPmtRtdc		(band_tdc.getInt( 4 , pmtTdcR )		);
			//hits[hit].setPmtLtfadc		(band_adc.getFloat( 6 , pmtAdcL )	);
			//hits[hit].setPmtRtfadc		(band_adc.getFloat( 6 , pmtAdcR )	);
			//hits[hit].setPmtLamp		(band_adc.getInt( 5 , pmtAdcL )		);
			//hits[hit].setPmtRamp		(band_adc.getInt( 5 , pmtAdcR )		);
			//hits[hit].setPmtLadc		(band_adc.getInt( 4 , pmtAdcL )		);
			//hits[hit].setPmtRadc		(band_adc.getInt( 4 , pmtAdcR )		);
			//hits[hit].setPmtLped		(band_adc.getInt( 7 , pmtAdcL )		);
			//hits[hit].setPmtRped		(band_adc.getInt( 7 , pmtAdcR )		);
		//}

		//If Edep maps are not empty, use them to calibrate Edep in MeV for each bar
		//For MC Edep other file is used with constants
		//if ( !bar_edep.empty() ) {
			//Calibration values are in channel/MeV
			//cout << "applying energycorrection before value: " << hits[hit].getEdep() << endl;
			//	hits[hit].setEdep		(band_hits.getEnergy	(hit)	/ bar_edep[band_hits.getBarKey	(hit)	]		);
			//		cout << "applying energycorrection after value: " << hits[hit].getEdep() << endl;
		//}

		// After all corrections to positions, set the pathlength and status:
		//hits[hit].setStatus		(band_hits.getStatus		(hit)			);
		//hits[hit].setDL			( TVector3(hits[hit].getX(),hits[hit].getY(),hits[hit].getZ()) );
		// Save how many neutron hits we have
		mult++;
	}
	return;
}
