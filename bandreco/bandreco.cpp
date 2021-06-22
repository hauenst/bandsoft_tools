#include "bandreco.h"

void BANDReco::Clear(){
	candidate_pmts.clear();
	candidate_bars.clear();
	return;
}

void BANDReco::Print(){
	map<int,Bar>::const_iterator bar_i;
	for( bar_i = candidate_bars.begin() ; bar_i != candidate_bars.end() ; ++bar_i ){
		int Bar_ID = bar_i->first;
		Bar bar = bar_i->second;
	
		std::cout << "Bar ID: " << Bar_ID << " " << bar.Bar_ID << "\n";
		std::cout << "\tS L C: " << bar.sector << " " << bar.layer << " " << bar.component << "\n";
		std::cout << "\tEdep: " << bar.Edep << "\n";
		std::cout << "\tAdcL, AdcR: " << bar.AdcL << " " << bar.AdcR << "\n";
		std::cout << "\tAmpL, AmpR: " << bar.AmpL << " " << bar.AmpR << "\n";
		std::cout << "\tToF, ToFFadc: " << bar.Tof << " " << bar.TofFtdc << "\n";
		std::cout << "\tx,y,z: " << bar.X << " " << bar.Y << " " << bar.Z << "\n";
		std::cout << "\txFtdc: " << bar.XFtdc << "\n";
		std::cout << "\tTdiff: " << bar.Tdiff << " " << bar.TdiffFtdc << "\n";

		std::cout << "\tPmtLadc: " 		<< bar.left.adc << "\n";
		std::cout << "\tPmtLamp: " 		<< bar.left.amp << "\n";
		std::cout << "\tPmtLftdc: " 		<< bar.left.ftdc << "\n";
		std::cout << "\tPmtLftdc_corr: " 	<< bar.left.ftdc_corr << "\n";
		std::cout << "\tPmtLped: " 		<< bar.left.ped << "\n";
		std::cout << "\tPmtLtdc: " 		<< bar.left.tdc << "\n";
		std::cout << "\tPmtLtdc_corr: " 	<< bar.left.tdc_corr << "\n";
		std::cout << "\tPmtRadc: " 		<< bar.right.adc << "\n";
		std::cout << "\tPmtRamp: " 		<< bar.right.amp << "\n";
		std::cout << "\tPmtRftdc: " 		<< bar.right.ftdc << "\n";
		std::cout << "\tPmtRftdc_corr: " 	<< bar.right.ftdc_corr << "\n";
		std::cout << "\tPmtRped: " 		<< bar.right.ped << "\n";
		std::cout << "\tPmtRtdc: " 		<< bar.right.tdc << "\n";
		std::cout << "\tPmtRtdc_corr: " 	<< bar.right.tdc_corr << "\n\n";
	}
	

	return;
}

void BANDReco::setPeriod( const int period ){
	(period == 0) ? SPRING2019 = true : FALL2019_WINTER2020 = true;
	
	// Set the vertex offset of the target:
	if( SPRING2019 ){
		BAND_MOTHER_OFFSET = -3; // [cm]
	}
	else if(FALL2019_WINTER2020){
		BAND_MOTHER_OFFSET = -3; // [cm]
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

	loaded_TW = true;
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

void BANDReco::readLayerOffset(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;

	// Load Layer offset file:
	int LAYERREF = -1;
	f.open(path+"/layer_offsets.txt");
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
			if( tdc_mean != 0 && ftdc_mean != 0 ) LAYERREF = sector*100 + component;
			TDCLayer[BARID] = tdc_mean;
			FTDCLayer[BARID] = ftdc_mean;
		}
	}
	f.close();
	
	// Now what we need to do is go through and overwrite the offsets to be the same in each
	// layer based on the reference we have in the file:
	map<int,double>::iterator bar_it;
	for( bar_it = TDCLayer.begin() ; bar_it != TDCLayer.end() ; ++bar_it ){
		int Bar_ID = bar_it->first;
		int layer = int((Bar_ID - int(Bar_ID/100)*100)/10);
		TDCLayer[Bar_ID] 	= TDCLayer[LAYERREF+layer*10];
		FTDCLayer[Bar_ID] 	= FTDCLayer[LAYERREF+layer*10];
	}

	return;
}

void BANDReco::readGlobalOffset(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;
	f.open(path+"/global_offsets_tdc.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer, component;
			double bkg, amp, mean, sigma, integral, flag;
			sector = layer = component = 
				bkg = amp = mean = sigma = integral = flag = 0;
			ss >> sector >> layer >> component >>
				bkg >> amp >> mean >> sigma >> integral >> flag;

			int BAR_ID = sector*100+layer*10+component;

			TDCGlobal[BAR_ID] = mean;
			TDCToFRes[BAR_ID] = sigma;
		}
	}
	f.close();

	f.open(path+"/global_offsets_fadc.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer, component;
			double bkg, amp, mean, sigma, integral, flag;
			sector = layer = component = 
				bkg = amp = mean = sigma = integral = flag = 0;
			ss >> sector >> layer >> component >>
				bkg >> amp >> mean >> sigma >> integral >> flag;

			int BAR_ID = sector*100+layer*10+component;

			FTDCGlobal[BAR_ID] = mean;
			FTDCToFRes[BAR_ID] = sigma;
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

void BANDReco::readStatus(){
	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	std::string line;
	std::ifstream f;

	// Load the Edep calibration:
	f.open(path+"/status_list.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer, component;
			int status;
			sector = layer = component =
				status = 0;
			ss >> sector >> layer >> component >> status;

			int BARID = sector*100 + layer*10 + component;
			STATUS[BARID] = status;
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
	if( !loaded_TW ) return 0.;
	double A = *x;
	double f0 = p[0] + p[1]/pow(A,p[2]); 				// p[3]-p[5] not used
	double f1 = p[6]/A + p[7]/(exp((A-p[8])/p[9]) + p[10] );	// p[11] not used
	double f2 = p[12] + p[13]*A;					// p[14]-p[17] not used
	double f3 = p[18]*sin((A-p[19])/p[20]) + p[21]/pow(A,p[22]);	// p[23] not used
	if( p[20] == 0 ) f3 = 0;

	if( f0 != f0 || f1 != f1 || f2 != f2 || f3 != f3 ){ cerr << "issue with timewalk\n"; exit(-1); }

	return f0+f1+f2+f3;
}

bool BANDReco::check_bar(const Bar * this_bar){
	PMT left 	= this_bar->left;
	PMT right 	= this_bar->right;
	if( left.amp == 0 || right.amp == 0 		) return false;
	if( left.adc == 0 || right.adc == 0		) return false;
	if( left.ftdc == 0 || right.ftdc == 0		) return false;
	if( left.tdc_corr == 0 || right.tdc_corr == 0 	) return false;
	if( left.ftdc_corr == 0 || right.ftdc_corr == 0 ) return false;
	if( left.tdc == 0 || right.tdc == 0 		) return false;
	if( left.sector != right.sector			) return false;
	if( left.layer != right.layer			) return false;
	if( left.component != right.component		) return false;
	if( left.layer == 6 || right.layer == 6 	) return true;	// if it's veto, we don't want to check the PMT difference
	if( abs(left.PMT_ID-right.PMT_ID)>1		) return false;

	return true;
}
bool BANDReco::check_status( const Bar * this_bar){
	if( STATUS[this_bar->Bar_ID] == 0 ) return false;
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
		int sector 	= (int) band_adc->getInt( 0, row );
		int layer 	= (int) band_adc->getInt( 1, row );
		int component 	= (int) band_adc->getInt( 2, row );
		int order	= (int) band_adc->getInt( 3, row );
		int PMT_ID	= sector*1000 + layer*100 + component*10 + order;
		// filter out the photodiode in the data -- TODO FIX THIS FOR CALIBRATION
		if( PMT_ID == 6661 || PMT_ID == 6660 ) continue;

		int adc		= (int) band_adc->getInt( 4, row );
		int amp		= (int) band_adc->getInt( 5, row );
		double ftdc	= (double) band_adc->getFloat( 6, row );
		int ped		= (int) band_adc->getInt( 7, row );
		
		if( amp >= 4095 ) continue;
		amp -= ped; // subtract off pedestal since it's not in the ExtendedFADCFitter
		if( loaded_TW && amp < 250 ) continue;		// cut for TW fits

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
		int sector 	= (int) band_tdc->getInt( 0 , row );
		int layer 	= (int) band_tdc->getInt( 1 , row );
		int component 	= (int) band_tdc->getInt( 2 , row );
		int order	= (int) band_tdc->getInt( 3 , row ) - 2; // minus two because TDC order is 2,3 for PMTs L,R
		int PMT_ID	= sector*1000 + layer*100 + component*10 + order;

		int tdc_raw	= (int) band_tdc->getInt( 4 , row );
		double tdc	= tdc_raw * 0.02345 - phaseCorr;

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
			candidate.PMT_ID	= this_pmt.PMT_ID;
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
			if( check_status(&newbar) == false ) continue; // check the status table

			// Now we can check to make sure we have a valid bar and set the quantities correctly
			double tdiff		= (newbar.left.tdc_corr - newbar.right.tdc_corr) - TDCOffsets[Bar_ID];
			double tdiff_ftdc	= (newbar.left.ftdc_corr - newbar.right.ftdc_corr) - FTDCOffsets[Bar_ID];
			// no velocity max Tdiff check since we don't have a full bar
	
			// Create the ToF:
			double meantime		= 0.5*(newbar.left.tdc_corr + newbar.right.tdc_corr - TDCOffsets[Bar_ID]	) - TDCPaddle[Bar_ID] - TDCLayer[Bar_ID] - TDCGlobal[Bar_ID];
			double meantime_ftdc	= 0.5*(newbar.left.ftdc_corr + newbar.right.ftdc_corr - FTDCOffsets[Bar_ID]	) - FTDCPaddle[Bar_ID] - FTDCLayer[Bar_ID] - FTDCGlobal[Bar_ID];

			// Create the x,y,z:
			double x 	= GlobalX[Bar_ID];
			double x_ftdc 	= GlobalX[Bar_ID];
			double y 	= GlobalY[Bar_ID];
			double z 	= GlobalZ[Bar_ID] + BAND_MOTHER_OFFSET;
										// GlobalZ is w.r. to the ideal target position
										// SVT cart is shifted w.r. to the ideal target position

			// Create the Edep:
			double sec_length = bar_lengths[newbar.sector-1];
			double mu = 1e10;
			double adcL_corr = newbar.left.adc 	* exp( (sec_length/2.-x)/mu );
			double adcR_corr = newbar.right.adc 	* exp( (sec_length/2.-x)/mu );
			double ampL_corr = newbar.left.amp 	* exp( (sec_length/2.-x)/mu );
			double ampR_corr = newbar.right.amp 	* exp( (sec_length/2.-x)/mu );

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
			newbar.XFtdc		= x_ftdc;
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
			else if( candidate1.order == 1 ){ // candidate 1 is right, candidate 2 is left
				newbar.left 	= candidate2;
				newbar.right 	= candidate1;
			}
			else{ cerr << "what\n"; continue;}
			newbar.sector 		= sector;
			newbar.layer		= layer;
			newbar.component	= component;
			newbar.Bar_ID		= Bar_ID;
			if( check_bar(&newbar) == false ) continue;	// make sure all the values for left and right are set
			if( check_status(&newbar) == false ) continue; // check the status table

			// Now we can check to make sure we have a valid bar and set the quantities correctly
			double tdiff		= (newbar.left.tdc_corr - newbar.right.tdc_corr) - TDCOffsets[Bar_ID];
			double tdiff_ftdc	= (newbar.left.ftdc_corr - newbar.right.ftdc_corr) - FTDCOffsets[Bar_ID];
			double maxTdiff		= bar_lengths[sector-1] / TDCVelocity[Bar_ID];
			double maxTdiff_ftdc	= bar_lengths[sector-1] / FTDCVelocity[Bar_ID];
			if( maxTdiff == maxTdiff || maxTdiff_ftdc == maxTdiff_ftdc ){
				if( fabs(tdiff) > maxTdiff 		) continue;
				if( fabs(tdiff_ftdc) > maxTdiff_ftdc 	) continue;
			}
	
			// Create the ToF:
			double meantime		= 0.5*(newbar.left.tdc_corr + newbar.right.tdc_corr - TDCOffsets[Bar_ID]	) - TDCPaddle[Bar_ID] - TDCLayer[Bar_ID] - TDCGlobal[Bar_ID];
			double meantime_ftdc	= 0.5*(newbar.left.ftdc_corr + newbar.right.ftdc_corr - FTDCOffsets[Bar_ID]	) - FTDCPaddle[Bar_ID] - FTDCLayer[Bar_ID] - FTDCGlobal[Bar_ID];

			// Create the x,y,z:
			double x 	= GlobalX[Bar_ID] - (0.5 * tdiff * TDCVelocity[Bar_ID]);
			double x_ftdc 	= GlobalX[Bar_ID] - (0.5 * tdiff_ftdc * FTDCVelocity[Bar_ID]);
			double y = GlobalY[Bar_ID];
			double z = GlobalZ[Bar_ID] + BAND_MOTHER_OFFSET;
										// GlobalZ is w.r. to the ideal target position
										// SVT cart is shifted w.r. to the ideal target position

			// Create the Edep:
			double sec_length = bar_lengths[newbar.sector-1];
			double mu = 1e10;
			double adcL_corr = newbar.left.adc 	* exp( (sec_length/2.-x)/mu );
			double adcR_corr = newbar.right.adc 	* exp( (sec_length/2.-x)/mu );
			double ampL_corr = newbar.left.amp 	* exp( (sec_length/2.-x)/mu );
			double ampR_corr = newbar.right.amp 	* exp( (sec_length/2.-x)/mu );

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
			newbar.XFtdc		= x_ftdc;
			newbar.Y		= y;
			newbar.Z		= z;
			candidate_bars[Bar_ID] = newbar;
		}
	}

	return;
}

void BANDReco::storeHits( int& mult , bandhit * hits , const double starttime , const double vtx_z ){

	mult = 0;
	map<int,Bar>::const_iterator bar_it;
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
		hits[mult].setXFadc		(this_bar.XFtdc);
		hits[mult].setY			(this_bar.Y);
		hits[mult].setZ			(this_bar.Z - vtx_z );

		hits[mult].setRawLtdccorr	(this_bar.left.tdc_corr);					//TW corrected time
		hits[mult].setRawLtdc		(this_bar.left.tdc);						// 0.02345 conversion and phase corr
		hits[mult].setPmtLtdc		( (this_bar.left.tdc + this_bar.left.trigphase)/0.02345 );	// raw TDC channel

		hits[mult].setRawRtdccorr	(this_bar.right.tdc_corr);					//TW corrected time
		hits[mult].setRawRtdc		(this_bar.right.tdc);						// 0.02345 conversion and phase corr
		hits[mult].setPmtRtdc		( (this_bar.right.tdc + this_bar.right.trigphase)/0.02345 );	// raw TDC channel
		
		hits[mult].setRawLtfadc		(this_bar.left.ftdc_corr);	// any additional correction to FTDC but not used at the moment
		hits[mult].setRawRtfadc		(this_bar.right.ftdc_corr);
		hits[mult].setPmtLtfadc		(this_bar.left.ftdc);		// copy of above since there are no additional corrections
		hits[mult].setPmtRtfadc		(this_bar.right.ftdc);

		hits[mult].setRawLadc		(this_bar.AdcL);		// attenuation-corrected ADC
		hits[mult].setRawRadc		(this_bar.AdcR);
		hits[mult].setRawLamp		(this_bar.AmpL);		// attenuation-corrected AMP
		hits[mult].setRawRamp		(this_bar.AmpR);

		hits[mult].setPmtLamp		(this_bar.left.amp);
		hits[mult].setPmtRamp		(this_bar.right.amp);
		hits[mult].setPmtLadc		(this_bar.left.adc);
		hits[mult].setPmtRadc		(this_bar.right.adc);
		hits[mult].setPmtLped		(this_bar.left.ped);
		hits[mult].setPmtRped		(this_bar.right.ped);
	
		hits[mult].setStatus		(0.);
		hits[mult].setDL		( TVector3(hits[mult].getX(),hits[mult].getY(),hits[mult].getZ()) );

		mult++;
	}
	return;
}
