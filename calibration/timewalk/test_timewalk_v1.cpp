#include <iostream>
#include <fstream>
#include <sstream>

#include "reader.h"
#include "bank.h"

#include "TString.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"

using namespace std;

struct PMT{
	int PMT_ID 	= 0;
	int sector	= 0;
	int layer	= 0;
	int component	= 0;
	int order	= 0;
	int adc		= 0;
	int amp		= 0;
	double ftdc	= 0;
	int ped 	= 0;
	double tdc	= 0;

};
map<int,vector<double>> TWParamsAMP;
map<int,vector<double>> TWParamsADC;
map<int,vector<double>> TWParamsAMP_It0;
map<int,vector<double>> TWParamsADC_It0;
map<int,vector<double>> TWParamsAMP_It2;
map<int,vector<double>> TWParamsADC_It2;

double getTriggerPhase( long timeStamp ) ;
void readParameters(void);
double getTriggerPhase( long timeStamp ) ;
void fitTW(TH2D * hist , TCanvas * c, int cd, int s, int l, int co, int o, double cut, 
		double &par1, double &par2 , double &par3, double &par4, double &par5, double &par6,
		double &par1_err, double &par2_err , double &par3_err, double &par4_err, double &par5_err, double &par6_err);
void walkCorr(	vector<double> *adcs		,
		vector<double> *adcsErr		,
		vector<double> *times		,
		vector<double> *timesErr	,
		vector<double> *res		,
		vector<double> *resErr		,
		double widthCut			,
		TH2D * hist			);
int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE, int thres, int lastBin );
double wlk( double *x , double *p);

int iteration = -1;
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [calibration period] [iteration] [inputFile] \n";
		cerr << "\t\t[calibration peroid] -- 0 = SPRING2019 / 1 == FALL2019_WINTER2020\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}
	
	// Read in the TW parameters
	readParameters();

	int period = atoi(argv[1]);
	bool SPRING2019 = false, FALL2019_WINTER2020 = false;
	(period == 0) ? SPRING2019 = true : FALL2019_WINTER2020 = true;

	int REFID = -1;
	if( SPRING2019 ){
		REFID = 3660; // V16-A
	}
	else if( FALL2019_WINTER2020 ){
		REFID = 4161; // 116B-R
	}

	iteration = atoi(argv[2]);


	// Create histograms for storage
	TH2D h2_tadc[5][6][7][2];
	TH2D h2_tamp[5][6][7][2];
	TH2D h2_tadc_overfill[5][6][7][2];
	TH2D h2_tamp_overfill[5][6][7][2];
	TH1D h1_tdiff[5][6][7][2];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				for( int order = 0 ; order < 2 ; ++order ){
					h1_tdiff[sector][layer][comp][order]
							= TH1D(	Form("tdiff_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("tdiff_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),200,-2,2);
					h2_tadc[sector][layer][comp][order] 
							= TH2D(	Form("adc_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("adc_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),300,0,30000,80,-2,2);
					h2_tadc_overfill[sector][layer][comp][order] 
							= TH2D(Form("adcO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("adcO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),300,15000,45000,80,-2,2);
					h2_tamp[sector][layer][comp][order] 
							= TH2D(Form("amp_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("amp_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),250,0,5000,40,-2,2); // first iteration
								//Form("amp_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),500,0,5000,80,-2,2); // second iteration
					h2_tamp_overfill[sector][layer][comp][order] 
							= TH2D(Form("ampO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("ampO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),250,0,5000,80,-2,2); // first iteration
								//Form("ampO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),500,0,5000,80,-2,2); // second iteration
				}
			}
		}
	}
	


	// Setup hipo reading for this file
	TString inputFile = argv[3];
	hipo::reader reader;
	reader.open(inputFile);
	hipo::dictionary  factory;
	hipo::schema	  schema;
	reader.readDictionary(factory);
	hipo::bank	band_adc		(factory.getSchema("BAND::adc"		));
	hipo::bank	band_tdc		(factory.getSchema("BAND::tdc"		));
	hipo::bank  	run_config 		(factory.getSchema("RUN::config"	));
	hipo::event 	readevent;

	int event_counter = 0;
	// Loop over all events in file
	while(reader.next()==true){


		// Count events
		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		//if( event_counter > 1000000 ) break;
		event_counter++;

		// Load data structure for this event:
		reader.read(readevent);
		readevent.getStructure(band_adc);
		readevent.getStructure(band_tdc);
		readevent.getStructure(run_config);

		//Get Event number from RUN::config
		int eventnumber = run_config.getInt( 1 , 0 );
		long timestamp = run_config.getLong(4 , 0 );
		double phaseCorr = getTriggerPhase(timestamp);

		// Loop over all the ADC hits:
		map<int,vector<PMT>> pmts_adc;
		map<int,vector<PMT>> pmts_tdc;
		for( int row = 0 ; row < band_adc.getRows() ; ++row ){
			int sector 	= band_adc.getInt( 0, row );
			int layer 	= band_adc.getInt( 1, row );
			int component 	= band_adc.getInt( 2, row );
			int order	= band_adc.getInt( 3, row );
			int PMT_ID	= sector*1000 + layer*100 + component*10 + order;

			int adc		= band_adc.getInt( 4, row );
			int amp		= band_adc.getInt( 5, row );
			double ftdc	= band_adc.getFloat( 6, row );
			int ped		= band_adc.getInt( 7, row );

			//if( amp < 2000 && PMT_ID!=REFID ) continue; // second iteration

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
			pmts_adc[PMT_ID].push_back(this_pmt);
		}
		// Loop over all the TDC hits to find a match:
		for( int row = 0 ; row < band_tdc.getRows() ; ++row ){
			int sector 	= band_tdc.getInt( 0 , row );
			int layer 	= band_tdc.getInt( 1 , row );
			int component 	= band_tdc.getInt( 2 , row );
			int order	= band_tdc.getInt( 3 , row ) - 2;
			int PMT_ID	= sector*1000 + layer*100 + component*10 + order;

			double tdc	= band_tdc.getInt( 4 , row ) * 0.02345 - phaseCorr;

			// If we do not have a matching FADC hit for the PMT, we 
			// skip this event completely
			if( pmts_adc[PMT_ID].size() == 0 ) continue;

			PMT this_pmt;
			this_pmt.PMT_ID 	= PMT_ID;
			this_pmt.sector		= sector;
			this_pmt.layer		= layer;
			this_pmt.component	= component;
			this_pmt.order		= order;
			this_pmt.tdc		= tdc;
			pmts_tdc[PMT_ID].push_back(this_pmt);
		}	
		
		// First we need to check if our photodiode reference signal exists:
		if( pmts_adc.count(REFID) == 0 || pmts_tdc.count(REFID) == 0 ) continue;
			// Make sure for the reference, there is only 1 TDC
		if( pmts_tdc[REFID].size() != 1 ) continue;
		double TREF = pmts_tdc[REFID][0].tdc;

		// Now let's loop over all the FADC PMT hits to do an additional filter.
		map<int,vector<PMT>>::iterator it;
		for( it = pmts_adc.begin() ; it != pmts_adc.end() ; ++it){
			int PMT_ID = it->first;
			if( PMT_ID == REFID ) continue;

			// If we don't have a TDC, we can't do anything
			if( pmts_tdc[PMT_ID].size() == 0 ) continue;

			// If we have more than one ADC, let's not consider this
			if( it->second.size() > 1 ){
				cerr << "unexpected behavior with double adc hit for single pmt. exiting...\n";
				exit(-1);
			}
			
			// If we have more than one TDC, let's not deal with this for the moment
			if( pmts_tdc[PMT_ID].size() > 1 ){
				continue;
			}

			// If we only have 1 of each:
			if( it->second.size() == 1 && pmts_tdc[PMT_ID].size() == 1){
				PMT this_pmt 	= it->second[0];
				int sector 	= this_pmt.sector;
				int layer	= this_pmt.layer;
				int component	= this_pmt.component;
				int order	= this_pmt.order;
				int adc		= this_pmt.adc;
				int amp		= this_pmt.amp;

				double tdc = pmts_tdc[PMT_ID][0].tdc;

				double tdc_ampcorr;
				double tdc_adccorr;
				if( iteration >= 0 ){
					tdc_ampcorr = tdc 
						- ( TWParamsAMP[PMT_ID][0] + TWParamsAMP[PMT_ID][1] / pow(amp,TWParamsAMP[PMT_ID][2]) ) ;
					tdc_adccorr = tdc
						- ( TWParamsADC[PMT_ID][0] + TWParamsADC[PMT_ID][1] / pow(adc,TWParamsADC[PMT_ID][2]) ) ;
				}
				if( iteration >= 1 ){
					tdc_ampcorr -= ( TWParamsAMP_It0[PMT_ID][0]/amp + TWParamsAMP_It0[PMT_ID][1]/(TWParamsAMP_It0[PMT_ID][4] + exp((amp-TWParamsAMP_It0[PMT_ID][2])/TWParamsAMP_It0[PMT_ID][3]) ) );
					//tdc_ampcorr -= ( TWParamsAMP_It0[PMT_ID][0] + TWParamsAMP_It0[PMT_ID][1] / pow(amp,TWParamsAMP_It0[PMT_ID][2]) );  // first iteration
					//tdc_adccorr -= ( TWParamsADC_It0[PMT_ID][0] + TWParamsADC_It0[PMT_ID][1] / pow(adc,TWParamsADC_It0[PMT_ID][2]) );  // first iteration
					
				}
				if( iteration >= 2 && amp > 2000 ){
					tdc_ampcorr -= ( TWParamsAMP_It2[PMT_ID][0] + TWParamsAMP_It2[PMT_ID][1] * amp ) ; // second iteration
					tdc_adccorr -= ( TWParamsADC_It2[PMT_ID][0] + TWParamsADC_It2[PMT_ID][1] * adc ) ; // second iteration
				}
	

				h1_tdiff[sector-1][layer-1][component-1][order].Fill( tdc_ampcorr - tdc_adccorr );

				if( amp >= 4095 ){
					h2_tadc_overfill[sector-1][layer-1][component-1][order].Fill( adc , tdc_ampcorr - TREF );
					h2_tamp_overfill[sector-1][layer-1][component-1][order].Fill( amp , tdc_ampcorr - TREF );
				}
				else{
					h2_tadc[sector-1][layer-1][component-1][order].Fill( adc , tdc_ampcorr - TREF );
					h2_tamp[sector-1][layer-1][component-1][order].Fill( amp , tdc_ampcorr - TREF );
				}
			}
			
		}
	} // end loop over events

	// Create output canvas to draw on
	TFile * outFile = new TFile(Form("iter%i_timewalk.root",iteration),"RECREATE");
	TCanvas c0("c0","c0",700,900);
	TCanvas c1("c1","c1",700,900);
	c0.Print(Form("iter%i_adc_timewalk.pdf(",iteration));
	c1.Print(Form("iter%i_amp_timewalk.pdf(",iteration));
	ofstream outfile_adc[2];
	ofstream outfile_adcR[2];
	outfile_adc[0].open(Form("iter%i_paramsadcL.txt",iteration));
	outfile_adc[1].open(Form("iter%i_paramsadcR.txt",iteration));
	ofstream outfile_amp[2];
	outfile_amp[0].open(Form("iter%i_paramsampL.txt",iteration));
	outfile_amp[1].open(Form("iter%i_paramsampR.txt",iteration));
	std::setprecision(9);
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			TCanvas cSLC_adc(Form("adc_sector_%i_layer_%i",(sector+1),(layer+1)),Form("adc_sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			TCanvas cSLC_amp(Form("amp_sector_%i_layer_%i",(sector+1),(layer+1)),Form("amp_sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			cSLC_adc.Divide(4,7);
			cSLC_amp.Divide(4,7);
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				for( int order = 0 ; order < 2 ; ++order ){
					// set titles:
					cout << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << (order) << "\n";
					h2_tadc[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tadc[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h2_tadc_overfill[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tadc_overfill[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp_overfill[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp_overfill[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h1_tdiff[sector][layer][comp][order].Write();
					if( h2_tadc[sector][layer][comp][order].Integral() == 0 ||
						h2_tamp[sector][layer][comp][order].Integral() == 0 ){

							outfile_adc[order] << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << 
								0 << " " << 0 << " " << 0 << " " << 
								0 << " " << 0 << " " << 0 << " " << 
								0 << " " << 0 << " " << 0 << " " << 
								0 << " " << 0 << " " << 0 << "\n";

							outfile_amp[order] << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << 
								0 << " " << 0 << " " << 0 << " " << 
								0 << " " << 0 << " " << 0 << " " << 
								0 << " " << 0 << " " << 0 << " " << 
								0 << " " << 0 << " " << 0 << "\n";

							continue;
					}
					double par1, par2, par3, par4, par5, par6;
					double par1e, par2e, par3e, par4e, par5e, par6e;
					
					// Go to ADC canvas and draw first
					par1 = 0, par2 = 0, par3 = 0, par4 = 0, par5 = 0, par6 = 0;
					par1e = 0, par2e = 0, par3e = 0, par4e = 0, par5e = 0, par6e = 0;
					cSLC_adc.cd(comp*4+1 + order*2);
					h2_tadc[sector][layer][comp][order].Draw("COL");
					h2_tadc[sector][layer][comp][order].Write();
					h2_tadc_overfill[sector][layer][comp][order].Write();
						// draw TW:
					//fitTW( &h2_tadc[sector][layer][comp][order] , &cSLC_adc,  4*comp+1 + (order)*2 + 1, 
					//		(sector+1), (layer+1), (comp+1), order, 25000,   // first iteration
					//	       	par1 , par2 , par3, par4, par5, par6,
					//		par1e, par2e , par3e, par4e, par5e, par6e );
					outfile_adc[order] << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << 
						par1 << " " << par2 << " " << par3 << " " << par4 << " " << par5 << " " << par6 << " " << 
						par1e << " " << par2e << " " << par3e << " " << par4e << " " << par5e << " " << par6e << "\n";


					// Go to AMP canvas and draw first
					par1 = 0, par2 = 0, par3 = 0, par4 = 0, par5 = 0, par6 = 0;
					par1e = 0, par2e = 0, par3e = 0, par4e = 0, par5e = 0, par6e = 0;
					cSLC_amp.cd(comp*4+1 + order*2);
					h2_tamp[sector][layer][comp][order].Draw("COL");
					h2_tamp[sector][layer][comp][order].Write();
					h2_tamp_overfill[sector][layer][comp][order].Write();
						// draw TW:
					fitTW( &h2_tamp[sector][layer][comp][order] , &cSLC_amp,  4*comp+1 + (order)*2 + 1, 
							(sector+1), (layer+1), (comp+1), order, 4095,   // first iteration
						       	par1 , par2 , par3, par4, par5, par6,
							par1e, par2e , par3e, par4e, par5e, par6e );
					outfile_amp[order] << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << 
						par1 << " " << par2 << " " << par3 << " " << par4 << " " << par5 << " " << par6 << " " << 
						par1e << " " << par2e << " " << par3e << " " << par4e << " " << par5e << " " << par6e << "\n";

				}
			}
			cSLC_adc.Print(Form("iter%i_adc_timewalk.pdf",iteration));
			cSLC_amp.Print(Form("iter%i_amp_timewalk.pdf",iteration));
		}
	}
	c0.Print(Form("iter%i_adc_timewalk.pdf)",iteration));
	c1.Print(Form("iter%i_amp_timewalk.pdf)",iteration));
	outfile_adc[0].close();
	outfile_adc[1].close();
	outfile_amp[0].close();
	outfile_amp[1].close();
	outFile->Close();
	delete outFile;


	return 0;
}

double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}

void readParameters(void){

	std::string line;
	std::ifstream f;

	/*
	f.open("../../include/ccdb_tables/time_walk_amp_left.txt");
	f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, par_c1, par_a2, par_b2, par_c2 };
			TWParams[PMTID] = temp;
		}
	}
	f.close();

	f.open("../../include/ccdb_tables/time_walk_amp_right.txt");
	f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, par_c1, par_a2, par_b2, par_c2 };
			TWParams[PMTID] = temp;
		}
	}
	f.close();
	*/

	f.open("init_paramsampL.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;

			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, par_c1, 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsAMP[PMTID] = temp;
		}
	}
	f.close();
	f.open("init_paramsampR.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, par_c1, 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsAMP[PMTID] = temp;
		}
	}
	f.close();

	f.open("init_paramsadcL.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, par_c1, 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsADC[PMTID] = temp;
		}
	}
	f.close();
	f.open("init_paramsadcR.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, par_c1, 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsADC[PMTID] = temp;
		}
	}
	f.close();

	f.open("iter0_paramsampL.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;

			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, par_c1, par_a2, par_b2, 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsAMP_It0[PMTID] = temp;
		}
	}
	f.close();
	f.open("iter0_paramsampR.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, par_c1, par_a2, par_b2, 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsAMP_It0[PMTID] = temp;
		}
	}
	f.close();

	f.open("iter0_paramsadcL.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, par_c1, par_a2, par_b2, 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsADC_It0[PMTID] = temp;
		}
	}
	f.close();
	f.open("iter0_paramsadcR.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, par_c1, par_a2, par_b2, 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsADC_It0[PMTID] = temp;
		}
	}
	f.close();

	/*
	f.open("iter1_paramsampL.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
		//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;

			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, par_c1, 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsAMP_It2[PMTID] = temp;
		}
	}
	f.close();
	f.open("iter1_paramsampR.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, 0. , 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsAMP_It2[PMTID] = temp;
		}
	}
	f.close();

	f.open("iter1_paramsadcL.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 0;
			vector<double> temp{ par_a1, par_b1, 0. , 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsADC_It2[PMTID] = temp;
		}
	}
	f.close();
	f.open("iter1_paramsadcR.txt");
	//f.ignore(1000, '\n'); // skip the first line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			//double par1, par2, par3, par4, par5, par6, par1_2, par2_2, par3_2, par4_2, par5_2, par6_2;
			//ss >> sector >> layer >> component >>
			//	par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> 
			//	par1_2 >> par2_2 >> par3_2 >> par4_2 >> par5_2 >> par6_2;
			double par_a1, par_b1, par_c1, par_a2, par_b2, par_c2;
			ss >> sector >> layer >> component >> par_a1 >> par_b1 >> par_c1 >> par_a2 >> par_b2 >> par_c2;
			
			int PMTID = sector*1000 + layer*100 + component*10 + 1;
			vector<double> temp{ par_a1, par_b1, 0. , 0., 0., 0. };
			//vector<double> temp{ par1, par2, par3, par4, par5, par6 };
			TWParamsADC_It2[PMTID] = temp;
		}
	}
	f.close();
	*/
	return;
}

void fitTW(TH2D * hist , TCanvas * c, int cd, int s, int l, int co, int o, double cut, 
		double &par1, double &par2 , double &par3, double &par4, double &par5, double &par6,
		double &par1_err, double &par2_err , double &par3_err, double &par4_err, double &par5_err, double &par6_err){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	walkCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, 0, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);
	for(int bin = 1 ; bin <= dim ; ++bin ){
		cout << xs[bin-1] << " " << ys[bin-1] << " " << xErrs[bin-1] << " " << yErrs[bin-1] << "\n";
	}
	//g->GetHistogram()->SetMaximum(0.1);	// third iteration (check results)
	//g->GetHistogram()->SetMinimum(-0.1);	// third iteration (check results)
	//int max = (int)cut;
	//g->GetXaxis()->SetLimits(0,max);	// third iteration (check results)
	TF1 * model = new TF1("timeWalk",wlk,0,20000,2);
	//model->SetParameter(0,0);
	//model->SetParameter(1,1E12);
	//model->SetParameter(2,5);
	//model->SetParLimits(0,0,4E4);
	//model->SetParameter(0,1000);
	//model->SetParameter(1,-11);
	//model->SetParameter(2,-16E3);
	//model->SetParameter(3,4E3);
	//model->SetParameter(4,-40);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par1 = ptr->Parameter(0);
	par1_err = ptr->ParError(0);

	par2 = ptr->Parameter(1);
	par2_err = ptr->ParError(1);

	par3 = ptr->Parameter(2);
	par3_err = ptr->ParError(2);

	par4 = ptr->Parameter(3);
	par4_err = ptr->ParError(3);

	par5 = ptr->Parameter(4);
	par5_err = ptr->ParError(4);

	par6 = ptr->Parameter(5);
	par6_err = ptr->ParError(5);

	c->cd(cd);
	//gStyle->SetTitleW(0.6);
	g->SetTitle(Form("SLCO: %i %i %i %i",s,l,co,o));
	g->SetMarkerStyle(20);
	g->Draw("AP");
		// third iteration (check results):
	//TLine *l1 = new TLine(0,0.04,max,0.04);
	//TLine *l2 = new TLine(0,-0.04,max,-0.04);
	//l1->SetLineWidth(2);
	//l1->SetLineStyle(9);
	//l2->SetLineWidth(2);
	//l2->SetLineStyle(9);
	//l1->Draw("SAME");
	//l2->Draw("SAME");


	c->Modified(); c->Update();
	return;
}

void walkCorr(	vector<double> *adcs		,
		vector<double> *adcsErr		,
		vector<double> *times		,
		vector<double> *timesErr	,
		vector<double> *res		,
		vector<double> *resErr		,
		double widthCut			,
		TH2D * hist			){
	int currBin = 0;
	int maxBin = hist->GetXaxis()->GetNbins();
	while( currBin < maxBin ){
		double xPt, yPt, yEr, ySig, ySigEr;
		double currBin_x = hist->GetXaxis()->GetBinCenter(currBin);
		int nEvents = 8000;
		//if( currBin_x > 600 ) nEvents = 4000;	// first iteration
		//if( currBin_x > 1000 ) nEvents = 2000;	// first iteration
		//if( currBin_x < 1000 ) nEvents = 250;	// first iteration
		//if( currBin_x < 1000 ) nEvents = 2000;	// first iteration
		//int nEvents = 3000;
		//if( currBin_x > 3000 ) nEvents = 500;	// second iteration
		//if( currBin_x < 500 ) nEvents = 500;	// second iteration
		//int nEvents = 4000; // third iteration (check results)
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr, nEvents, maxBin );
		currBin += step ;
		if( xPt < widthCut) continue;
		if( yEr > 0.1) continue;
		adcs		->push_back(xPt);
		adcsErr		->push_back(0);
		times		->push_back(yPt);
		timesErr	->push_back(yEr);
		res		->push_back(ySig);
		resErr		->push_back(ySigEr);
	}
}

double wlk( double *x , double *p){
	double var = *x;
	//return p[0] + p[1]*var + p[2]*pow(var,2) + p[3]*pow(var,3) + p[4]*pow(var,4) ;
	//return p[0] + p[1]*var + p[2]*pow(var,2) + p[3]*pow(var,3) + p[4]*pow(var,4) + p[5]*pow(var,5);
	return p[0] + p[1]*var;
	return p[0]/var + p[1]/( p[4]+exp((var-p[2])/p[3]) );
	return p[0] + p[1] / pow(var,p[2]); // first iteration
	//return p[0] + p[1]*var; // second iteration

}

int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE, int thres, int lastBin ){
	int cnt = 0;
	int step = 0;
	TCanvas * trash = new TCanvas("trash");
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= lastBin) break;
		step+=1;
	}

	if( cnt >= thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);

		// Getting the mean x value in this range:
		hist->GetXaxis()->SetRange(bin,bin+step);
		x = hist->GetMean(1);
		hist->GetXaxis()->SetRange();
		TFitResultPtr f = pj->Fit("gaus","QESR","",-2,2);
		y = f->Parameter(1);
		yE = f->ParError(1);
		sig = f->Parameter(2);
		sigE = f->ParError(2);

		if( write ){
			pj->Write();
		}
		delete pj;
	}

	delete trash;

	return step;
}
