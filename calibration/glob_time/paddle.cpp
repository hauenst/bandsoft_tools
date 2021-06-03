#include <iostream>
#include <fstream>
#include <sstream>

#include "reader.h"
#include "bank.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TLine.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

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
	double ftdc_corr= 0;
	int ped 	= 0;
	double tdc	= 0;
	double tdc_corr	= 0;
	double trigphase= 0;
	int idx_fadc	= -1;
	int idx_tdc	= -1;
};


struct Bar{
	PMT left;
	PMT right;
	int sector	= 0;
	int layer	= 0;
	int component	= 0;
	int Bar_ID 	= 0;
	double Edep	= 0;
	double Tof	= 0;
	double TofFadc	= 0;
	double X	= 0;
	double Y	= 0;
	double Z	= 0;
	double Tdiff		= 0;
	double TdiffFadc	= 0;
};

map<int,vector<double>> TWParamsAMP;
map<int,double> TDCOffsets;
map<int,double> FADCOffsets;
map<int,double> TDCVelocity;
map<int,double> FADCVelocity;

double timewalk( double *x , double *p);
double getTriggerPhase( long timeStamp ) ;
void readParameters(void);
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
double bar_lengths[5] = {163.7,201.9,51.2,51.2,201.9};
void offsetFit( TH1D * hist , TCanvas * c , int cd, int layer, double &amp, double &mean, double &sigma );
bool check_bar(Bar * this_bar);

int period = -1;
bool SPRING2019 = false;
bool FALL2019_WINTER2020 = false;

int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [calibration period] [inputFile1] [inputFile2] ... \n";
		cerr << "\t\t[calibration peroid] -- 0 = SPRING2019 / 1 == FALL2019_WINTER2020\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}
	

	period = atoi(argv[1]);
	(period == 0) ? SPRING2019 = true : FALL2019_WINTER2020 = true;

	// Read in the TW parameters and LR offsets and Eff Vel
	readParameters();

	// Create histograms for storage
	TH1D h1_poff[5][6][7];
	TH1D h1_poff_fadc[5][6][7];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				h1_poff[sector][layer][comp]
						= TH1D(	Form("poff_%i_%i_%i",sector+1,layer+1,comp+1),
							Form("poff_%i_%i_%i",sector+1,layer+1,comp+1),8000,-200,200);
				h1_poff_fadc[sector][layer][comp]
						= TH1D(	Form("poff_fadc_%i_%i_%i",sector+1,layer+1,comp+1),
							Form("poff_fadc_%i_%i_%i",sector+1,layer+1,comp+1),8000,-200,200);
			}
		}
	}


	// Setup hipo reading for this file
	for( int fi = 2 ; fi < argc ; ++fi ){
		TString inputFile = argv[fi];
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
			//if( event_counter > 500000 ) break;
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
			map<int,vector<PMT>> raw_pmts_adc;
			map<int,vector<PMT>> raw_pmts_tdc;
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
				this_pmt.idx_fadc	= row;
				raw_pmts_adc[PMT_ID].push_back(this_pmt);
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
			
			// Now let's loop over all the FADC PMT hits to do an additional filter.
			// and then add a valid PMT to our bar list
			map<int,vector<PMT>>::iterator it;
			map<int,PMT> candidate_pmts;
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
					candidate.idx_fadc	= this_pmt.idx_fadc;
					candidate.tdc		= raw_pmts_tdc[PMT_ID][0].tdc;
					candidate.trigphase	= raw_pmts_tdc[PMT_ID][0].trigphase;
					candidate.idx_tdc	= raw_pmts_tdc[PMT_ID][0].idx_tdc;

					double tdc 		= candidate.tdc;
					double * tw_pars 	= &TWParamsAMP[PMT_ID][0];
					double amp		= candidate.amp;
					double tdc_ampcorr 	= tdc - timewalk(&amp,tw_pars);
					candidate.tdc_corr	= tdc_ampcorr;				// store TW corrected time
					candidate.ftdc_corr	= candidate.ftdc;			// no TW correction for FADC

					candidate_pmts[PMT_ID] = candidate;
				}
			}

			// Now that we have a single struct per PMT, let's loop over that to create our bars
			map<int,PMT>::iterator pmt_i;
			map<int,PMT>::iterator pmt_j;
			map<int,Bar> candidate_bars;
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
					double tdiff_ftdc	= (newbar.left.ftdc_corr - newbar.right.ftdc_corr) - FADCOffsets[Bar_ID];
					double meantime		= 0.5*(newbar.left.tdc_corr + newbar.right.tdc_corr - TDCOffsets[Bar_ID]	);
					double meantime_ftdc	= 0.5*(newbar.left.ftdc_corr + newbar.right.ftdc_corr - FADCOffsets[Bar_ID]	);
					newbar.Tof 		= meantime;
					newbar.TofFadc 		= meantime_ftdc;
					newbar.Tdiff		= tdiff;
					newbar.TdiffFadc	= tdiff_ftdc;
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
					double tdiff_ftdc	= (newbar.left.ftdc_corr - newbar.right.ftdc_corr) - FADCOffsets[Bar_ID];
					double maxTdiff		= bar_lengths[sector-1] / TDCVelocity[Bar_ID];
					double maxTdiff_ftdc	= bar_lengths[sector-1] / FADCVelocity[Bar_ID];
					if( fabs(tdiff) > maxTdiff 		) continue;
					if( fabs(tdiff_ftdc) > maxTdiff_ftdc 	) continue;

					double meantime		= 0.5*(newbar.left.tdc_corr + newbar.right.tdc_corr - TDCOffsets[Bar_ID]	);
					double meantime_ftdc	= 0.5*(newbar.left.ftdc_corr + newbar.right.ftdc_corr - FADCOffsets[Bar_ID]	);
					newbar.Tof 		= meantime;
					newbar.TofFadc 		= meantime_ftdc;
					newbar.Tdiff		= tdiff;
					newbar.TdiffFadc	= tdiff_ftdc;
					candidate_bars[Bar_ID] = newbar;
				}
			}

			// Now let's loop over all of the bars to fill our histograms!
			map<int,Bar>::iterator bar_it;
			for( bar_it = candidate_bars.begin() ; bar_it != candidate_bars.end() ; ++bar_it){
				int Bar_ID 	= bar_it->first;
				Bar this_bar 	= bar_it->second;

				// Grab the reference and make sure it exists:
				int this_ref 	= 200 + this_bar.layer*10 + 1; // sector 2, layer X, component 1
				if( check_bar(&candidate_bars[this_ref]) == false ) continue;
				double tdc_ref 	= candidate_bars[this_ref].Tof;
				double ftdc_ref = candidate_bars[this_ref].TofFadc;

				int sector 	= this_bar.sector;
				int layer 	= this_bar.layer;
				int component 	= this_bar.component;
				double mean_tdc 	= this_bar.Tof;
				double mean_ftdc 	= this_bar.TofFadc;

				h1_poff[sector-1][layer-1][component-1]		.Fill( mean_tdc - tdc_ref 	);
				h1_poff_fadc[sector-1][layer-1][component-1]	.Fill( mean_ftdc - ftdc_ref 	);
				
			}
		} // end loop over events
	} // end loop over files

	// Create output canvas to draw on
	TFile * outFile = new TFile("off_paddle.root","RECREATE");
	TCanvas c0("c0","c0",700,900);
	c0.Print("off_paddle.pdf(");
	ofstream outfile;
	outfile.open("paddle_offsets.txt");
	std::setprecision(9);
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			TCanvas cSLC(Form("sector_%i_layer_%i",(sector+1),(layer+1)),Form("sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			cSLC.Divide(2,7);
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				// set titles:
				cout << (sector+1) << " " << (layer+1) << " " << (comp+1)  << "\n";
				h1_poff[sector][layer][comp].SetTitle(Form("TDC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_poff_fadc[sector][layer][comp].SetTitle(Form("FADC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_poff[sector][layer][comp].Write();
				h1_poff_fadc[sector][layer][comp].Write();

				if( h1_poff[sector][layer][comp].Integral() == 0 || 
					h1_poff_fadc[sector][layer][comp].Integral() == 0 ||
					(sector+1 == 2 && comp+1 == 1 ) ){
						outfile << (sector+1) << "\t" << (layer+1) << "\t" << (comp+1) << "\t"
							<< 0 << "\t" << 0 << "\t" << 0 << "\t"
							<< 0 << "\t" << 0 << "\t" << 0 << "\n";
						continue;
				}

				double amp, mean, sigma;
				amp = mean = sigma = 0.;
				offsetFit( &h1_poff[sector][layer][comp] 	, &cSLC , 2*comp+1 , (layer+1) , amp, mean, sigma );
				outfile << (sector+1) << "\t" << (layer+1) << "\t" << (comp+1) << "\t"
					<< amp << "\t" << mean << "\t" << sigma << "\t";


				offsetFit( &h1_poff_fadc[sector][layer][comp]   , &cSLC , 2*comp+2 , (layer+1) , amp, mean, sigma );
				outfile << amp << "\t" << mean << "\t" << sigma << "\n";


				
			}
			cSLC.Print("off_paddle.pdf");
		}
	}
	outfile.close();
	c0.Print("off_paddle.pdf)");
	outFile->Close();
	delete outFile;
	

	return 0;
}

bool check_bar(Bar * this_bar){
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

	// Load L-R offset file:
	f.open(path+"/lr_offsets.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int sector, layer ,component;
			double tdc_off, tdc_veff, tdc_width;
			double fadc_off, fadc_veff, fadc_width;
			tdc_off = tdc_veff = tdc_width = 
			fadc_off = fadc_veff = fadc_width = 0;
			ss >> sector >> layer >> component >>
				tdc_off >> tdc_veff >> tdc_width >>
				fadc_off >> fadc_veff >> fadc_width;

			int BARID = sector*100 + layer*10 + component;
			TDCOffsets[BARID] = tdc_off;
			TDCVelocity[BARID] = tdc_veff;
			FADCOffsets[BARID] = fadc_off;
			FADCVelocity[BARID] = fadc_veff;
		}
	}
	f.close();

	return;
}


double timewalk( double *x , double *p){
	double A = *x;
	double f0 = p[0] + p[1]/pow(A,p[2]); 				// p[3]-p[5] not used
	double f1 = p[6]/A + p[7]/(exp((A-p[8])/p[9]) + p[10] );	// p[11] not used
	double f2 = p[12] + p[13]*A;					// p[14]-p[17] not used
	double f3 = p[18]*sin((A-p[19])/p[20]) + p[21]/pow(A,p[22]);	// p[23] not used
	if( p[20] == 0 ) f3 = 0;

	if( f0 != f0 || f1 != f1 || f2 != f2 || f3 != f3 ){ cerr << "issue with timewalk\n"; exit(-1); }

	return f0+f1+f2+f3;
}

void offsetFit( TH1D * hist , TCanvas * c , int cd, int layer, double &amp, double &mean, double &sigma ){

	c->cd(cd);
	
	// get the initial mean/sigma and zoom in around the peak:
	mean = hist->GetMean();
	sigma = hist->GetStdDev();
	cout << "\tinit: " << mean << " " << sigma << " " << mean-5*sigma << " " << mean+5*sigma << "\n";
	hist->GetXaxis()->SetRangeUser( mean - 5*sigma , mean + 5*sigma );

	// update the mean/sigma based on zoomed-in histogram and then get a min/max x for fitting:
	mean = hist->GetMean();
	sigma = hist->GetStdDev();
	double xpeak = hist->GetXaxis()->GetBinCenter( hist->GetMaximumBin() );
	double xmin = xpeak - 3.00*hist->GetStdDev();
	double xmax = xpeak + 0.75*hist->GetStdDev();
	if( layer == 6 ){
		xmin = xpeak - 3.00*hist->GetStdDev();
		xmax = xpeak + 3.00*hist->GetStdDev();
	}
	cout << "\tfit: " << xmin << " " << xmax << "\n";
	TFitResultPtr ptr = hist->Fit("gaus","QESR","",xmin,xmax);

	amp 	= ptr->Parameter(0);
	mean 	= ptr->Parameter(1);
	sigma 	= ptr->Parameter(2);

	// if bad fit:
	if( sigma > 1 ){
		amp = hist->GetMaximum();
		mean = hist->GetMean();
		sigma = hist->GetStdDev();
	}

	// Do final zooming based on refined mean/sigma:
	hist->GetXaxis()->SetRangeUser( mean - 10*sigma , mean + 10*sigma );

	hist->Draw();
	c->Modified(); c->Update();
	return;
}
