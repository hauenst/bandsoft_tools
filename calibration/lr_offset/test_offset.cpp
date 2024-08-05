#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "reader.h"
#include "bank.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TLine.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

// TODO: we can check the slope on each side to make sure they match if I want to go higher

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
void intercept( TH1D * hist , double avg , int sector, double& lower , double& upper , double& al, double& au );
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
double bar_lengths[5] = {163.7,201.9,51.2,51.2,201.9};
void offsetFit( TH1D * hist , TCanvas * c , int cd, int sector, double &offset, double &veff, double &width );

int period = -1;
bool SPRING2019 = false;
bool FALL2019_WINTER2020 = false;

int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [calibration period] [inputFile] \n";
		cerr << "\t\t[calibration peroid] -- 0 = SPRING2019 / 1 == FALL2019_WINTER2020\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}
	

	period = atoi(argv[1]);
	(period == 0) ? SPRING2019 = true : FALL2019_WINTER2020 = true;

	// Read in the TW parameters
	readParameters();

	// Create histograms for storage
	TH1D h1_tdiff[5][6][7];
	TH1D h1_tdiff_fadc[5][6][7];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				h1_tdiff[sector][layer][comp]
						= TH1D(	Form("tdiff_%i_%i_%i",sector+1,layer+1,comp+1),
							Form("tdiff_%i_%i_%i",sector+1,layer+1,comp+1),1600,-40,40);
				h1_tdiff_fadc[sector][layer][comp]
						= TH1D(	Form("tdiff_fadc_%i_%i_%i",sector+1,layer+1,comp+1),
							Form("tdiff_fadc_%i_%i_%i",sector+1,layer+1,comp+1),1600,-40,40);
			}
		}
	}
	


	// Setup hipo reading for this file
	TString inputFile = argv[2];
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

			if( amp < 250 || amp >= 4095 ) continue; 

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
		
		// Now let's loop over all the FADC PMT hits to do an additional filter.
		// and then add a valid PMT to our bar list
		map<int,vector<PMT>>::iterator it;
		map<int,Bar> bars;
		for( it = pmts_adc.begin() ; it != pmts_adc.end() ; ++it){
			int PMT_ID = it->first;

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
				double amp	= this_pmt.amp;
				double tdc 	= pmts_tdc[PMT_ID][0].tdc;

				double * tw_pars 	= &TWParamsAMP[PMT_ID][0];
				double tdc_ampcorr 	= tdc - timewalk(&amp,tw_pars);

				// Since we have now a valid PMT with a TDC and FADC,
				// create that and store it for Bar pairing
				int Bar_ID = (PMT_ID)/10;
				PMT valid_pmt 		= this_pmt;
				valid_pmt.tdc 		= tdc;
				valid_pmt.tdc_corr 	= tdc_ampcorr;
				valid_pmt.ftdc_corr	= this_pmt.ftdc;
				

				if( order == 0 ){
					valid_pmt.tdc_corr	-= TDCOffsets[Bar_ID];
					valid_pmt.ftdc_corr	-= FADCOffsets[Bar_ID];
					bars[Bar_ID].left 	= valid_pmt;
				}
				else	bars[Bar_ID].right 	= valid_pmt;
			}
		}

		// Now let's loop over all of the bars!
		map<int,Bar>::iterator bar_it;
		for( bar_it = bars.begin() ; bar_it != bars.end() ; ++bar_it){
			int Bar_ID 	= bar_it->first;
			Bar this_bar 	= bar_it->second;
			// Grab the individual PMTs:
			PMT left 	= this_bar.left;
			PMT right 	= this_bar.right;
			if( left.amp == 0 || right.amp == 0 		) continue; // no pair
			if( left.adc == 0 || right.adc == 0		) continue;
			if( left.ftdc == 0 || right.ftdc == 0		) continue;
			if( left.tdc_corr == 0 || right.tdc_corr == 0 	) continue;
			if( left.tdc == 0 || right.tdc == 0 		) continue;
			if( left.sector != right.sector			) continue;
			if( left.layer != right.layer			) continue;
			if( left.component != right.component		) continue;
			if( fabs(left.PMT_ID-right.PMT_ID)>1		) continue;

			// Set the bar information now that we have a valid bar pairing:
			this_bar.Bar_ID 	= Bar_ID;
			this_bar.sector 	= left.sector;
			this_bar.layer 		= left.layer;
			this_bar.component 	= left.component;

			int sector 	= this_bar.sector;
			int layer 	= this_bar.layer;
			int component 	= this_bar.component;
			h1_tdiff[sector-1][layer-1][component-1]	.Fill( left.tdc_corr - right.tdc_corr );
			h1_tdiff_fadc[sector-1][layer-1][component-1]	.Fill( left.ftdc_corr - right.ftdc_corr );
			
		}
	} // end loop over events

	// Create output canvas to draw on
	TFile * outFile = new TFile("offsets.root","RECREATE");
	TCanvas c0("c0","c0",700,900);
	c0.Print("offsets.pdf(");
	ofstream outfile;
	outfile.open("lr_offsets.txt");
	std::setprecision(9);
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			TCanvas cSLC(Form("sector_%i_layer_%i",(sector+1),(layer+1)),Form("sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			cSLC.Divide(2,7);
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				// set titles:
				cout << (sector+1) << " " << (layer+1) << " " << (comp+1)  << "\n";
				h1_tdiff[sector][layer][comp].SetTitle(Form("TDC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_tdiff_fadc[sector][layer][comp].SetTitle(Form("FADC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_tdiff[sector][layer][comp].Write();
				h1_tdiff_fadc[sector][layer][comp].Write();

				if( h1_tdiff[sector][layer][comp].Integral() == 0 || 
					h1_tdiff_fadc[sector][layer][comp].Integral() == 0 ){
						outfile << (sector+1) << "\t" << (layer+1) << "\t" << (comp+1) << "\t"
							<< 0 << "\t" << 0 << "\t" << 0 << "\t"
							<< 0 << "\t" << 0 << "\t" << 0 << "\n";
						continue;
				}

				double offset, veff, width;
				offsetFit( &h1_tdiff[sector][layer][comp] 	, &cSLC , 2*comp+1 , (sector+1) , offset, veff, width );
				
				outfile << (sector+1) << "\t" << (layer+1) << "\t" << (comp+1) << "\t"
					<< offset << "\t" << veff << "\t" << width << "\t";

				offsetFit( &h1_tdiff_fadc[sector][layer][comp]   , &cSLC , 2*comp+2 , (sector+1) , offset, veff, width );

				outfile << offset << "\t" << veff << "\t" << width << "\n";
				
				/*

						outfile_amp << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << 
							0 << " " << 0 << " " << 0 << " " << 
							0 << " " << 0 << " " << 0 << " " << 
							0 << " " << 0 << " " << 0 << " " << 
							0 << " " << 0 << " " << 0 << "\n";

						continue;
				}
				double par1, par2, par3, par4, par5, par6;
				double par1e, par2e, par3e, par4e, par5e, par6e;
				
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
			*/
			}
			cSLC.Print("offsets.pdf");
		}
	}
	outfile.close();
	c0.Print("offsets.pdf)");
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

	std::string path = string(getenv("BANDSOFT_TOOLS_DIR"))+"/include/calibrations";
	if( SPRING2019 ) path += "/spring2019";
	else if( FALL2019_WINTER2020 ) path += "/fall2019";
	else{ cerr << "cannot load file\n"; exit(-1); }
	cout << path << "\n";
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

	// Load R PMT file:
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

void offsetFit( TH1D * hist , TCanvas * c , int cd, int sector, double &offset, double &veff, double &width ){
	double avg = 0;
	int ctr = 0;
	for( int bin = 1 ; bin <= hist->GetXaxis()->GetNbins(); ++bin ){
		int content = hist->GetBinContent(bin);
		if( content > 1 ){
			avg += content;
			++ctr;
		}
	}
	avg /= ctr;

	double xLower = 0;
	double xUpper = 0;
	double avgl = 0;
	double avgu = 0;
	c->cd(cd);
	intercept( hist , avg , sector, xLower , xUpper , avgl, avgu );
	
	offset = (xLower+xUpper)/2.;
	width = (xUpper-xLower);
	veff = bar_lengths[sector-1]/(width/2.);


	c->cd(cd);
	hist->Draw();

	TLine * lLine = new TLine(xLower,0,xLower,hist->GetMaximum());
	TLine * uLine = new TLine(xUpper,0,xUpper,hist->GetMaximum());
	TLine * mLine = new TLine(offset,0,offset,hist->GetMaximum());

	TLine * aLine = new TLine(xLower,avg,xUpper,avg);
	TLine * hLine = new TLine(xLower,avgl,xUpper,avgl);
	TLine * iLine = new TLine(xLower,avgu,xUpper,avgu);

	lLine->SetLineColor(2);
	lLine->SetLineWidth(1);
	lLine->SetLineStyle(2);
	uLine->SetLineColor(2);
	uLine->SetLineWidth(1);
	uLine->SetLineStyle(2);
	hLine->SetLineColor(4);
	hLine->SetLineWidth(1);
	hLine->SetLineStyle(2);
	iLine->SetLineColor(4);
	iLine->SetLineWidth(1);
	iLine->SetLineStyle(2);
	aLine->SetLineColor(1);
	aLine->SetLineWidth(1);
	aLine->SetLineStyle(2);
	mLine->SetLineColor(2);
	mLine->SetLineWidth(1);
	lLine->Draw("SAME");
	uLine->Draw("SAME");
	hLine->Draw("SAME");
	iLine->Draw("SAME");
	mLine->Draw("SAME");
	aLine->Draw("SAME");
	c->SetLogy();
	c->Modified(); c->Update();
	return;
}
void intercept( TH1D * hist , double avg , int sector, double& lower , double& upper , double& al, double& au ){
	
	al = 0, au = 0;
	if( sector == 1  ){
		al = 0.2*avg;
		au = 0.9*avg;
	}
	else if( sector == 3 || sector == 4){
		al = 0.4*avg;
		au = 0.9*avg;
	}
	else if( sector == 2 || sector == 5){
		al = 0.2*avg;
		au = 0.6*avg;
	}
	else if( sector == 5 ){
		al = 0.4*avg;
		au = 0.8*avg;
	}
	else{ cerr << "invalid\n"; exit(-1); }

	double y1l = al;
	double y2l = au;
	int bin1l = hist->FindFirstBinAbove( y1l );
	int bin2l = hist->FindFirstBinAbove( y2l );
	double x1l =  hist->GetXaxis()->GetBinCenter(bin1l);
	double x2l =  hist->GetXaxis()->GetBinCenter(bin2l);

	TFitResultPtr ptrl = hist->Fit("pol1","QESR+","",x1l,x2l);
	
	double p0l = ptrl->Parameter(0);
	double p1l = ptrl->Parameter(1);
	lower = - p0l / p1l;

	//y1l = hist->GetBinContent(bin1l);
	//y2l = hist->GetBinContent(bin2l);

	//double ml = (y2l-y1l)/(x2l-x1l);
	//double bl = y1l - ml*x1l;
	//lower = -bl / ml;

	double y1u = al;
	double y2u = au;
	int bin1u = hist->FindLastBinAbove( y1u );
	int bin2u = hist->FindLastBinAbove( y2u );
	double x1u =  hist->GetXaxis()->GetBinCenter(bin1u);
	double x2u =  hist->GetXaxis()->GetBinCenter(bin2u);

	TFitResultPtr ptru = hist->Fit("pol1","QESR+","",x2u,x1u);

	double p0u = ptru->Parameter(0);
	double p1u = ptru->Parameter(1);
	upper = - p0u / p1u;

	//y1u = hist->GetBinContent(bin1u);
	//y2u = hist->GetBinContent(bin2u);

	//double mu = (y2u-y1u)/(x2u-x1u);
	//double bu = y1u - mu*x1u;
	//upper = -bu / mu;

	//cout << avg << " " << y1l << " " << y2l << " " << y1u << " " << y2u << "\n";
	//cout << x1l << " " << x2l << " " << x1u << " " << x2u << "\n";
	//cout << lower << " " << upper << "\n\n";

	//double thres = avg*0.3;
	//lower = hist->GetXaxis()->GetBinCenter(
	// 		hist->FindFirstBinAbove( thres ) );
	//upper = hist->GetXaxis()->GetBinCenter(
	//			hist->FindLastBinAbove( thres ) );

	return;
}
