#include <iostream>
#include <fstream>
#include <sstream>

#include "reader.h"
#include "bank.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "bandreco.h"
#include "constants.h"


const int maxNeutrons	= 200;

using namespace std;


void offsetFit( TH1D * hist , TCanvas * c , int cd, int layer, double &amp, double &mean, double &sigma );

int period = -1;

int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [calibration period] [inputFile1] [inputFile2] ... \n";
		cerr << "\t\t[calibration peroid] -- 0 = SPRING2019 / 1 == FALL2019_WINTER2020\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}

	period = atoi(argv[1]); // TODO just get it from run number directly

	// Initialize our BAND reconstruction engine:
	BANDReco * BAND = new BANDReco();
	BAND->setPeriod(period);
	BAND->readTW();			// TW calibration values for each PMT
	BAND->readLROffset();		// (L-R) offsets for each bar
	BAND->readPaddleOffset();	// bar offsets relative to bar 2X7 in each layer X
	BAND->readLayerOffset();	// layer offsets relative to layer 5
	BAND->readGeometry();		// geometry table for each bar
	BAND->readEnergyCalib();	// energy calibration for Adc->MeVee 
	BAND->readStatus();		// status table for 0,1 = bad,good bar


	// Create histograms for storage
	TH1D h1_loff[5][6][7];
	TH1D h1_loff_fadc[5][6][7];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < BAND->slc[layer][sector] ; ++comp ){
				h1_loff[sector][layer][comp]
						= TH1D(	Form("loff_%i_%i_%i",sector+1,layer+1,comp+1),
							Form("loff_%i_%i_%i",sector+1,layer+1,comp+1),8000,-200,200);
				h1_loff_fadc[sector][layer][comp]
						= TH1D(	Form("loff_fadc_%i_%i_%i",sector+1,layer+1,comp+1),
							Form("loff_fadc_%i_%i_%i",sector+1,layer+1,comp+1),8000,-200,200);
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
		hipo::bank	event_info		(factory.getSchema("REC::Event"		));
		hipo::event 	readevent;

		int event_counter = 0;
		// Loop over all events in file
		while(reader.next()==true){
			BAND->Clear();

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if( event_counter > 500000 ) break;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(band_adc);
			readevent.getStructure(band_tdc);
			readevent.getStructure(run_config);
			readevent.getStructure(event_info);


			// Form the PMTs and Bars for BAND:
			BAND->createPMTs( &band_adc, &band_tdc, &run_config );
			BAND->createBars();

			// Now let's loop over all of the bars to fill our histograms!
			map<int,Bar> candidate_bars = BAND->getCandidateBars();
			map<int,Bar>::iterator bar_it;
			for( bar_it = candidate_bars.begin() ; bar_it != candidate_bars.end() ; ++bar_it){
				int Bar_ID 	= bar_it->first;
				Bar this_bar 	= bar_it->second;

				int sector 	= this_bar.sector;
				int layer 	= this_bar.layer;
				int component 	= this_bar.component;
				if( sector!=2 || component!=7 ) continue; // only align paddles 2X1 relative to 257

				double mean_tdc 	= this_bar.Tof;
				double mean_ftdc 	= this_bar.TofFtdc;

	
				// get the reference bar:
				int ref_ID = 200 + 50 + 7;
				if( candidate_bars.count(ref_ID) == 0 ) continue;
				double reftime 		= candidate_bars[ref_ID].Tof;
				double reftime_ftdc 	= candidate_bars[ref_ID].TofFtdc;

				h1_loff[sector-1][layer-1][component-1]		.Fill( mean_tdc - 	reftime		);
				h1_loff_fadc[sector-1][layer-1][component-1]	.Fill( mean_ftdc - 	reftime_ftdc	);
			}

		} // end loop over events
	} // end loop over files

	// Create output canvas to draw on
	TFile * outFile = new TFile("off_layer.root","RECREATE");
	TCanvas c0("c0","c0",700,900);
	c0.Print("off_layer.pdf(");
	ofstream outfile;
	outfile.open("layer_offsets.txt");
	std::setprecision(9);
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			TCanvas cSLC(Form("sector_%i_layer_%i",(sector+1),(layer+1)),Form("sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			cSLC.Divide(2,7);
			for( int comp = 0 ; comp < BAND->slc[layer][sector] ; ++comp ){
				// set titles:
				cout << (sector+1) << " " << (layer+1) << " " << (comp+1)  << "\n";
				h1_loff[sector][layer][comp].SetTitle(Form("TDC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_loff_fadc[sector][layer][comp].SetTitle(Form("FADC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_loff[sector][layer][comp].Write();
				h1_loff_fadc[sector][layer][comp].Write();

				if( h1_loff[sector][layer][comp].Integral() == 0 || 
					h1_loff_fadc[sector][layer][comp].Integral() == 0 ||
					(sector+1 == 2 && comp+1 == 7 && layer+1==5 ) ){
						outfile << (sector+1) << "\t" << (layer+1) << "\t" << (comp+1) << "\t"
							<< 0 << "\t" << 0 << "\t" << 0 << "\t"
							<< 0 << "\t" << 0 << "\t" << 0 << "\n";
						continue;
				}

				double amp, mean, sigma;
				amp = mean = sigma = 0.;
				offsetFit( &h1_loff[sector][layer][comp] 	, &cSLC , 2*comp+1 , (layer+1) , amp, mean, sigma );
				outfile << (sector+1) << "\t" << (layer+1) << "\t" << (comp+1) << "\t"
					<< amp << "\t" << mean << "\t" << sigma << "\t";


				offsetFit( &h1_loff_fadc[sector][layer][comp]   , &cSLC , 2*comp+2 , (layer+1) , amp, mean, sigma );
				outfile << amp << "\t" << mean << "\t" << sigma << "\n";


				
			}
			cSLC.Print("off_layer.pdf");
		}
	}
	outfile.close();
	c0.Print("off_layer.pdf)");
	outFile->Close();
	delete outFile;
	

	return 0;
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
