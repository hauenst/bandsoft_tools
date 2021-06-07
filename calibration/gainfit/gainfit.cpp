#include <iostream>

#include "reader.h"
#include "bank.h"

#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"

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

double getTriggerPhase( long timeStamp ) ;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [inputFile] \n\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}

	// Create histograms for storage
	TH1D h_adcspec[5][6][7][2];
	TH1D h_ampspec[5][6][7][2];
	TH1D h_adcspec_overfill[5][6][7][2];
	TH1D h_ampspec_overfill[5][6][7][2];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				for( int order = 0 ; order < 2 ; ++order ){
					h_adcspec[sector][layer][comp][order] 
							= TH1D(Form("adc_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),Form("adc_%i_%i_%i_%i",sector+1,layer+1,comp+1,order+1),300,0,30000);
					h_adcspec_overfill[sector][layer][comp][order] 
							= TH1D(Form("amp_%i_%i_%i_%i",sector,layer,comp,order),Form("amp_%i_%i_%i_%i",sector,layer,comp,order),300,0,30000);
					h_ampspec[sector][layer][comp][order] 
							= TH1D(Form("adcO_%i_%i_%i_%i",sector,layer,comp,order),Form("adcO_%i_%i_%i_%i",sector,layer,comp,order),250,0,5000);
					h_ampspec_overfill[sector][layer][comp][order] 
							= TH1D(Form("ampO_%i_%i_%i_%i",sector,layer,comp,order),Form("ampO_%i_%i_%i_%i",sector,layer,comp,order),250,0,5000);
				}
			}
		}
	}
	


	// Setup hipo reading for this file
	TString inputFile = argv[1];
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

			int tdc	= band_tdc.getInt( 4 , row ) * 0.02345 - phaseCorr;

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
		map<int,vector<PMT>>::iterator it;
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
				int amp		= this_pmt.amp;

				if( amp >= 4095 ){
					h_adcspec_overfill[sector-1][layer-1][component-1][order].Fill( adc );
					h_ampspec_overfill[sector-1][layer-1][component-1][order].Fill( amp );
				}
				else{
					h_adcspec[sector-1][layer-1][component-1][order].Fill( adc );
					h_ampspec[sector-1][layer-1][component-1][order].Fill( amp );
				}
			}
			
		}
	} // end loop over events

	// Create output canvas to draw on
	TCanvas c0("c0","c0",700,900);
	TCanvas c1("c1","c1",700,900);
	c0.Print("adc_cosmics.pdf(");
	c1.Print("amp_cosmics.pdf(");
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			TCanvas cSLC_adc(Form("adc_sector_%i_layer_%i",(sector+1),(layer+1)),Form("adc_sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			TCanvas cSLC_amp(Form("amp_sector_%i_layer_%i",(sector+1),(layer+1)),Form("amp_sector_%i_layer_%i",(sector+1),(layer+1)),700,900);
			cSLC_adc.Divide(2,7);
			cSLC_amp.Divide(2,7);
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				for( int order = 0 ; order < 2 ; ++order ){
					cSLC_adc.cd(comp*2+1 + order);
					h_adcspec[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h_adcspec[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h_adcspec[sector][layer][comp][order].Draw();
					h_adcspec_overfill[sector][layer][comp][order].SetLineColor(2);
					h_adcspec_overfill[sector][layer][comp][order].Draw("SAME");

					cSLC_amp.cd(comp*2+1 + order);
					h_ampspec[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h_ampspec[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h_ampspec[sector][layer][comp][order].Draw();
					h_ampspec_overfill[sector][layer][comp][order].SetLineColor(2);
					h_ampspec_overfill[sector][layer][comp][order].Draw("SAME");
				}
			}
			cSLC_adc.Print("adc_cosmics.pdf");
			cSLC_amp.Print("amp_cosmics.pdf");
		}
	}
	c0.Print("adc_cosmics.pdf)");
	c1.Print("amp_cosmics.pdf)");



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
