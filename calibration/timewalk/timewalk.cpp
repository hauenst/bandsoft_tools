#include <iostream>

#include "reader.h"
#include "bank.h"

#include "TString.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
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

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [calibration period] [inputFile] \n";
		cerr << "\t\t[calibration peroid] -- 0 = SPRING2019 / 1 == FALL2019_WINTER2020\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}




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


	// Create histograms for storage
	TH2D h2_tadc[5][6][7][2];
	TH2D h2_tamp[5][6][7][2];
	TH2D h2_tadc_overfill[5][6][7][2];
	TH2D h2_tamp_overfill[5][6][7][2];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < slc[layer][sector] ; ++comp ){
				for( int order = 0 ; order < 2 ; ++order ){
					h2_tadc[sector][layer][comp][order] 
							= TH2D(	Form("adc_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("adc_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),300,0,30000,350,140,175);
					h2_tadc_overfill[sector][layer][comp][order] 
							= TH2D(Form("adcO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("adcO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),300,15000,45000,350,140,175);
					h2_tamp[sector][layer][comp][order] 
							= TH2D(Form("amp_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("amp_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),250,0,5000,350,140,175);
					h2_tamp_overfill[sector][layer][comp][order] 
							= TH2D(Form("ampO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),
								Form("ampO_%i_%i_%i_%i",sector+1,layer+1,comp+1,order),250,0,5000,350,140,175);
				}
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
		//if( event_counter > 50000 ) break;
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
			
			//if( amp > 1500 && PMT_ID != REFID ) continue;

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
				if( amp >= 4095 ){
					h2_tadc_overfill[sector-1][layer-1][component-1][order].Fill( adc , tdc - TREF );
					h2_tamp_overfill[sector-1][layer-1][component-1][order].Fill( amp , tdc - TREF );
				}
				else{
					h2_tadc[sector-1][layer-1][component-1][order].Fill( adc , tdc - TREF );
					h2_tamp[sector-1][layer-1][component-1][order].Fill( amp , tdc - TREF );
				}
			}
			
		}
	} // end loop over events

	// Create output canvas to draw on
	TFile * outFile = new TFile("init_timewalk.root","RECREATE");
	TCanvas c0("c0","c0",700,900);
	TCanvas c1("c1","c1",700,900);
	c0.Print("init_adc_timewalk.pdf(");
	c1.Print("init_amp_timewalk.pdf(");
	ofstream outfile_adc[2];
	ofstream outfile_adcR[2];
	outfile_adc[0].open("init_paramsadcL.txt");
	outfile_adc[1].open("init_paramsadcR.txt");
	ofstream outfile_amp[2];
	outfile_amp[0].open("init_paramsampL.txt");
	outfile_amp[1].open("init_paramsampR.txt");
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
					//cout << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << (order) << "\n";
					h2_tadc[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tadc[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h2_tadc_overfill[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tadc_overfill[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp_overfill[sector][layer][comp][0].SetTitle(Form("Sector %i, Layer %i, Component %i, L PMT",(sector+1),(layer+1),(comp+1)));
					h2_tamp_overfill[sector][layer][comp][1].SetTitle(Form("Sector %i, Layer %i, Component %i, R PMT",(sector+1),(layer+1),(comp+1)));
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
					fitTW( &h2_tadc[sector][layer][comp][order] , &cSLC_adc,  4*comp+1 + (order)*2 + 1, 
							(sector+1), (layer+1), (comp+1), order, 0,
						       	par1 , par2 , par3, par4, par5, par6,
							par1e, par2e , par3e, par4e, par5e, par6e );
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
							(sector+1), (layer+1), (comp+1), order, 0, 
						       	par1 , par2 , par3, par4, par5, par6,
							par1e, par2e , par3e, par4e, par5e, par6e );
					outfile_amp[order] << (sector+1) << " " << (layer+1) << " " << (comp+1) << " " << 
						par1 << " " << par2 << " " << par3 << " " << par4 << " " << par5 << " " << par6 << " " << 
						par1e << " " << par2e << " " << par3e << " " << par4e << " " << par5e << " " << par6e << "\n";

				}
				//break;
			}
			cSLC_adc.Print("init_adc_timewalk.pdf");
			cSLC_amp.Print("init_amp_timewalk.pdf");
		}
	}
	c0.Print("init_adc_timewalk.pdf)");
	c1.Print("init_amp_timewalk.pdf)");
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


void fitTW(TH2D * hist , TCanvas * c, int cd, int s, int l, int co, int o, double cut, 
		double &par1, double &par2 , double &par3, double &par4, double &par5, double &par6,
		double &par1_err, double &par2_err , double &par3_err, double &par4_err, double &par5_err, double &par6_err){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	walkCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, cut, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);
	TF1 * model = new TF1("timeWalk",wlk,0,20000,3);
	model->SetParameter(0,150);
	model->SetParameter(1,80);
	model->SetParameter(2,0.5);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par1 = ptr->Parameter(0);
	par1_err = ptr->ParError(0);

	par2 = ptr->Parameter(1);
	par2_err = ptr->ParError(1);

	par3 = ptr->Parameter(2);
	par3_err = ptr->ParError(2);

	//par4 = ptr->Parameter(3);
	//par4_err = ptr->ParError(3);

	//par5 = ptr->Parameter(4);
	//par5_err = ptr->ParError(4);

	//par6 = ptr->Parameter(5);
	//par6_err = ptr->ParError(5);

	c->cd(cd);
	//gStyle->SetTitleW(0.6);
	g->SetTitle(Form("SLCO: %i %i %i %i",s,l,co,o));
	g->SetMarkerStyle(20);
	g->Draw("AP");
	TLine * line = new TLine(4095,0,4059,200);
	line->Draw("SAME");
	c->Modified(); c->Update();

	//delete model;
	//delete g;
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
		int nEvents = 400;
		//if( currBin_x > 3000 ) nEvents = 200;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr , nEvents, maxBin );
		currBin += step ;
		if( xPt < widthCut) continue;
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
	return p[0] + p[1] / pow(var,p[2]);
	//return p[0] + p[1]*(var) + p[2]*pow(var,2) + p[3]*pow(var,3) + p[4]*pow(var,4) + p[5]*pow(var,5);
	//return p[0] + p[1] / pow(var,p[2]);
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
		TFitResultPtr f = pj->Fit("gaus","QESR","",-300,300);
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
