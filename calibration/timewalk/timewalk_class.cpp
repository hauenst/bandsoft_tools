#include <iostream>
#include <fstream>
#include <sstream>

#include "reader.h"
#include "bank.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "bandreco.h"
#include "constants.h"



using namespace std;


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

int period = -1;

int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [calibration period] [inputFile1] [inputFile2] ... \n";
		cerr << "\t\t[calibration peroid] -- 0 = SPRING2019 / 1 == FALL2019_WINTER2020\n";
		cerr << "\t\t[inputFile] = ____.hipo\n\n";
		return -1;
	}
	
	// Get the period:
	period = atoi(argv[1]); // TODO just get it from run number directly
	bool SPRING2019 = false, FALL2019_WINTER2020 = false;
	(period == 0) ? SPRING2019 = true : FALL2019_WINTER2020 = true;
	// Get the reference photodiode channel:
	int REFID = -1;
	if( SPRING2019 ){
		REFID = 3660; // V16-A
	}
	else if( FALL2019_WINTER2020 ){
		REFID = 4161; // 116B-R
	}

	// Initialize our BAND reconstruction engine:
	BANDReco * BAND = new BANDReco();
	BAND->setPeriod(period);


	// Create histograms for storage
	TH2D h2_tadc[5][6][7][2];
	TH2D h2_tamp[5][6][7][2];
	TH2D h2_tadc_overfill[5][6][7][2];
	TH2D h2_tamp_overfill[5][6][7][2];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < BAND->slc[layer][sector] ; ++comp ){
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

			// Now let's loop over all the PMTs to fill our histograms!
			map<int,PMT> candidate_pmts = BAND->getCandidatePMTs();
			map<int,PMT>::const_iterator pmt_it;

			for( pmt_it = candidate_pmts.begin() ; pmt_it != candidate_pmts.end() ; ++pmt_it){
				int PMT_ID	= pmt_it->first;
				PMT this_pmt 	= pmt_it->second;

				double tdc = 	this_pmt.tdc;
				double amp = 	this_pmt.amp;
				double adc = 	this_pmt.adc;
				int sector = 	this_pmt.sector;
				int layer =	this_pmt.layer;
				int order = 	this_pmt.order;
				int component = this_pmt.component;

				// get the reference pmt:
				if( candidate_pmts.count(REFID) == 0 ) continue;
				double TREF 		= candidate_pmts[REFID].tdc;
				h2_tamp[sector-1][layer-1][component-1][order].Fill( amp , tdc - TREF );
				h2_tadc[sector-1][layer-1][component-1][order].Fill( adc , tdc - TREF );

			}

		} // end loop over events
	} // end loop over files

	
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
			for( int comp = 0 ; comp < BAND->slc[layer][sector] ; ++comp ){
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
