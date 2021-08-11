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
#include "TLine.h"


const int maxNeutrons	= 200;

using namespace std;


void offsetFit( TH1D * hist , TCanvas * c , int cd, int sector, double &offset, double &veff, double &width );
double bar_lengths[5] = {163.7,201.9,51.2,51.2,201.9};
void intercept( TH1D * hist , double avg , int sector, double& lower , double& upper , double& al, double& au );

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
	BAND->readTW();
	BAND->readStatus();


	// Create histograms for storage
	TH1D h1_tdiff[5][6][7];
	TH1D h1_tdiff_fadc[5][6][7];
	for( int sector = 0 ; sector < 5 ; ++sector ){
		for( int layer = 0 ; layer < 6 ; ++layer ){
			for( int comp = 0 ; comp < BAND->slc[layer][sector] ; ++comp ){
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
				double tdiff		= this_bar.Tdiff;
				double tdiff_ftdc	= this_bar.TdiffFtdc;
	

				h1_tdiff[sector-1][layer-1][component-1]	.Fill(  tdiff 		);
				h1_tdiff_fadc[sector-1][layer-1][component-1]	.Fill(  tdiff_ftdc 	);
			}

		} // end loop over events
	} // end loop over files

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
			for( int comp = 0 ; comp < BAND->slc[layer][sector] ; ++comp ){
				// set titles:
				cout << (sector+1) << " " << (layer+1) << " " << (comp+1)  << "\n";
				h1_tdiff[sector][layer][comp].SetTitle(Form("TDC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_tdiff_fadc[sector][layer][comp].SetTitle(Form("FADC Sector %i, Layer %i, Component %i",(sector+1),(layer+1),(comp+1)));
				h1_tdiff[sector][layer][comp].Write();
				h1_tdiff_fadc[sector][layer][comp].Write();

				if( h1_tdiff[sector][layer][comp].Integral() == 0 || 
					h1_tdiff_fadc[sector][layer][comp].Integral() == 0 || (layer+1) == 6 ){
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
