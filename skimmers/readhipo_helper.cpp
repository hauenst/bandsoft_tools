#include "readhipo_helper.h"


double FADC_GLOBSHIFT[600] = {0.};
double FADC_RUNBYRUNSHIFT[100000] = {0.};
double TDC_GLOBSHIFT[600] = {0.};
double TDC_RUNBYRUNSHIFT[100000] = {0.};

void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun){
	
	if( band_hits.getRows() > maxNeutrons ) return; // not interested in events with more than max # BAND hits for now
	for( int hit = 0 ; hit < band_hits.getRows() ; hit++ ){
		//if( band_hits.getStatus(hit) != 0 ) continue;	// not interested in an event that has a veto hit

		// Set the hits information
		hits[hit].setSector		(band_hits.getSector		(hit)			);
		hits[hit].setLayer		(band_hits.getLayer		(hit)			);
		hits[hit].setComponent		(band_hits.getComponent		(hit)			);
		hits[hit].setBarID		(band_hits.getBarKey		(hit)			);
		hits[hit].setEdep		(band_hits.getEnergy		(hit)			);
		hits[hit].setTof		(band_hits.getTime		(hit) - starttime	);
		hits[hit].setTofFadc		(band_hits.getTimeFadc		(hit) - starttime	);
		hits[hit].setTdiff		(band_hits.getDifftimeTdc	(hit)			);
		hits[hit].setTdiffFadc		(band_hits.getDifftimeFadc	(hit)			);
		hits[hit].setX			(band_hits.getX			(hit)			);
		hits[hit].setY			(band_hits.getY			(hit)			);
		hits[hit].setZ			(band_hits.getZ			(hit)			);
		hits[hit].setStatus		(band_hits.getStatus		(hit)			);


		// Using the band hit struct, get the raw hit PMT information to use later
		int rawhit_idxL = band_hits.getLpmtindex(hit);
		int rawhit_idxR = band_hits.getRpmtindex(hit);
		// 	Get the raw hit information corresponding to the band hit above
		hits[hit].setRawLtdc		(band_rawhits.getFloat( 7 , rawhit_idxL ) 		);	
        	hits[hit].setRawRtdc		(band_rawhits.getFloat( 7 , rawhit_idxR ) 		);
        	hits[hit].setRawLtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxL ) 		);
        	hits[hit].setRawRtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxR ) 		);
        	hits[hit].setRawLtfadc		(band_rawhits.getFloat( 8 , rawhit_idxL ) 		);
        	hits[hit].setRawRtfadc		(band_rawhits.getFloat( 8 , rawhit_idxR ) 		);
        	hits[hit].setRawLamp		(band_rawhits.getFloat( 6 , rawhit_idxL )		);
        	hits[hit].setRawRamp		(band_rawhits.getFloat( 6 , rawhit_idxR )		);
        	hits[hit].setRawLadc		(band_rawhits.getFloat( 5 , rawhit_idxL )		);
        	hits[hit].setRawRadc		(band_rawhits.getFloat( 5 , rawhit_idxR )		);

		// Using the rawhit struct, get the raw PMT information to use later
		int pmtTdcL	= band_rawhits.getInt( 10 , rawhit_idxL );
		int pmtAdcL	= band_rawhits.getInt( 11 , rawhit_idxL );
		int pmtTdcR	= band_rawhits.getInt( 10 , rawhit_idxR );
		int pmtAdcR	= band_rawhits.getInt( 11 , rawhit_idxR );
		//	Get the raw pmt information corresponding to the band hit above
		hits[hit].setPmtLtdc		(band_tdc.getInt( 4 , pmtTdcL )		);
		hits[hit].setPmtRtdc		(band_tdc.getInt( 4 , pmtTdcR )		);
		hits[hit].setPmtLtfadc		(band_adc.getFloat( 6 , pmtAdcL )	);
		hits[hit].setPmtRtfadc		(band_adc.getFloat( 6 , pmtAdcR )	);
		hits[hit].setPmtLamp		(band_adc.getInt( 5 , pmtAdcL )		);
		hits[hit].setPmtRamp		(band_adc.getInt( 5 , pmtAdcR )		);
		hits[hit].setPmtLadc		(band_adc.getInt( 4 , pmtAdcL )		); 
		hits[hit].setPmtRadc		(band_adc.getInt( 4 , pmtAdcR )		);
		hits[hit].setPmtLped		(band_adc.getInt( 7 , pmtAdcL )		); 
		hits[hit].setPmtRped		(band_adc.getInt( 7 , pmtAdcR )		); 

		// Save how many neutron hits we have
		mult++;
	}
	
}

int getRunNumber( string filename ){
	string parsed = filename.substr( filename.find("inc_") );
	string moreparse = parsed.substr(4,6);
	cout << "\t*Intepreted run number from file name: " << stoi(moreparse) << "\n";
        return stoi(moreparse);
}

void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime ){
	if( eventInfo.getRows() != 1 ){ 
		cerr << "getEventInfo::NotImplementedFunction\n"; 
		exit(-1); 
	}
	//integrated_charge       += (double)eventInfo.getBCG(0); 	// not calibrated currently
	//livetime 		= (double)eventInfo.getLT(0);		// not calibrated currently
	starttime		= (double)eventInfo.getSTTime(0);
	return;
}
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status ){
	pid 		= particles.getPid(0);
	momentum 	= particles.getV3P(0);
	vertex		= particles.getV3v(0);
	time		= particles.getVt(0);
	charge		= particles.getCharge(0);
	beta		= particles.getBeta(0);
	chi2pid		= particles.getChi2pid(0);
	status		= particles.getStatus(0);
	return;
}
void getProtonInfo( BParticle particles, double pid[maxProtons], TVector3 momentum[maxProtons], TVector3 vertex[maxProtons],
			double time[maxProtons], double charge[maxProtons], double beta[maxProtons], double chi2pid[maxProtons], double status[maxProtons] , int& multiplicity ){
	// Takes the first proton in the bank
	multiplicity = 0;
	for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
		pid[multiplicity] 		= particles.getPid(row);
		charge[multiplicity]		= particles.getCharge(row);
		if( charge[multiplicity] == 1 ){
			momentum 	[multiplicity]	= particles.getV3P(row);
			vertex		[multiplicity]	= particles.getV3v(row);
			time		[multiplicity]	= particles.getVt(row);
			beta		[multiplicity]	= particles.getBeta(row);
			chi2pid		[multiplicity]	= particles.getChi2pid(row);
			status		[multiplicity]	= particles.getStatus(row);
			multiplicity ++;
		}
	}
}
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot){
			// Anymore fiducial cuts that are needed for the electron:
			// 	E/p cut
			// 	vertex cut
			// 	minimum momentum  cut / maximum  momentum  cut
			// 	electron ToF cut
			// 	minimum W cut
			//	chi2 cut
			//	pcal energy cut?
			//

	if( pid != 11 || charge != -1 ) return false;
	//if( lV < 2 || lW < 2 ) return false;
	//if( momentum.Mag() < 1 || momentum.Mag() > 4.2 ) return false;
	//if( chi2pid == 0 || chi2pid > 1 ) return false;
	//if( E_tot/momentum.Mag() > 0.4 || E_tot/momentum.Mag() < 0.1 ) return false;
	//if( vertex.X() < -2 || vertex.X() > 2) return false;
	//if( vertex.Y() < -2 || vertex.Y() > 2) return false;
	//if( vertex.Z() < -7 || vertex.Z() > 2) return false;
	//if( time < 15 ) return false;
	//if( momentum.Mag() < 2 || momentum.Mag() > 10.6 ) return false;
	//if( E_tot / momentum.Mag() < 0.15 || E_tot / momentum.Mag() > 0.3 ) return false;

	return true;
}
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult ){
	//if( momentum.Mag() == 0 || momentum.Mag() > 4.2 ) return false;
	//if( beta <= 0 || beta > 1 ) return false;
	//if( mult != 1 ) return false;
	//if( momentum.Theta() == 0 ) return false;
	//if( beta <= 0 ) return false;
	//if( del_vertex.Z() < 0 || del_vertex.Z() > 5 ) return false;
	//if( chi2pid < -3 || chi2pid > 3 ) return false;
	
	
	return true;
}

bool pointsToBand(double theta,double phi,double z_m){
	//double z = z_m*100; // from m to cm
	double z = z_m;

	// Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
	double thickness  = 7.2;                                // thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

	// Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = (zUpst + 5*thickness) - z_m;

	double rho   = zDown/cos(theta);
	double xDown = rho*sin(theta)*cos(phi);
	double yDown = rho*sin(theta)*sin(phi);

	double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	// Sector boundaries
	double topSec1  = globalY + 13*thickness;
	double topSec2  = globalY + 10*thickness;
	double topSec34 = globalY +  3*thickness;
	double topSec5  = globalY -  3*thickness;
	double downSec5 = globalY -  5*thickness;

	if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

	if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
			(yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
			(yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
			(yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
	  )
		return 1;

	return 0;
}

void LoadGlobalShift( string filename_tdc , string filename_fadc ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename_fadc);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_GLOBSHIFT[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();

	f.open(filename_tdc);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		TDC_GLOBSHIFT[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

void LoadRunByRunShift(){
	ifstream f;
	int runnum;
	double pol0, height, mean, sig, temp;

	f.open("../include/runByrun_offset_fadc-10082019.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
	f.open("../include/runByrun_offset_tdc_firstIter-04132020.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		TDC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

