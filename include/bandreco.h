#ifndef __BANDRECO_H__
#define __BANDRECO_H__
#include <iostream>
#include <fstream>
#include <sstream>

#include "bank.h"
#include "bandhit.h"
using std::map;
using std::vector;
using std::cerr;
using std::cout;
using std::string;

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
	double tdc	= 0;	// Phase-correction & 0.02345 conversion
	double tdc_corr	= 0;	// TW corrected
	double trigphase= 0;
	int idx_ftdc	= -1;
	int idx_tdc	= -1;
};

struct Bar{
	PMT left;
	PMT right;
	int sector	= 0;
	int layer	= 0;
	int component	= 0;
	int Bar_ID 	= 0;
	double Edep	= 0;//
	double AdcL	= 0;//
	double AdcR	= 0;//
	double AmpL	= 0;//
	double AmpR	= 0;//
	double Tof	= 0;//
	double TofFtdc	= 0;//
	double X	= 0;
	double XFtdc	= 0;
	double Y	= 0;
	double Z	= 0;
	double Tdiff		= 0;//
	double TdiffFtdc	= 0;//
};


class BANDReco{
	
	public:
		// Constructor
		BANDReco(){};
		// Destructor
		~BANDReco(){};

		void Clear();
		void Print();

		void setPeriod(const int period);

		map<int,Bar> getCandidateBars(void){return candidate_bars;};
		map<int,PMT> getCandidatePMTs(void){return candidate_pmts;};

		void createPMTs( const hipo::bank * band_adc , const hipo::bank * band_tdc , const hipo::bank * run_config  );
		void createBars( );
		void storeHits( int& mult , bandhit * hits , const double starttime , const double vtx_z );

		void readTW();
		void readLROffset();
		void readPaddleOffset();
		void readLayerOffset();
		void readGlobalOffset();
		void readGeometry();
		void readEnergyCalib();
		void readStatus();
		
		// Hard-coded parameters:
		const double bar_lengths[5] = {163.7,201.9,51.2,51.2,201.9};
		const int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

	private:
		bool SPRING2019 = false;
		bool FALL2019_WINTER2020 = false;
		bool CALIBRATION = false;
		double BAND_MOTHER_OFFSET = 0.; // [cm]

		double timewalk( const double *x , const double *p);
		double getTriggerPhase( const long timeStamp ) ;
		bool check_bar(const Bar * this_bar);
		bool check_status( const Bar * this_bar);

		map<int,vector<double>> TWParamsAMP;
		map<int,double> TDCOffsets;
		map<int,double> FTDCOffsets;
		map<int,double> TDCVelocity;
		map<int,double> FTDCVelocity;
		map<int,double> TDCPaddle;
		map<int,double> FTDCPaddle;
		map<int,double> TDCLayer;
		map<int,double> FTDCLayer;
		map<int,double> GlobalX;
		map<int,double> GlobalY;
		map<int,double> GlobalZ;
		map<int,double> ADCtoMEV;
		map<int,int> STATUS;
		map<int,double> TDCGlobal; 
		map<int,double> FTDCGlobal; 
		map<int,double> TDCToFRes; 
		map<int,double> FTDCToFRes; 

		// Storage containers for each event:
		map<int,PMT> candidate_pmts;
		map<int,Bar> candidate_bars;

};
#endif
