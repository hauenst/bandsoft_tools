#ifndef __READHIPO_HELPER_H__
#define __READHIPO_HELPER_H__

#include <fstream>
#include <string>
#include <vector>

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"
#include "bank.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"
#include "genpart.h"
#include "TRandom3.h"

using namespace std;

const int maxProtons	= 100;
const int maxNeutrons	= 200;
const int maxGens	= 10;
const int maxParticles	= 100;
const int maxScinHits = 100;
const int maxEcalhits = 100;
const double time_thresBANDhit = 300;
const double adctoMeVee_data = 1100; //conversion for data
const double adctoMeVee_sim = 1E4;//conversion for simulation
const double VERTEX_OFFSET = -3; // [cm]
const double BAND_OFFSET = -3; //{cm]


class shiftsReader {
	public:
		void LoadInitBar( string filename );
		void LoadInitBarFadc( string filename );
		void LoadInitRun( string filename );
		void LoadInitRunFadc( string filename );
		void LoadEffVel( string filename_S6200 , string filename_S6291 );
		void LoadLrOff( string filename_S6200 , string filename_S6291 );
		double * getInitBar(void);
		double * getInitBarFadc(void);
		double * getInitRun(void);
		double * getInitRunFadc(void);

		double * getFadcEffVel	(int Runno);
		double * getTdcEffVel	(int Runno);
		double * getFadcLrOff	(int Runno);
		double * getTdcLrOff	(int Runno);
	private:
		double InitBar[670] = {0.};
		double InitBarFadc[670] = {0.};
		double InitRun[100000] = {0.};
		double InitRunFadc[100000] = {0.};
		// Spring 2019 6200-6291 GeV/c constants
		double S6200FadcEffVel		[670] = {0.};
		double S6200TdcEffVel		[670] = {0.};
		double S6200TdcLrOffsets	[670] = {0.};
		double S6200FadcLrOffsets	[670] = {0.};
		// Spring 2019 6291-INF constants
		double S6291FadcEffVel		[670] = {0.};
		double S6291TdcEffVel		[670] = {0.};
		double S6291TdcLrOffsets	[670] = {0.};
		double S6291FadcLrOffsets	[670] = {0.};
};


int getRunNumber( string filename );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun, std::map<int,double> &bar_x, std::map<int,double> &bar_y, std::map<int,double> &bar_z,
			std::map<int,double> &bar_edep , int hotfix=0, double* s6200_fadc_lroffset=NULL , double* s6200_tdc_lroffset=NULL,
			double* s6291_fadc_lroffset=NULL,	double* s6291_tdc_lroffset=NULL,
			double* s6200_fadc_effvel=NULL,	double* s6200_tdc_effvel=NULL,
	       		double* s6291_fadc_effvel=NULL,	double* s6291_tdc_effvel=NULL	);
void getElectronInfo( BParticle particles, BCalorimeter calorimeter, BScintillator scintillator, hipo::bank DC_Track, hipo::bank DC_Traj,
			int pbankIndex, clashit &electron,
			double starttime , int thisRun , double Ebeam );
void getTaggedInfo( clashit eHit, bandhit nHit[maxNeutrons], taghit tag[maxNeutrons] ,
		double Ebeam , int nMult );
bool goodNeutronEvent(bandhit hits[maxNeutrons], int nMult, int& leadindex, int mcdataselect);
void getMcInfo( hipo::bank gen_particles , hipo::bank gen_info , genpart mcParts[maxGens] ,
		double &starttime, double &weight, double &Ebeam , int &genMult );
void getScinHits( BScintillator scintillator, int pindex[maxScinHits], int detid[maxScinHits], double energy[maxScinHits], double time[maxScinHits],
			    TVector3 posVector[maxScinHits], double path[maxScinHits], int status[maxScinHits], int posIndex[maxParticles], int posMult, int &scinHits);
void getParticleInfo( BParticle particles, int pid[maxParticles], TVector3 momentum[maxParticles], TVector3 vertex[maxParticles],
								double time[maxParticles], int charge[maxParticles], double beta[maxParticles], double chi2pid[maxParticles], int status[maxParticles] , int index[maxParticles], int& multiplicity );
void getBANDBarGeometry(string filename, std::map<int,double> &bar_x, std::map<int,double> &bar_y, std::map<int,double> &bar_z);
void getBANDEdepCalibration(string filename, std::map<int,double> &bar_edep);
//Parametrization from Giovanni (GWU) and FX based on double pion analysis in RGA
void smearRGA(TVector3 &vpar);
void recalculate_clashit_kinematics(clashit &input_ehit, double Ebeam, TVector3 &smeared_electron);

#endif
