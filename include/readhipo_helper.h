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

using namespace std;

const int maxProtons	= 100;
const int maxNeutrons	= 200;
const int maxGens	= 10;
const double thresBANDhit = 5.;
const double time_thresBANDhit = 300;
const double adctoMeVee_data = 2300; //conversion for data
const double adctoMeVee_sim = 1E4;//conversion for simulation
const double VERTEX_OFFSET = -3; // [cm]

class shiftsReader {
	public:
		void LoadInitBar( string filename );
		void LoadInitBarFadc( string filename );
		void LoadInitRun( string filename );
		void LoadInitRunFadc( string filename );
		double * getInitBar(void);
		double * getInitBarFadc(void);
		double * getInitRun(void);
		double * getInitRunFadc(void);
	private:
		double InitBar[600] = {0.};
		double InitBarFadc[600] = {0.};
		double InitRun[100000] = {0.};
		double InitRunFadc[100000] = {0.};
};


int getRunNumber( string filename );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun);
void getElectronInfo( BParticle particles, BCalorimeter calorimeter, BScintillator scintillator, hipo::bank DC_Track, hipo::bank DC_Traj,
			clashit &electron,
			double starttime , int thisRun , double Ebeam );
void getTaggedInfo( clashit eHit, bandhit nHit[maxNeutrons], taghit tag[maxNeutrons] ,
		double Ebeam , int nMult );
bool goodNeutronEvent(bandhit hits[maxNeutrons], int nMult, int& leadindex, int mcdataselect);
void getMcInfo( hipo::bank gen_particles , hipo::bank gen_info , genpart mcParts[maxGens] ,
		double &starttime, double &weight, double &Ebeam , int &genMult );
#endif
