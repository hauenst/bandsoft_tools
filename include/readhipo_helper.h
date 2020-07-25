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

using namespace std;

const int maxProtons	= 100;
const int maxNeutrons	= 200;


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
void getElectronInfo( BParticle particles, BCalorimeter calorimeter, BScintillator scintillator,
			clashit &electron,
			double starttime , int thisRun , double Ebeam );
void getTaggedInfo( clashit eHit, bandhit nHit[maxNeutrons], taghit tag[maxNeutrons] ,
		double Ebeam , int nMult );
#endif
