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

using namespace std;

const int maxProtons	= 100;
const int maxNeutrons	= 15;

// Need to add theta, phi, etc..

class bandhit {
	public:
		void Clear();
		
		void setSector(int)			
		void setLayer(int)			
		void setComponent(int);
		void setBarID(int);
		void setEdep(double);
		void setTof(double);
		void setTofFadc(double);
		
		int getSector		(void)		{return	sector		;}
		int getLayer		(void)		{return layer		;}
		int getComponent	(void)		{return component	;}
		double getBarID		(void)		{return barID		;}
		double getEdep		(void)		{return edep		;}
		double getTof		(void)		{return tof		;}
		double getTofFadc	(void)		{return tof_fadc	;}
	private:
		int sector;
		int layer;
		int component;
		int barID;
	
		double edep;
		double tof;
		double tof_fadc;

		double dL;
		
		double adcL;
		double adcR;
		double ampL;
		double ampR;
	
		double timeL;
		double timeR;
		double timecorrL;
		double timecorrR;
		double time_fadcL;
		double time_fadcR;
};


int getRunNumber( string filename );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex, 
			double& time, int& charge, double& beta, double& chi2pid, int& status );
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getProtonInfo( BParticle particles, double pid[maxProtons], TVector3 momentum[maxProtons], TVector3 vertex[maxProtons],
			double time[maxProtons], double charge[maxProtons], double beta[maxProtons], double chi2pid[maxProtons], double status[maxProtons] , int& multiplicity );
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult );
void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, int& mult, double id[maxNeutrons], double edep[maxNeutrons],
			double time[maxNeutrons], double timefadc[maxNeutrons], TVector3 path[maxNeutrons] , double starttime , int thisRun);
bool pointsToBand(double theta,double phi,double z_m);

void LoadGlobalShift( string filename_tdc , string filename_fadc );
void LoadRunByRunShift();
#endif
