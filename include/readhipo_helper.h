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

using namespace std;

const int maxProtons	= 100;
const int maxNeutrons	= 15;

int getRunNumber( string filename );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex, 
			double& time, int& charge, double& beta, double& chi2pid, int& status );
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getProtonInfo( BParticle particles, double pid[maxProtons], TVector3 momentum[maxProtons], TVector3 vertex[maxProtons],
			double time[maxProtons], double charge[maxProtons], double beta[maxProtons], double chi2pid[maxProtons], double status[maxProtons] , int& multiplicity );
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult );
void getNeutronInfo( BBand band_hits, int& mult, double id[maxNeutrons], double edep[maxNeutrons],
			double time[maxNeutrons], TVector3 path[maxNeutrons] , double starttime , int thisRun);
bool pointsToBand(double theta,double phi,double z_m);

void LoadGlobalShift( string filename_tdc , string filename_fadc );
void LoadRunByRunShift();
#endif
