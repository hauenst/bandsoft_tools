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
#include "particles.h"
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
//const double time_thresBANDhit = 300;


int getRunNumber( string filename );
void getEventInfo( hipo::bank eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, BCalorimeter calorimeter, hipo::bank scintillator, hipo::bank DC_Track, hipo::bank DC_Traj, hipo::bank cherenkov,
			int pbankIndex,
			clashit &electron,
			double starttime , int thisRun , double Ebeam );
void getTaggedInfo( clashit eHit, bandhit nHit[maxNeutrons], taghit tag[maxNeutrons] ,
		double Ebeam , int nMult );
bool goodNeutronEvent(bandhit hits[maxNeutrons], int nMult, int& leadindex, int mcdataselect, int& nPass );
void getMcInfo( hipo::bank gen_particles , hipo::bank gen_info , genpart mcParts[maxGens] ,
		double &starttime, double &weight, double &Ebeam , int &genMult );
void getParticleInfo( hipo::bank claspart, particles part[maxParticles], hipo::bank scintillator ,int& multiplicity );
//Parametrization from Giovanni (GWU) and FX based on double pion analysis in RGA
void smearRGA(TVector3 &vpar);
void recalculate_clashit_kinematics(clashit &input_ehit, double Ebeam, TVector3 &smeared_electron);

#endif
