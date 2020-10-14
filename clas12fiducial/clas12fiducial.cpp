#include "clas12fiducial.h"

#include "TF1.h"
#include "TString.h"
#include "TMath.h"

#include <fstream>
#include <unistd.h>

using namespace std;

clas12fiducial::clas12fiducial() {

	ReadMomentumParameters("upper_momentum_fit.dat", "lower_momentum_fit.dat");
	InitMomentumFunctions();
	InitFiducialFunctions();

}

clas12fiducial::~clas12fiducial() {
}

void clas12fiducial::ReadMomentumParameters(TString upperParFile, TString lowerParFile) {

	ifstream upperFile(upperParFile);
	ifstream lowerFile(lowerParFile);
	
	int sector, fidPar;
	double momPar0, momPar1, momPar2;

	while(upperFile >> sector >> fidPar >> momPar0 >> momPar1 >> momPar2) {

		momentumPar[0][sector-1][fidPar][0] = momPar0;
		momentumPar[0][sector-1][fidPar][1] = momPar1;
		momentumPar[0][sector-1][fidPar][2] = momPar2;

	}

	while(lowerFile >> sector >> fidPar >> momPar0 >> momPar1 >> momPar2) {

		momentumPar[1][sector-1][fidPar][0] = momPar0;
		momentumPar[1][sector-1][fidPar][1] = momPar1;
		momentumPar[1][sector-1][fidPar][2] = momPar2;

	}


}

void clas12fiducial::InitMomentumFunctions() {

	TString momentum_functions[3] = {"[0] + exp([1]*x+[2])", "[0] - exp([1]*x+[2])", "[0] + exp([1]*x+[2])"};
	TString boundary[2] = {"upper","lower"};

	for(int sector = 0; sector < 6; sector++) {
		for(int fidPar = 0; fidPar < 3; fidPar++) {
			for(int i = 0; i < 2; i++) {
				momentumFunction[i][sector][fidPar] = new TF1(Form("%sMomFunc_%i_%i", boundary[i].Data(), sector, fidPar), momentum_functions[fidPar], 0., 12.);
				for(int momPar = 0; momPar < 3; momPar++) {
					momentumFunction[i][sector][fidPar]->SetParameter(momPar, momentumPar[i][sector][fidPar][momPar]);
				}	
			}
		}
	}
}


void clas12fiducial::InitFiducialFunctions() {

	TString fiducial_functions[2] = {"[0] - exp([1]*x+[2])", "[0] + exp([1]*x+[2])"};
	TString boundary[2] = {"upper","lower"};

	for(int sector = 0; sector < 6; sector++) {
		for(int i = 0; i < 2; i++) {
			fiducialFunction[i][sector] = new TF1(Form("%sFidFunc_%i", boundary[i].Data(), sector), fiducial_functions[i], 0., 60.);
		}
	}
}

void clas12fiducial::SetFiducialParameters(double pe) {

	for(int sector = 0; sector < 6; sector++) {
		for(int i = 0; i < 2; i++) {
			for(int fidPar = 0; fidPar < 3; fidPar++) {
				double par = GetFiducialParameter(i, sector, fidPar, pe);
				fiducialFunction[i][sector]->SetParameter(fidPar, par);
			}
		}
	}
}


double clas12fiducial::GetFiducialParameter(int boundary, int sector, int par, double pe) {

	return momentumFunction[boundary][sector][par]->Eval(pe);

}


int clas12fiducial::GetElectronAcceptance(double theta, double phi, double pe) {

	theta*=TMath::RadToDeg();
	phi*=TMath::RadToDeg();

	SetFiducialParameters(pe);

	for(int sector = 0; sector < 6; sector++) {

		double maxPhi = fiducialFunction[0][sector]->Eval(theta);
		double minPhi = fiducialFunction[1][sector]->Eval(theta);

		if (phi > minPhi && phi < maxPhi) {
			return (sector + 1);
		}
	}

	// if it's not already in a sector, 
	// check for values of (phi - 360) 
/*
	 for(int sector = 0; sector < 6; sector++) {

                double maxPhi = fiducialFunction[0][sector]->Eval(theta);
                double minPhi = fiducialFunction[1][sector]->Eval(theta);

                if ((phi - 360.) > minPhi && (phi - 360.) < maxPhi) {
                        return (sector + 1);
                }
        }
*/

	return -1;
}


