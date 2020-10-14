#ifndef CLASFID_HH
#define CLASFID_HH

class TString;
class TF1;

class clas12fiducial{

public:

	clas12fiducial();
	~clas12fiducial();

	int GetElectronAcceptance(double, double, double);


private:

	void ReadMomentumParameters(TString, TString);

	void InitMomentumFunctions();
	void InitFiducialFunctions();
	void SetFiducialParameters(double);

	double GetFiducialParameter(int, int, int, double);
	
	double momentumPar[2][6][3][3];	

	TF1* fiducialFunction[2][6];
	TF1* momentumFunction[2][6][3];

};

#endif

