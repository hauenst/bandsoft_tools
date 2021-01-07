#ifndef E_PID_HH
#define E_PID_HH

#include "TF1.h"
#include "TCanvas.h"

class clashit;
class e_pid{

public:

	e_pid();
	~e_pid();

	bool isElectron(clashit * eHit);
	bool isElectronLoose(clashit * eHit);

	void setParamsRGB(double Ebeam);
	void fillParams();
	void setIntervalEpcal(double newInterval);
	void setIntervalMom(double newInterval);
	void setColor(int newColor);

	void drawEpcal(int sector, TCanvas * myCanvas);
	void drawMom(int sector, TCanvas * myCanvas);

private:
	
	bool SF_Epcal_Cut(int sector, double Epcal, double SF);
	bool SF_Mom_Cut(int sector, double p, double SF);
	bool SFpcal_SFecin_Cut(double p, double Epcal, double Eecin);

	bool SF_Cut(int sector, double x, double SF, double params[6][6], double interval);
	double meanSF(int sector, double x, double params[6][6]);
	double sigmaSF(int sector, double x, double params[6][6]);
	double maxSF(int sector, double x, double params[6][6], double interval);
	double minSF(int sector, double x, double params[6][6], double interval);
	double FF(double x, double sf1, double sf2, double sf3);

	char paramFileNameEpcal[100];
	double paramsEpcal[6][6];
	double intervalEpcal;
	char paramFileNameMom[100];
	double paramsMom[6][6];
	double intervalMom;
	int color;
};

#endif

