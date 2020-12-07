#ifndef E_PID_HH
#define E_PID_HH

class clashit;
class e_pid{

public:

	e_pid();
	~e_pid();

	void setParamsRGB(double Ebeam);
	void setParamsRGA();	
	void setInterval(double newInterval);
	bool isElectron(clashit * eHit);
 	
private:
	
	bool SF_Epcal_Cut(int sector, double SF, double Epcal);
	bool SFpcal_SFecin_Cut(double p, double Epcal, double Eecin);
	double meanSF(int sector, double Epcal);
	double sigmaSF(int sector, double Epcal);
	double maxSF(int sector, double Epcal);
	double minSF(int sector, double Epcal);
	double FF(double E, double sf1, double sf2, double sf3, double sf4);

	double params[6][8];
	double interval;
};

#endif

