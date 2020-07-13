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
		bandhit();
		~bandhit();

		void Clear();
		void LoadShifts(){return;} // EMPTY FOR NOW -- TODO load shifts and calculate momentum
		
		void setSector		(int	iSector		)		{Sector		= iSector	; return;}
		void setLayer		(int	iLayer		)		{Layer		= iLayer	; return;}
		void setComponent	(int	iComponent	)		{Component	= iComponent	; return;}
		void setBarID		(int	iBarID		)		{BarID		= iBarID	; return;}

		void setEdep		(double iEdep		)		{Edep		= iEdep		; return;}
		void setTof		(double iTof		)		{Tof		= iTof		; return;}
		void setTofFadc		(double iTofFadc	)		{TofFadc	= iTofFadc	; return;}
		void setTdiff		(double iTdiff		)		{Tdiff		= iTdiff	; return;}
		void setTdiffFadc	(double iTdiffFadc	)		{TdiffFadc	= iTdiffFadc	; return;}
		void setX		(double iX		)		{X		= iX		; return;}
		void setY		(double iY		)		{Y		= iY		; return;}
		void setZ		(double iZ		)		{Z		= iZ		; return;}
                                                                                                                  
		void setRawLtdc		(double iRawLtdc	)		{RawLtdc	= iRawLtdc	; return;}
		void setRawRtdc		(double iRawRtdc	)		{RawRtdc	= iRawRtdc	; return;}
		void setRawLtdccorr	(double iRawLtdccorr	)		{RawLtdccorr	= iRawLtdccorr	; return;}
		void setRawRtdccorr	(double iRawRtdccorr	)		{RawRtdccorr	= iRawRtdccorr	; return;}
		void setRawLtfadc	(double iRawLtfadc	)		{RawLtfadc	= iRawLtfadc	; return;}
		void setRawRtfadc	(double iRawRtfadc	)		{RawRtfadc	= iRawRtfadc	; return;}
		void setRawLamp		(double iRawLamp	)		{RawLamp	= iRawLamp	; return;}
		void setRawRamp		(double iRawRamp	)		{RawRamp	= iRawRamp	; return;}
		void setRawLadc		(double iRawLadc	)		{RawLadc	= iRawLadc	; return;}
		void setRawRadc		(double iRawRadc	)		{RawRadc	= iRawRadc	; return;}
                                                                                                                  
		void setPmtLtdc		(double iPmtLtdc	)		{PmtLtdc	= iPmtLtdc	; return;}
		void setPmtRtdc		(double iPmtRtdc	)		{PmtRtdc	= iPmtRtdc	; return;}
		void setPmtLtfadc	(double iPmtLtfadc	)		{PmtLtfadc	= iPmtLtfadc	; return;}
		void setPmtRtfadc	(double iPmtRtfadc	)		{PmtRtfadc	= iPmtRtfadc	; return;}
		void setPmtLamp		(double iPmtLamp	)		{PmtLamp	= iPmtLamp	; return;}
		void setPmtRamp		(double iPmtRamp	)		{PmtRamp	= iPmtRamp	; return;}
		void setPmtLadc		(double iPmtLadc	)		{PmtLadc	= iPmtLadc	; return;}
		void setPmtRadc		(double iPmtRadc	)		{PmtRadc	= iPmtRadc	; return;}
		void setPmtLped		(double iPmtLped	)		{PmtLped	= iPmtLped	; return;}
		void setPmtRped		(double iPmtRped	)		{PmtRped	= iPmtRped	; return;}

	
		int	getSector	(void)		{return	Sector		;}
		int 	getLayer	(void)		{return Layer		;}
		int 	getComponent	(void)		{return Component	;}
		double 	getBarID	(void)		{return BarID		;}
		double 	getEdep		(void)		{return Edep		;}
		double 	getTof		(void)		{return Tof		;}
		double 	getTofFadc	(void)		{return TofFadc		;}
		double	getTdiff	(void)		{return Tdiff		;}
		double	getTdiffFadc	(void)		{return TdiffFadc	;}
		double	getX		(void)		{return X		;}
		double	getY		(void)		{return Y		;}
		double	getZ		(void)		{return Z		;}
                                                                             
		double 	getRawLtdc	(void)		{return	RawLtdc		;}	
		double 	getRawRtdc	(void)		{return	RawRtdc		;}	
		double 	getRawLtdccorr	(void)		{return	RawLtdccorr	;}
		double 	getRawRtdccorr	(void)		{return	RawRtdccorr	;}
		double 	getRawLtfadc	(void)		{return	RawLtfadc	;}
		double 	getRawRtfadc	(void)		{return	RawRtfadc	;}
		double 	getRawLamp	(void)		{return	RawLamp		;}	
		double 	getRawRamp	(void)		{return	RawRamp		;}	
		double 	getRawLadc	(void)		{return	RawLadc		;}	
		double 	getRawRadc	(void)		{return	RawRadc		;}	
                                                                             
		double 	getPmtLtdc	(void)		{return	PmtLtdc		;}	
		double 	getPmtRtdc	(void)		{return	PmtRtdc		;}	
		double 	getPmtLtfadc	(void)		{return	PmtLtfadc	;}
		double 	getPmtRtfadc	(void)		{return	PmtRtfadc	;}
		double 	getPmtLamp	(void)		{return	PmtLamp		;}	
		double 	getPmtRamp	(void)		{return	PmtRamp		;}	
		double 	getPmtLadc	(void)		{return	PmtLadc		;}	
		double 	getPmtRadc	(void)		{return	PmtRadc		;}	
		double 	getPmtLped	(void)		{return	PmtLped		;}
		double 	getPmtRped	(void)		{return PmtRped		;}

		ClassDef(bandhit,1);
	private:
		int	Sector		;
		int 	Layer		;
		int 	Component	;
		double 	BarID		;
		double 	Edep		;
		double 	Tof		;
		double 	TofFadc		;
		double 	Tdiff		;
		double 	TdiffFadc	;
		double 	X		;
		double 	Y		;
		double 	Z		;

		double 	RawLtdc		;
		double 	RawRtdc		;
		double 	RawLtdccorr	;
		double 	RawRtdccorr	;
		double 	RawLtfadc	;
		double 	RawRtfadc	;
		double 	RawLamp		;
		double 	RawRamp		;
		double 	RawLadc		;
		double 	RawRadc		;
			     	
		double 	PmtLtdc		;
		double 	PmtRtdc		;
		double 	PmtLtfadc	;
		double 	PmtRtfadc	;
		double 	PmtLamp		;
		double 	PmtRamp		;
		double 	PmtLadc		;
		double 	PmtRadc		;
		double 	PmtLped		;
		double 	PmtRped		;
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

void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun);
bool pointsToBand(double theta,double phi,double z_m);

void LoadGlobalShift( string filename_tdc , string filename_fadc );
void LoadRunByRunShift();
#endif
