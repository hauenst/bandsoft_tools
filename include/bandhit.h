#ifndef __BANDHIT_H__
#define __BANDHIT_H__

#include "TObject.h"
#include "TVector3.h"
#include "constants.h"

class bandhit : public TObject {
	public:
		bandhit();
		~bandhit();


		void Clear();
		
		// Some custom get functions
		TVector3 getDL		(void)		{return TVector3(X,Y,Z)	;}
		double getBeta		(void)		{return TVector3(X,Y,Z).Mag() / Tof / cAir; }
		TVector3 getMomentumN	(void){
			double beta = TVector3(X,Y,Z).Mag() / Tof / cAir;
			double mom = mN / sqrt(1./pow(beta,2) - 1.);
			TVector3 momN; momN.SetMagThetaPhi(mom,TVector3(X,Y,Z).Theta(),TVector3(X,Y,Z).Phi() );
			return momN;
		}


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
		void setStatus		(double iStatus		)		{Status		= iStatus	; return;}
                                                                                                                  
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
		double	getStatus	(void)		{return Status		;}
                                                                             
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

		ClassDef(bandhit,3);
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
		int	Status		;

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

#endif
