#ifndef __CLASHIT_H__
#define __CLASHIT_H__

#include "TObject.h"
#include "constants.h"

class clashit : public TObject {
	public:
		clashit();
		~clashit();


		void Clear();
		
		void setSector		(int	iSector		)		{Sector		= iSector	; return;}
		void setPID		(int	iPID		)		{PID		= iPID		; return;}
		void setCharge		(int	iCharge		)		{Charge		= iCharge	; return;}
		void setStatus		(int	iStatus		)		{Status		= iStatus	; return;}

		void setTime		(double	iTime		)		{Time		= iTime		; return;}
		void setBeta		(double	iBeta		)		{Beta		= iBeta		; return;}
		void setChi2		(double	iChi2		)		{Chi2		= iChi2		; return;}
		void setEtot		(double	iEtot		)		{Etot		= iEtot		; return;}
		void setEpcal		(double	iEpcal		)		{Epcal		= iEpcal	; return;}
		void setEoP		(double	iEoP		)		{EoP		= iEoP		; return;}
		void setTimeScint	(double	iTimeScint	)		{TimeScint	= iTimeScint	; return;}
		void setPathScint	(double	iPathScint	)		{PathScint	= iPathScint	; return;}
		void setU		(double	iU		)		{U		= iU		; return;}
		void setV		(double	iV		)		{V		= iV		; return;}
		void setW		(double	iW		)		{W		= iW		; return;}
		void setVtx		(double	iVtx		)		{Vtx		= iVtx		; return;}
		void setVty		(double	iVty		)		{Vty		= iVty		; return;}
		void setVtz		(double	iVtz		)		{Vtz		= iVtz		; return;}

		void setMomentum	(double	iMomentum	)		{Momentum	= iMomentum	; return;}
		void setTheta		(double	iTheta		)		{Theta		= iTheta	; return;}
		void setPhi		(double	iPhi		)		{Phi		= iPhi		; return;}

		void setQ		(double	iQ		)		{Q		= iQ		; return;}
		void setThetaQ		(double	iThetaQ		)		{ThetaQ		= iThetaQ	; return;}
		void setPhiQ		(double	iPhiQ		)		{PhiQ		= iPhiQ		; return;}

		void setQ2		(double	iQ2		)		{Q2		= iQ2		; return;}
		void setOmega		(double	iOmega		)		{Omega		= iOmega	; return;}
		void setXb		(double	iXb		)		{Xb		= iXb		; return;}
		void setW2		(double	iW2		)		{W2		= iW2		; return;}


		int	getSector	(void)		{return Sector		;}
		int	getPID		(void)		{return PID		;}
		int	getCharge	(void)		{return Charge		;}
		int	getStatus	(void)		{return Status		;}
                                                                
		double	getTime		(void)		{return Time		;}
		double	getBeta		(void)		{return Beta		;}
		double	getChi2		(void)		{return Chi2		;}
		double	getEtot		(void)		{return Etot		;}
		double	getEpcal	(void)		{return Epcal		;}
		double	getEoP		(void)		{return EoP		;}
		double	getTimeScint	(void)		{return TimeScint	;}
		double	getPathScint	(void)		{return PathScint	;}
		double	getU		(void)		{return U		;}
		double	getV		(void)		{return V		;}
		double	getW		(void)		{return W		;}
		double	getVtx		(void)		{return Vtx		;}
		double	getVty		(void)		{return Vty		;}
		double	getVtz		(void)		{return Vtz		;}
                                                                
		double	getMomentum	(void)		{return Momentum	;}
		double	getTheta	(void)		{return Theta		;}
		double	getPhi		(void)		{return Phi		;}
                                                               
		double	getQ		(void)		{return Q		;}
		double	getThetaQ	(void)		{return ThetaQ		;}
		double	getPhiQ		(void)		{return PhiQ		;}
                                                              
		double	getQ2		(void)		{return Q2		;}
		double	getOmega	(void)		{return Omega		;}
		double	getXb		(void)		{return Xb		;}
		double	getW2		(void)		{return W2		;}

		ClassDef(clashit,2);
	private:
		int	Sector		;
		int 	PID		;
		int	Charge		;
		int	Status		;

		double	Time		;
		double	Beta		;
		double	Chi2		;
		double	Etot		;
		double	Epcal		;
		double	EoP		;
		double	TimeScint	;
		double	PathScint	;
		double 	U		;
		double	V		;
		double  W		;
		double	Vtx		;
		double 	Vty		;
		double 	Vtz		;

		double 	Momentum	;
		double	Theta		;
		double 	Phi		;
		
		double 	Q		;
		double	ThetaQ		;
		double 	PhiQ		;

		double 	Q2		;
		double 	Omega		;
		double	Xb		;
		double	W2		;
};



#endif
