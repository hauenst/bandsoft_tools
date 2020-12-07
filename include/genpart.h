#ifndef __GENPART_H__
#define __GENPART_H__

#include "TObject.h"
#include "constants.h"
#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

class genpart : public TObject {
	public:
		genpart();
		~genpart();


		void Clear();
		
		void setPID		(int	iPID		)		{PID		= iPID		; return;}

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


		int	getPID		(void)		{return PID		;}

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

		ClassDef(genpart,1);
	private:
		int 	PID		;

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

		clashit	electron	;
		bandhit	neutron		;
		taghit	tag		;
};



#endif
