#ifndef __TAGHIT_H__
#define __TAGHIT_H__

#include "TObject.h"
#include "TVector3.h"
#include "constants.h"

class taghit : public TObject {
	public:
		taghit();
		~taghit();

		void Clear();

		void setMomentumE	(TVector3	iMomentumE	)	{MomentumE	= iMomentumE	; return;}
		void setMomentumN	(TVector3	iMomentumN	)	{MomentumN	= iMomentumN	; return;}
		void setMomentumQ	(TVector3	iMomentumQ	)	{MomentumQ	= iMomentumQ	; return;}
		void setMomentumB	(TVector3	iMomentumB	)	{MomentumB	= iMomentumB	; return;}

		void setPhiNQ		(double		iPhiNQ		)	{PhiNQ		= iPhiNQ	; return;}
		void setThetaNQ		(double		iThetaNQ	)	{ThetaNQ	= iThetaNQ	; return;}
		void setWp		(double		iWp		)	{Wp		= iWp		; return;}
		void setAs		(double		iAs		)	{As		= iAs		; return;}
		void setPt		(TVector3	iPt		)	{Pt		= iPt		; return;}
		void setXp		(double		iXp		)	{Xp		= iXp		; return;}
		void setXp_WP		(double		iXp_WP		)	{Xp_WP		= iXp_WP	; return;}
		void setXp_Bj		(double		iXp_Bj		)	{Xp_Bj		= iXp_Bj	; return;}
		void setXp_PRC		(double		iXp_PRC		)	{Xp_PRC		= iXp_PRC	; return;}

		TVector3 getMomentumE	(void)	{return MomentumE;}
		TVector3 getMomentumN	(void)	{return MomentumN;}
		TVector3 getMomentumQ	(void)	{return MomentumQ;}
		TVector3 getMomentumB	(void)	{return MomentumB;}

		double getPhiNQ		(void)	{return PhiNQ	;}
		double getThetaNQ	(void)	{return ThetaNQ	;}
		double getWp		(void)	{return Wp	;}
		double getAs		(void)	{return As	;}
		TVector3 getPt		(void)	{return Pt	;}
		double getXp		(void)	{return Xp	;}
		double getXp_WP		(void)	{return Xp_WP	;}
		double getXp_Bj		(void)	{return Xp_Bj	;}
		double getXp_PRC	(void)	{return Xp_PRC	;}

		ClassDef(taghit,3);
	private:
		TVector3	MomentumE;
		TVector3	MomentumN;
		TVector3	MomentumQ;
		TVector3	MomentumB;

		double		PhiNQ;
		double		ThetaNQ;
		double		Wp;
		double		As;
		TVector3	Pt;
		double		Xp;
		double		Xp_WP;
		double		Xp_Bj;
		double		Xp_PRC;
};



#endif
