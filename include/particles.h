#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include "TObject.h"
#include "constants.h"

class particles : public TObject {
	public:
		particles();
		~particles();

		void Clear();
		void Print();

		void setPID		(int	iPID		)		{PID		= iPID		; return;}
		void setCharge		(int	iCharge		)		{Charge		= iCharge	; return;}
		void setStatus		(int	iStatus		)		{Status		= iStatus	; return;}

		void setTime		(double	iTime		)		{Time		= iTime		; return;}
		void setBeta		(double	iBeta		)		{Beta		= iBeta		; return;}
		void setChi2		(double	iChi2		)		{Chi2		= iChi2		; return;}
		void setPindex		(double iPindex		)		{Pindex		= iPindex	; return;}

		void setVtx		(double	iVtx		)		{Vtx		= iVtx		; return;}
		void setVty		(double	iVty		)		{Vty		= iVty		; return;}
		void setVtz		(double	iVtz		)		{Vtz		= iVtz		; return;}

		void setMomentum	(double	iMomentum	)		{Momentum	= iMomentum	; return;}
		void setTheta		(double	iTheta		)		{Theta		= iTheta	; return;}
		void setPhi		(double	iPhi		)		{Phi		= iPhi		; return;}

		// Scintillator banks
		void setScint_detector	(std::vector<int>    iScint_detector		)	{Scint_detector		= iScint_detector		; return; }
		void setScint_status	(std::vector<int>    iScint_status		)	{Scint_status		= iScint_status			; return; }
		void setScint_sector	(std::vector<int>    iScint_sector		)	{Scint_sector		= iScint_sector			; return; }
		void setScint_layer	(std::vector<int>    iScint_layer		)	{Scint_layer		= iScint_layer			; return; }
		void setScint_component	(std::vector<int>    iScint_component		)	{Scint_component	= iScint_component		; return; }
		void setScint_Edep	(std::vector<double> iScint_Edep		)	{Scint_Edep		= iScint_Edep			; return; }
		void setScint_time	(std::vector<double> iScint_time		)	{Scint_time		= iScint_time			; return; }
		void setScint_path	(std::vector<double> iScint_path		)	{Scint_path		= iScint_path			; return; }
		void setScint_chi2	(std::vector<double> iScint_chi2		)	{Scint_chi2		= iScint_chi2			; return; }
		void setScint_x		(std::vector<double> iScint_x			)	{Scint_x		= iScint_x			; return; }
		void setScint_y		(std::vector<double> iScint_y			)	{Scint_y		= iScint_y			; return; }
		void setScint_z		(std::vector<double> iScint_z			)	{Scint_z		= iScint_z			; return; }

		int	getPID		(void)		{return PID		;}
		int	getCharge	(void)		{return Charge		;}
		int	getStatus	(void)		{return Status		;}

		double	getTime		(void)		{return Time		;}
		double	getBeta		(void)		{return Beta		;}
		double	getChi2		(void)		{return Chi2		;}
		double 	getPindex	(void)		{return Pindex		;}
		double	getVtx		(void)		{return Vtx		;}
		double	getVty		(void)		{return Vty		;}
		double	getVtz		(void)		{return Vtz		;}

		double	getMomentum	(void)		{return Momentum	;}
		double	getTheta	(void)		{return Theta		;}
		double	getPhi		(void)		{return Phi		;}

		std::vector<int>    getScint_detector		(void)		{return Scint_detector	;}
		std::vector<int>    getScint_status		(void)		{return Scint_status	;}
		std::vector<int>    getScint_sector		(void)		{return Scint_sector	;}
		std::vector<int>    getScint_layer		(void)		{return Scint_layer	;}
		std::vector<int>    getScint_component		(void)		{return Scint_component	;}
		std::vector<double> getScint_Edep		(void)		{return Scint_Edep	;}
		std::vector<double> getScint_time		(void)		{return Scint_time	;}
		std::vector<double> getScint_path		(void)		{return Scint_path	;}
		std::vector<double> getScint_chi2		(void)		{return Scint_chi2	;}
		std::vector<double> getScint_x			(void)		{return Scint_x		;}
		std::vector<double> getScint_y			(void)		{return Scint_y		;}
		std::vector<double> getScint_z			(void)		{return Scint_z		;}

		ClassDef(particles,1);

	private:
		int 	PID		;
		int	Charge		;
		int	Status		;

		double	Time		;
		double	Beta		;
		double	Chi2		;
		int     Pindex		;
		double	Vtx		;
		double 	Vty		;
		double 	Vtz		;

		double 	Momentum	;
		double	Theta		;
		double 	Phi		;

		std::vector<int>	Scint_detector	;
		std::vector<int>	Scint_status	;
		std::vector<int>	Scint_sector	;
		std::vector<int>	Scint_layer	;
		std::vector<int>	Scint_component	;
		std::vector<double>	Scint_Edep	;
		std::vector<double>	Scint_time	;
		std::vector<double>	Scint_path	;
		std::vector<double>	Scint_chi2	;
		std::vector<double>	Scint_x		;
		std::vector<double>	Scint_y		;
		std::vector<double>	Scint_z		;
};



#endif
