#include "taghit.h"

ClassImp(taghit);

taghit::taghit(){}	// Empty constructor
taghit::~taghit(){}	// Empty destructor

void taghit::Clear(){
	MomentumE.Clear();
	MomentumN.Clear();
	MomentumQ.Clear();
	MomentumB.Clear();

	PhiNQ		= 0;
	ThetaNQ		= 0;
	Wp		= 0;
	As		= 0;
	Pt.Clear();
	Xp		= 0;
	Xp_WP		= 0;
	Xp_Bj		= 0;
	Xp_PRC		= 0;
}
