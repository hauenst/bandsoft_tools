#include "genpart.h"

ClassImp(genpart);

genpart::genpart(){}	// Empty constructor
genpart::~genpart(){}	// Empty destructor

void genpart::Clear(){
	PID		= 0;

	Vtx		= 0;
	Vty		= 0;
	Vtz		= 0;

	Momentum	= 0;
	Theta		= 0;
	Phi		= 0;

	Q		= 0;
	ThetaQ		= 0;
	PhiQ		= 0;

	Q2		= 0;
	Omega		= 0;
	Xb		= 0;
	W2		= 0;
}
